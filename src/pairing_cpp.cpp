// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <utility>
#include <vector>
using namespace Rcpp;
using namespace RcppParallel;

// CPU compute kernels for run_pairing(..., backend = "cpp"). For one chunk of TCRalpha
// chains, these compute the same pairs as the GPU kernels in
// inst/python/pairing_new/pairing_core.py (MadHype.chunk() and TShell.chunk()), in
// parallel with RcppParallel. All other steps are done in R (R/utils_run_pairing.R).
//
// The GPU kernels compare every alpha chain with every beta chain using dense matrix
// products. Here, the sparsity of the data is used instead: each chain is only in a few
// wells, so for each alpha chain we only visit the beta chains in the same wells. Pairs
// with no shared wells are handled separately (see each kernel).
//
// Inputs are sparse matrices (Matrix::dgCMatrix) without explicit zeros:
//   alpha_t: wells x alpha chains, i.e. the wells (and read fractions) of each alpha chain
//   beta:    beta chains x wells,  i.e. the beta chains (and read fractions) in each well
// Alpha chains [start, end) (0-based) are processed. Returned "i" is the position in the
// chunk (0-based) and "j" is the beta chain (0-based), as for the GPU kernels.

namespace {

// read-only view of the slots of a dgCMatrix
struct SparseView {
  const RVector<int> p;
  const RVector<int> i;
  const RVector<double> x;
  SparseView(const IntegerVector& p, const IntegerVector& i, const NumericVector& x) : p(p), i(i), x(x) {}
};

// Count the wells shared by alpha chain a and every beta chain (cnt), and optionally sum the
// products of their read fractions (prod). Beta chains with at least one shared well are added
// to 'touched'. cnt (and prod) must be all zeros for the beta chains before the call.
inline void accumulate_shared_wells(const SparseView& alpha_t, const SparseView& beta, int a,
                                    std::vector<int>& cnt, std::vector<double>* prod,
                                    std::vector<int>& touched) {
  touched.clear();
  for (int pa = alpha_t.p[a]; pa < alpha_t.p[a + 1]; ++pa) {
    const int w = alpha_t.i[pa];
    const double xa = alpha_t.x[pa];
    for (int pb = beta.p[w]; pb < beta.p[w + 1]; ++pb) {
      const int j = beta.i[pb];
      if (cnt[j]++ == 0) touched.push_back(j);
      if (prod) (*prod)[j] += xa * beta.x[pb];
    }
  }
}

// ---------------------------------------------------------------------------------------
// MAD-HYPE: an alpha/beta pair with wij wells with both chains, wi = wa - wij wells with only
// the alpha chain, and wj = wb - wij wells with only the beta chain is called if
// wj < mdh(wij, wi). For pairs with no shared wells (wij = 0), the pair is called if
// wb < mdh(0, wa), so these are found from the beta chains sorted by wb (beta_by_wb).
struct MadhypeWorker : public Worker {
  const SparseView& alpha_t;
  const SparseView& beta;
  const RVector<int> wa, wb, beta_by_wb;
  const RMatrix<int> mdh;
  const int start, n_beta;
  std::vector<std::vector<std::pair<int, int> > >& pairs; // per alpha chain: (beta, wij)

  MadhypeWorker(const SparseView& alpha_t, const SparseView& beta, const IntegerVector& wa,
                const IntegerVector& wb, const IntegerVector& beta_by_wb, const IntegerMatrix& mdh,
                int start, int n_beta, std::vector<std::vector<std::pair<int, int> > >& pairs)
    : alpha_t(alpha_t), beta(beta), wa(wa), wb(wb), beta_by_wb(beta_by_wb), mdh(mdh),
      start(start), n_beta(n_beta), pairs(pairs) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<int> cnt(n_beta, 0);
    std::vector<int> touched;
    for (std::size_t k = begin; k < end; ++k) {
      const int a = start + (int)k;
      accumulate_shared_wells(alpha_t, beta, a, cnt, NULL, touched);
      std::vector<std::pair<int, int> >& out = pairs[k];
      for (std::size_t t = 0; t < touched.size(); ++t) {
        const int j = touched[t];
        const int wij = cnt[j];
        if (wb[j] - wij < mdh(wij, wa[a] - wij)) out.push_back(std::make_pair(j, wij));
      }
      const int thres_no_overlap = mdh(0, wa[a]);
      for (std::size_t q = 0; q < beta_by_wb.length(); ++q) {
        const int j = beta_by_wb[q];
        if (wb[j] >= thres_no_overlap) break;
        if (cnt[j] == 0) out.push_back(std::make_pair(j, 0));
      }
      std::sort(out.begin(), out.end()); // by beta chain, as for the GPU
      for (std::size_t t = 0; t < touched.size(); ++t) cnt[touched[t]] = 0;
    }
  }
};

// ---------------------------------------------------------------------------------------
// T-SHELL: for each alpha chain, the 3 beta chains with the highest Pearson correlation of
// read fractions across all n_wells wells. With prod = sum of products of read fractions over
// shared wells, mean = mean read fraction, and norm = L2 norm of (read fraction - mean):
//   r = (prod - n_wells * mean_a * mean_b) / (norm_a * norm_b)
// For pairs with no shared wells, prod = 0, so r is highest for the beta chains with the lowest
// mean_b / norm_b (beta_order). Ties are broken by choosing the beta chain with the higher
// index, which gives the same result as the (stable) argsort used on the GPU.
// If filter_before_top3, r = -1 for pairs with wij <= 2 or loss fraction sum >= 0.5 (computed
// in float32, as on the GPU), so pairs with no shared wells all have r = -1.

struct Candidate {
  double r;
  int j;
};
// is candidate a ranked above candidate b?
inline bool ranks_above(const Candidate& a, const Candidate& b) {
  return a.r > b.r || (a.r == b.r && a.j > b.j);
}
// keep the best n_top candidates in 'top' (best first)
inline void add_candidate(std::vector<Candidate>& top, std::size_t n_top, const Candidate& c) {
  if (top.size() == n_top && !ranks_above(c, top.back())) return;
  std::vector<Candidate>::iterator pos = top.begin();
  while (pos != top.end() && !ranks_above(c, *pos)) ++pos;
  top.insert(pos, c);
  if (top.size() > n_top) top.pop_back();
}

struct TshellWorker : public Worker {
  const SparseView& alpha_t;
  const SparseView& beta;
  const RVector<double> mean_a, norm_a, mean_b, norm_b;
  const RVector<int> wa, wb, beta_order;
  const bool filter_before_top3;
  const int start, n_beta, n_wells, n_top;
  RVector<int> out_i, out_j;
  RVector<double> out_r;

  TshellWorker(const SparseView& alpha_t, const SparseView& beta, const NumericVector& mean_a,
               const NumericVector& norm_a, const NumericVector& mean_b, const NumericVector& norm_b,
               const IntegerVector& wa, const IntegerVector& wb, const IntegerVector& beta_order,
               bool filter_before_top3, int start, int n_beta, int n_wells, int n_top,
               IntegerVector& out_i, IntegerVector& out_j, NumericVector& out_r)
    : alpha_t(alpha_t), beta(beta), mean_a(mean_a), norm_a(norm_a), mean_b(mean_b), norm_b(norm_b),
      wa(wa), wb(wb), beta_order(beta_order), filter_before_top3(filter_before_top3), start(start),
      n_beta(n_beta), n_wells(n_wells), n_top(n_top), out_i(out_i), out_j(out_j), out_r(out_r) {}

  inline bool filtered_out(int a, int j, int wij) const {
    if (wij <= 2) return true;
    const float fwij = (float)wij, fwa = (float)wa[a], fwb = (float)wb[j];
    const float loss_a_frac = (fwb - fwij) / (fwij + (fwb - fwij) + (fwa - fwij));
    const float loss_b_frac = (fwa - fwij) / (fwij + (fwb - fwij) + (fwa - fwij));
    return (loss_a_frac + loss_b_frac) >= 0.5f;
  }

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<int> cnt(n_beta, 0);
    std::vector<double> prod(n_beta, 0.0);
    std::vector<int> touched;
    std::vector<Candidate> top;
    top.reserve(n_top + 1);
    for (std::size_t k = begin; k < end; ++k) {
      const int a = start + (int)k;
      accumulate_shared_wells(alpha_t, beta, a, cnt, &prod, touched);
      top.clear();
      // pairs with shared wells
      for (std::size_t t = 0; t < touched.size(); ++t) {
        const int j = touched[t];
        Candidate c;
        c.j = j;
        if (filter_before_top3 && filtered_out(a, j, cnt[j])) {
          c.r = -1.0;
        } else {
          c.r = (prod[j] - n_wells * mean_a[a] * mean_b[j]) / (norm_a[a] * norm_b[j]);
        }
        add_candidate(top, n_top, c);
      }
      // pairs with no shared wells: only the best n_top of these can be in the top n_top
      int n_added = 0;
      if (filter_before_top3) {
        for (int j = n_beta - 1; j >= 0 && n_added < n_top; --j) {
          if (cnt[j] > 0) continue;
          Candidate c = {-1.0, j};
          add_candidate(top, n_top, c);
          ++n_added;
        }
      } else {
        for (int q = 0; q < n_beta && n_added < n_top; ++q) {
          const int j = beta_order[q];
          if (cnt[j] > 0) continue;
          Candidate c = {(0.0 - n_wells * mean_a[a] * mean_b[j]) / (norm_a[a] * norm_b[j]), j};
          add_candidate(top, n_top, c);
          ++n_added;
        }
      }
      // write out from lowest to highest correlation, as for the GPU
      for (int t = 0; t < n_top; ++t) {
        const Candidate& c = top[n_top - 1 - t];
        const std::size_t pos = k * n_top + t;
        out_i[pos] = (int)k;
        out_j[pos] = c.j;
        out_r[pos] = c.r;
      }
      for (std::size_t t = 0; t < touched.size(); ++t) {
        cnt[touched[t]] = 0;
        prod[touched[t]] = 0.0;
      }
    }
  }
};

void check_chunk(int start, int end, int n_alpha) {
  if (start < 0 || end > n_alpha || start >= end) stop("Invalid chunk of alpha chains.");
}

} // namespace

// [[Rcpp::export]]
List madhype_chunk_cpp(S4 alpha_t, S4 beta, IntegerVector wa, IntegerVector wb, IntegerMatrix mdh,
                       IntegerVector beta_by_wb, int start, int end) {
  IntegerVector a_p = alpha_t.slot("p"), a_i = alpha_t.slot("i"), b_p = beta.slot("p"), b_i = beta.slot("i");
  NumericVector a_x = alpha_t.slot("x"), b_x = beta.slot("x");
  IntegerVector b_dim = beta.slot("Dim");
  const int n_alpha = a_p.length() - 1, n_beta = b_dim[0];
  check_chunk(start, end, n_alpha);
  SparseView alpha_view(a_p, a_i, a_x), beta_view(b_p, b_i, b_x);

  const int n_chunk = end - start;
  std::vector<std::vector<std::pair<int, int> > > pairs(n_chunk);
  MadhypeWorker worker(alpha_view, beta_view, wa, wb, beta_by_wb, mdh, start, n_beta, pairs);
  parallelFor(0, n_chunk, worker, 16);

  std::size_t n_pairs = 0;
  for (int k = 0; k < n_chunk; ++k) n_pairs += pairs[k].size();
  IntegerVector out_i(n_pairs), out_j(n_pairs), out_wij(n_pairs);
  std::size_t pos = 0;
  for (int k = 0; k < n_chunk; ++k) {
    for (std::size_t t = 0; t < pairs[k].size(); ++t, ++pos) {
      out_i[pos] = k;
      out_j[pos] = pairs[k][t].first;
      out_wij[pos] = pairs[k][t].second;
    }
  }
  return List::create(_["i"] = out_i, _["j"] = out_j, _["wij"] = out_wij);
}

// [[Rcpp::export]]
List tshell_chunk_cpp(S4 alpha_t, S4 beta, NumericVector mean_a, NumericVector norm_a,
                      NumericVector mean_b, NumericVector norm_b, IntegerVector wa, IntegerVector wb,
                      IntegerVector beta_order, bool filter_before_top3, int start, int end) {
  IntegerVector a_p = alpha_t.slot("p"), a_i = alpha_t.slot("i"), b_p = beta.slot("p"), b_i = beta.slot("i");
  NumericVector a_x = alpha_t.slot("x"), b_x = beta.slot("x");
  IntegerVector b_dim = beta.slot("Dim");
  const int n_alpha = a_p.length() - 1, n_beta = b_dim[0], n_wells = b_dim[1];
  check_chunk(start, end, n_alpha);
  SparseView alpha_view(a_p, a_i, a_x), beta_view(b_p, b_i, b_x);

  const int n_chunk = end - start;
  const int n_top = std::min(3, n_beta);
  IntegerVector out_i(n_chunk * n_top), out_j(n_chunk * n_top);
  NumericVector out_r(n_chunk * n_top);
  TshellWorker worker(alpha_view, beta_view, mean_a, norm_a, mean_b, norm_b, wa, wb, beta_order,
                      filter_before_top3, start, n_beta, n_wells, n_top, out_i, out_j, out_r);
  parallelFor(0, n_chunk, worker, 16);
  return List::create(_["i"] = out_i, _["j"] = out_j, _["r"] = out_r);
}
