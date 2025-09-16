// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// all_mat1: n_tcrs1 x 60 (integer indices into rows of submat)
// all_mat2: n_tcrs2 x 60 (integer indices into cols of submat)
// submat  : R x C numeric substitution matrix
// tcr2_equals_tcr1: if true, fill diagonal+lower triangle with NA like in your R code
//
// Returns: n_tcrs1 x n_tcrs2 numeric matrix
//
// [[Rcpp::export]]
NumericMatrix tcrdist_cpp(const IntegerMatrix& all_mat1,
                          const IntegerMatrix& all_mat2,
                          const IntegerMatrix& submat,
                          const bool tcr2_equals_tcr1) {

  const int n1 = all_mat1.nrow();
  const int n2 = all_mat2.nrow();
  const int K  = all_mat1.ncol();        // expected 60
  if (all_mat2.ncol() != K) {
    stop("all_mat1 and all_mat2 must have the same number of columns (K).");
  }

  const int nrow_sub = submat.nrow();
  const int ncol_sub = submat.ncol();

  NumericMatrix out(n1, n2);

  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      if (tcr2_equals_tcr1 && i <= j) {
        out(i, j) = NA_REAL;
        continue;
      }
      double s = 0.0;
      // Sum over k: submat[ all_mat1[i,k], all_mat2[j,k] ]
      // Convert 1-based R indices to 0-based C++ indices.
      for (int k = 0; k < K; ++k) {
        int r = all_mat1(i, k) - 1;
        int c = all_mat2(j, k) - 1;
        s += submat(r, c);
      }
      out(i, j) = s;
    }
  }

  return out;
}





// [[Rcpp::export]]
NumericMatrix tcrdist_parallel(const IntegerMatrix& all_mat1,
                               const IntegerMatrix& all_mat2,
                               const NumericMatrix& submat,
                               const bool tcr2_equals_tcr1) {
  int n1 = all_mat1.nrow();
  int n2 = all_mat2.nrow();

  struct TcrWorker : public Worker {
    const RMatrix<int> mat1;
    const RMatrix<int> mat2;
    const RMatrix<double> submat;
    const bool tcr2_equals_tcr1;
    RMatrix<double> out;

    int K, nrow_sub;

    TcrWorker(const IntegerMatrix& all_mat1,
              const IntegerMatrix& all_mat2,
              const NumericMatrix& submat,
              bool tcr2_equals_tcr1,
              NumericMatrix& out)
      : mat1(all_mat1), mat2(all_mat2), submat(submat),
        tcr2_equals_tcr1(tcr2_equals_tcr1), out(out) {
      K = mat1.ncol();
      nrow_sub = submat.nrow();
    }

    void operator()(std::size_t begin, std::size_t end) {
      const int n2 = mat2.nrow();

      for (std::size_t idx = begin; idx < end; ++idx) {
        int i = idx / n2;   // row index
        int j = idx % n2;   // col index

        if (tcr2_equals_tcr1 && i <= j) {
          out(i, j) = NA_REAL;
          continue;
        }

        double s = 0.0;
        for (int k = 0; k < K; ++k) {
          int r = mat1(i, k) - 1;
          int c = mat2(j, k) - 1;
          s += submat(r, c);
        }
        out(i, j) = s;
      }
    }
  };

  NumericMatrix out(n1, n2);

  TcrWorker worker(all_mat1, all_mat2, submat, tcr2_equals_tcr1, out);

  parallelFor(0, n1 * n2, worker);

  return out;
}

