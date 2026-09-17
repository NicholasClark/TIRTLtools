// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// Minimal CPU compute kernel for TCRdist(..., backend = "cpp"). This mirrors
// tcrdist_chunk() in inst/python/TCRdist_new/TCRdist_core.py -- given a chunk
// of already-encoded TCRs, sum substitution-matrix penalties over all encoded
// features for every pair -- but runs in parallel on the CPU (RcppParallel)
// instead of calling out to python/cupy/mlx.
//
// tcr1_enc, tcr2_enc: n x k integer matrices of encoded TCR features (0-indexed
//   codes into submat, as produced by .encode_tcrs_new() in R/TCRdist.R)
// submat: square integer substitution/penalty matrix indexed by the encoded
//   feature values
//
// Returns a dense (n1 x n2) integer matrix of summed penalties.

namespace {
// mat1/mat2 here are K x n (features x TCRs) -- the transpose of the n x K
// tcr1_enc/tcr2_enc encoding -- so that, for a fixed TCR index, all K encoded
// features are contiguous in memory (R matrices are column-major, so walking
// k for a fixed row of an n x K matrix would stride by n on every access).
struct TCRdistChunkWorker : public Worker {
  const RMatrix<int> mat1_t;
  const RMatrix<int> mat2_t;
  const RMatrix<int> submat;
  RMatrix<int> out;
  const int K;
  const int n2;

  TCRdistChunkWorker(const IntegerMatrix& tcr1_enc_t, const IntegerMatrix& tcr2_enc_t,
                      const IntegerMatrix& submat, IntegerMatrix& out)
    : mat1_t(tcr1_enc_t), mat2_t(tcr2_enc_t), submat(submat), out(out),
      K(tcr1_enc_t.nrow()), n2(tcr2_enc_t.ncol()) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t idx = begin; idx < end; ++idx) {
      int i = idx / n2;
      int j = idx % n2;
      const int* col1 = mat1_t.column(i).begin();
      const int* col2 = mat2_t.column(j).begin();
      int s = 0;
      for (int k = 0; k < K; ++k) {
        s += submat(col1[k], col2[k]);
      }
      out(i, j) = s;
    }
  }
};
}

// [[Rcpp::export]]
IntegerMatrix tcrdist_chunk_cpp(const IntegerMatrix& tcr1_enc,
                                const IntegerMatrix& tcr2_enc,
                                const IntegerMatrix& submat) {
  int n1 = tcr1_enc.nrow();
  int n2 = tcr2_enc.nrow();
  int K = tcr1_enc.ncol();
  if (tcr2_enc.ncol() != K) {
    stop("tcr1_enc and tcr2_enc must have the same number of columns.");
  }

  // Transpose once per chunk (cheap: O(n*K)) so the hot (i,j) loop below can
  // read each TCR's encoded features contiguously (see TCRdistChunkWorker).
  IntegerMatrix tcr1_enc_t(K, n1);
  IntegerMatrix tcr2_enc_t(K, n2);
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < n1; ++i) tcr1_enc_t(k, i) = tcr1_enc(i, k);
    for (int j = 0; j < n2; ++j) tcr2_enc_t(k, j) = tcr2_enc(j, k);
  }

  IntegerMatrix out(n1, n2);
  TCRdistChunkWorker worker(tcr1_enc_t, tcr2_enc_t, submat, out);
  parallelFor(0, (std::size_t)n1 * (std::size_t)n2, worker);

  return out;
}
