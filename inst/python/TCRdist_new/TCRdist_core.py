# Minimal GPU/CPU compute kernel for TCRdist_new().
#
# All TCR encoding, chunking, sparsification, and output assembly happens in R
# (see R/TCRdist_new.R). This module does only the one step that benefits from
# running on a GPU: summing substitution-matrix penalties over encoded TCR
# features for a chunk of TCRs.

import numpy as np


def tcrdist_chunk(tcr1_enc, tcr2_enc, submat, backend="numpy"):
    """
    Compute pairwise TCRdist (summed substitution penalties) between two chunks
    of already-encoded TCRs.

    tcr1_enc: integer array, shape (n1, n_features) -- encoded TCR features for chunk 1
    tcr2_enc: integer array, shape (n2, n_features) -- encoded TCR features for chunk 2
    submat:   square integer substitution/penalty matrix indexed by the encoded feature values
    backend:  "numpy", "cupy", or "mlx"

    Returns a dense (n1, n2) int32 numpy array of summed penalties.
    """
    if backend == "cupy":
        import cupy as mx
    elif backend == "mlx":
        import mlx.core as mx
    else:
        import numpy as mx

    a = mx.array(tcr1_enc)
    b = mx.array(tcr2_enc)
    s = mx.array(submat)

    result = mx.sum(s[a[:, None, :], b[None, :, :]], axis=2)

    if backend == "cupy":
        result = mx.asnumpy(result)
    else:
        result = np.array(result)

    return result.astype(np.int32)
