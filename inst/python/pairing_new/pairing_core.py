# Minimal GPU compute kernels for run_pairing().
#
# All data preparation, chunking, progress messages, and output assembly happens
# in R (see R/utils_run_pairing.R). This module does only the steps that run on a
# GPU (cupy for NVIDIA GPUs, mlx for Apple silicon), for one chunk of TCRalpha
# chains at a time:
#   - MAD-HYPE: count the wells shared by each alpha/beta pair and keep the pairs
#     that pass the MAD-HYPE look-up table
#   - T-SHELL: correlate the read fractions of each alpha/beta pair across wells
#     and keep the 3 best-correlated beta chains for each alpha chain
#
# The arithmetic (in float32) is the same as in inst/python/pairing/pairing_all_backends.py,
# so the results are identical.

import numpy as np


def _load_backend(backend):
    if backend == "cupy":
        import cupy as mx
    elif backend == "mlx":
        import mlx.core as mx
    else:
        raise ValueError("backend must be 'cupy' or 'mlx', not '" + str(backend) + "'")
    return mx


def check_backend(backend):
    """Raise an error if the GPU array library for 'backend' can't be loaded."""
    _load_backend(backend)
    return True


def _to_numpy(x):
    if hasattr(x, "get"):  # cupy
        return x.get()
    return np.array(x)  # mlx


class MadHype:
    """
    MAD-HYPE pairs for chunks of alpha chains.

    beta_bin: float32 array (n_beta, n_wells): 1 if the beta chain is in the well, otherwise 0
    b_total:  int32 array (n_beta,): number of wells with each beta chain
    mdh:      int32 array (n_wells + 1, n_wells + 1): MAD-HYPE look-up table. An alpha/beta pair
              with wij wells with both chains, wi wells with only alpha, and wj wells with only
              beta is called if wj < mdh[wij, wi].
    backend:  "cupy" or "mlx"
    """

    def __init__(self, beta_bin, b_total, mdh, backend):
        self.mx = _load_backend(backend)
        mx = self.mx
        self.beta_bin_t = mx.array(beta_bin).T  # (n_wells, n_beta)
        self.b_total = mx.array(b_total)[None, :]  # (1, n_beta)
        self.mdh = mx.array(mdh)

    def chunk(self, alpha_bin, a_total):
        """
        alpha_bin: float32 array (n_chunk, n_wells): 1 if the alpha chain is in the well, otherwise 0
        a_total:   int32 array (n_chunk,): number of wells with each alpha chain

        Returns a dict with one value per pair: "i" (row of alpha_bin), "j" (row of beta_bin),
        both 0-based, and "wij" (number of wells with both chains).
        """
        mx = self.mx
        a_total = mx.array(a_total)[:, None]  # (n_chunk, 1)
        overlaps = mx.matmul(mx.array(alpha_bin), self.beta_bin_t)  # wij for every pair
        wj = self.b_total - overlaps
        wi = a_total - overlaps
        mask = wj < self.mdh[overlaps.astype(mx.int32), wi.astype(mx.int32)]
        i, j = np.nonzero(_to_numpy(mask))
        wij = _to_numpy(overlaps[mx.array(i), mx.array(j)])
        return {"i": i.astype(np.int32), "j": j.astype(np.int32), "wij": wij}


class TShell:
    """
    T-SHELL pairs (3 best-correlated beta chains for each alpha chain) for chunks of alpha chains.

    beta:     float32 array (n_beta, n_wells): read fraction of each beta chain in each well
    b_total:  int32 array (n_beta,): number of wells with each beta chain
    filter_before_top3: if True, before choosing the top 3, set the correlation to -1 for pairs
              with <= 2 wells with both chains or with a loss fraction sum >= 0.5
    backend:  "cupy" or "mlx"
    """

    def __init__(self, beta, b_total, filter_before_top3, backend):
        self.mx = _load_backend(backend)
        mx = self.mx
        beta = mx.array(beta)
        beta_scaled = beta - mx.mean(beta, axis=1, keepdims=True)
        beta_scaled = beta_scaled / mx.linalg.norm(beta_scaled, ord=2, axis=1, keepdims=True)
        self.beta_scaled_t = beta_scaled.T  # (n_wells, n_beta)
        self.filter_before_top3 = filter_before_top3
        if filter_before_top3:
            self.beta_bin_t = (beta > 0).T.astype(mx.float32)
            self.b_total = mx.array(b_total)[None, :]

    def chunk(self, alpha, a_total):
        """
        alpha:   float32 array (n_chunk, n_wells): read fraction of each alpha chain in each well
        a_total: int32 array (n_chunk,): number of wells with each alpha chain

        Returns a dict with 3 values for each alpha chain (lowest correlation first): "i" (row of alpha),
        "j" (row of beta), both 0-based, and "r" (correlation).
        """
        mx = self.mx
        alpha = mx.array(alpha)
        alpha_scaled = alpha - mx.mean(alpha, axis=1, keepdims=True)
        alpha_scaled = alpha_scaled / mx.linalg.norm(alpha_scaled, ord=2, axis=1, keepdims=True)
        cors = mx.matmul(alpha_scaled, self.beta_scaled_t)
        if self.filter_before_top3:
            wij = mx.matmul((alpha > 0).astype(mx.float32), self.beta_bin_t)
            wa = mx.array(a_total)[:, None]
            wb = self.b_total
            loss_a_frac = (wb - wij) / (wij + (wb - wij) + (wa - wij))
            loss_b_frac = (wa - wij) / (wij + (wb - wij) + (wa - wij))
            mask = mx.logical_or(wij <= 2, (loss_a_frac + loss_b_frac) >= 0.5)
            cors = mx.where(mask, -1, cors)
        top3 = mx.argsort(cors, axis=1)[:, -3:]
        r = mx.take_along_axis(cors, top3, axis=1)
        n_chunk, n_top = top3.shape
        return {
            "i": np.repeat(np.arange(n_chunk, dtype=np.int32), n_top),
            "j": _to_numpy(top3).flatten().astype(np.int32),
            "r": _to_numpy(r).flatten(),
        }
