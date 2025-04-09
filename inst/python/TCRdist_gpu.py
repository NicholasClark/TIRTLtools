# load libraries/packages --------------

import utils
gpu = utils.check_gpu()
module = utils.check_cupy_or_mlx()

import numpy as np
import pandas as pd
import time
import scipy
import os

if gpu == "nvidia" and module == "cupy":
    print("Loading cupy to perform TCRdist")
    import cupy as mx #cuda python backend, connect to T4 runtime or other with GPU
elif gpu == "apple" and module == "mlx":
    print("Loading mlx to perform TCRdist")
    import mlx.core as mx #use this for apple silicon
else:
    print("Loading numpy to perform TCRdist")
    import numpy as mx #use this for CPU only

# helper functions --------------

# function to pad center with gaps until target length
def pad_center(seq, target_length):
   seq_length = len(seq)
   if seq_length >= target_length:
       return seq[:target_length]
   else:
       total_padding = target_length - seq_length
       first_half = seq[:seq_length // 2]
       second_half = seq[seq_length // 2:]
       return first_half + ['_'] * total_padding + second_half

def process_TCRs(tcr, params_vec, n_max=np.inf):
    n_max = min(n_max, tcr.shape[0])
    cdr3amat = np.array([pad_center(list(seq), 29) for seq in tcr['cdr3a'][:n_max] ])
    cdr3amatint = np.vectorize(params_vec.get)(cdr3amat)
    cdr3bmat = np.array([pad_center(list(seq), 29) for seq in tcr['cdr3b'][:n_max] ])
    cdr3bmatint = np.vectorize(params_vec.get)(cdr3bmat)
    
    cols_to_use = slice(3, -2) #truncate CDR3s
    
    encoded = np.column_stack([
        np.vectorize(params_vec.get)(tcr['va'][:n_max]),
        cdr3amatint[:,cols_to_use],
        np.vectorize(params_vec.get)(tcr['vb'][:n_max]),
        cdr3bmatint[:,cols_to_use]
    ])
    tcrs=mx.array(encoded).astype(mx.uint8)
    return(tcrs)

#### Naive no-loop TCRdist function for GPU with sparse output -- either sparse csr_matrix or pandas dataframe
#### Returns only edges with TCRdist less than or equal to cutoff (default = 90)
#### Returns dataframe with 3 columns: 'row' (row_index), 'col' (column index), and 'TCRdist' (TCRdist value)
#### Can also return a sparse matrix of type scipy.sparse.csr_matrix 
def TCRdist_inner(tcr1, tcr2, submat, tcrdist_cutoff=90, 
                  ch1=0, ch2=0, output = "edge_list",
                  only_lower_tri = True,
                  compare_to_self = False):
    result = mx.sum(submat[tcr1[:, None, :], tcr2[ None,:, :]],axis=2)
    ### set values with TCRdist == 0 to negative 1 so that they are not lost when converted to sparse matrix
    mask = result == 0
    mask = mask*(-1)
    result = result+mask
    ## set values greater than the cutoff to zero
    less_or_equal = mx.less_equal(result, tcrdist_cutoff)
    result = result*less_or_equal
    if mx.__name__ == "cupy":
        result = mx.asnumpy(result)
    if output in ["sparse", "both", "edge_list"]:
        #score_dtype = np.int16
        score_dtype = np.int32
        #if tcrdist_cutoff <= 255:
        #    score_dtype = np.uint16
        ### convert matrix to sparse (gets rid of all zero values)
        result_sparse = scipy.sparse.csr_matrix(result, dtype = score_dtype)
    if output in ["edge_list", "both"]:
        ### convert matrix to dataframe with indices and TCRdist values
        coo_mat = result_sparse.tocoo()
        df = pd.DataFrame({'edge1_0index': coo_mat.row+ch1, 'edge2_0index': coo_mat.col+ch2, 'TCRdist': coo_mat.data})
        if compare_to_self:
            if only_lower_tri:
                df = df[df['edge1_0index'] > df['edge2_0index']]
            else:
                df = df[df['edge1_0index'] != df['edge2_0index']]
        ### Change values of -1 back to zero (see comment above)
        df['TCRdist'] = df['TCRdist'].replace(-1, 0)
    if output == "both":
        return(result_sparse, df)
    elif output == "edge_list":
        return(df)
    elif output == "sparse":
        return(result_sparse)
    

### TCRdist function for GPU with batching for tcr1 and tcr2 lists and sparse output
### Returns a pandas data frame of all edges with TCRdist less than cutoff (default = 90)
### Returns dataframe with 3 columns: 'row' (row_index), 'col' (column index), and 'TCRdist' (TCRdist value)
def TCRdist_batch(tcr1, submat, params_df, tcr2=None, tcrdist_cutoff=90, chunk_size=1000, chunk_size_col = None, print_chunk_size=1000, print_res = True, only_lower_tri = True):
    #chunk_size = np.int64(chunk_size)
    #print_chunk_size = np.int64(print_chunk_size)
    compare_to_self = False
    submat = mx.array(submat, dtype = mx.uint8)
    params_vec = dict(zip(params_df["feature"], params_df["value"]))
    tcr1_mx = process_TCRs(tcr1, params_vec=params_vec)
    tcr1 = tcr1.copy()
    tcr1['tcr_index'] = range(len(tcr1))
    if tcr2 is None:
        compare_to_self = True
        tcr2_mx = tcr1_mx
    else:
        tcr2 = tcr2.copy()
        tcr2_mx = process_TCRs(tcr2, params_vec=params_vec)
        tcr2['tcr_index'] = range(len(tcr2))
    n1 = tcr1_mx.shape[0]
    n2 = tcr2_mx.shape[0]
    num_chunks1 = np.int64(np.ceil(n1//chunk_size))
    if chunk_size_col is None:
        chunk_size_col = chunk_size
    else:
        chunk_size_col = min(chunk_size_col, n2)
    num_chunks2 = np.int64(np.ceil(n2//chunk_size_col))
    if print_res:
        print('total number of chunks (rows):', num_chunks1)
        print('total number of chunks (cols):', num_chunks2)
    start_time = time.time()
    res_list = []
    for ch in range(0, n1, chunk_size):
        if print_res:
            if ch % print_chunk_size == 0:
                print('Processing chunk (rows)', ch)
        chunk_end = min(ch + chunk_size, n1)
        row_range1 = slice(ch, chunk_end)
        tcr1_tmp = tcr1_mx[row_range1,:]
        for ch2 in range(0, n2, chunk_size_col):
            if compare_to_self and ch < ch2 and only_lower_tri:
                continue
            chunk_end2 = min(ch2 + chunk_size_col, n2)
            row_range2 = slice(ch2, chunk_end2)
            tcr2_tmp = tcr2_mx[row_range2,:]
            edges_tmp = TCRdist_inner(tcr1=tcr1_tmp, tcr2=tcr2_tmp, submat=submat,
                                       tcrdist_cutoff=tcrdist_cutoff,
                                       ch1=ch, ch2=ch2, output="edge_list", only_lower_tri = only_lower_tri,
                                       compare_to_self = compare_to_self)
            res_list.append(edges_tmp)
    res = pd.concat(res_list)
    res.reset_index(inplace=True)
    res = res.drop('index', axis=1)
    end_time = time.time()
    if print_res:
        res
        print(f"Time taken: {end_time - start_time:.6f} seconds")
    if compare_to_self:
        res_dict = {
            'TCRdist_df': res,
            'tcr1': tcr1
        }
    else:
        res_dict = {
            'TCRdist_df': res,
            'tcr1': tcr1,
            'tcr2': tcr2
        }
    return(res_dict)
