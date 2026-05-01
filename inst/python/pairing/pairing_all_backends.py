import utils
import numpy as np
import pandas as pd
from datetime import datetime
import sys
import os
## Note: either cupy, mlx.core, or numpy is loaded as "mx"

### if backend is specified, load it
# if backend == "cupy":
#   import cupy as mx
# elif backend == "mlx":
#   import mlx.core as mx
# elif backend == "numpy":
#   import numpy as mx
# else: ### if backend is not specified, check and load correct one


gpu = utils.check_gpu()
module = utils.check_cupy_or_mlx()

if gpu == "nvidia" and module == "cupy":
    print("Loading cupy")
    import cupy as mx #cuda python backend, connect to T4 runtime or other with GPU
elif gpu == "apple" and module == "mlx":
    print("Loading mlx")
    import mlx.core as mx #use this for apple silicon
else:
    print("Loading numpy")
    import numpy as mx #use this for CPU only

def pairing(prefix, folder_out, bigmas, bigmbs, mdh, backend="auto", filter_before_top3 = False, read_files = False, write_files = True, chunk_size = 500) :
  if read_files:
    print("start load:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    bigmas = mx.array(np.loadtxt(os.path.join(folder_out, prefix+'_bigmas.tsv'), delimiter='\t', dtype=np.float32))
    bigmbs = mx.array(np.loadtxt(os.path.join(folder_out, prefix+'_bigmbs.tsv'), delimiter='\t', dtype=np.float32))
    mdh = mx.array(np.loadtxt(os.path.join(folder_out, prefix+'_mdh.tsv'), delimiter='\t', dtype=np.int32))
  else:
    bigmas = mx.array(bigmas)
    bigmbs = mx.array(bigmbs)
    mdh = mx.array(mdh)
  mdh = madhyper_process(prefix = prefix, folder_out = folder_out, bigmas = bigmas, bigmbs = bigmbs, mdh = mdh, write_files = write_files, chunk_size = chunk_size)
  corr = correlation_process(prefix = prefix, folder_out = folder_out, bigmas = bigmas, bigmbs = bigmbs, filter_before_top3 = filter_before_top3, write_files = write_files, chunk_size = chunk_size)
  res = {"mdh": mdh, "corr": corr}
  return(res)


def top3_indices(arr):
    # Get indices of top 3 values per row using argsort (descending)
    top3 = mx.argsort(arr, axis=1)[:, -3:]  # shape: (500, 3)
    
    # Build row indices: [0,0,0, 1,1,1, 2,2,2, ...]
    row_indices = mx.repeat(mx.arange(arr.shape[0]), 3)  # shape: (1500,)
    
    # Flatten col indices
    col_indices = top3.flatten()  # shape: (1500,)
    
    # Stack into (1500, 2) pairs
    result = mx.stack([row_indices, col_indices], axis=1)
    
    return result
  
def to_numpy(x):
  if hasattr(x, 'get'):  # CuPy
      return(x.get())
  return(np.array(x))  # NumPy or MLX

def madhyper_process(prefix, folder_out, bigmas, bigmbs, mdh, write_files = False, chunk_size = 500):
    #print("start load:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    #bigmas = mx.array(np.loadtxt(os.path.join(folder_out, prefix+'_bigmas.tsv'), delimiter='\t', dtype=np.float32))
    #bigmbs = mx.array(np.loadtxt(os.path.join(folder_out, prefix+'_bigmbs.tsv'), delimiter='\t', dtype=np.float32))
    #mdh = mx.array(np.loadtxt(os.path.join(folder_out, prefix+'_mdh.tsv'), delimiter='\t', dtype=np.int32))
    rowinds_bigmas=mx.arange(bigmas.shape[0])
    rowinds_bigmbs=mx.arange(bigmbs.shape[0])
    results = []
    #chunk_size =500  # Define your chunk size
    chunk_size = int(chunk_size)
    n_wells = bigmas.shape[1]  # Assuming n_wells is the number of columns in bigmas
    print('total number of chunks', bigmas.shape[0]//chunk_size)
    total_chunks = (bigmas.shape[0] // chunk_size)+1
    b_total = mx.sum(bigmbs > 0, axis=1,keepdims=True)
    bigmbs=(bigmbs > 0).T.astype(mx.float32)
    print("start time for MAD-HYPE:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for ch in range(0, bigmas.shape[0], chunk_size):
        percent_complete = int((ch // chunk_size + 1) / total_chunks * 100)
        # Print progress only on 5%, 10%, etc.
        if percent_complete % 10 == 0 and percent_complete > 0:
            print(f'Progress: {ch} ({percent_complete}%)')
        chunk_end = min(ch + chunk_size, bigmas.shape[0])
        row_range = slice(ch, chunk_end)
        a_total = mx.sum(bigmas[row_range] > 0, axis=1,keepdims=True)
        overlaps = mx.matmul((bigmas[row_range] > 0).astype(mx.float32), bigmbs) #optimized
        mask_condition=-(overlaps.T - b_total).T < mdh[(overlaps).astype(mx.int16), -(overlaps - a_total).astype(mx.int16)]# mad hype only
        pairs = mx.array(np.argwhere(to_numpy(mask_condition)))
        result = {
              'alpha_nuc': 1+to_numpy(rowinds_bigmas[row_range][pairs[:, 0]]),
              'beta_nuc': 1+to_numpy(rowinds_bigmbs[pairs[:, 1]]),
              'wij': to_numpy(overlaps[pairs[:, 0], pairs[:, 1]]),
              'wa': to_numpy(a_total[:,0][pairs[:, 0]]),
              'wb': to_numpy(b_total[:,0][pairs[:, 1]])
              }
        results.append(result)
#result is a list of dictionaries, each dictionary contains the results for a chunk of rows. You can convert it to a pandas DataFrame like this: 
    print("end time for MAD-HYPE:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

#make pandas dataframe from each element of results and concatenate them
    results_df = pd.concat([pd.DataFrame(result) for result in results])
    results_df = results_df.reset_index(drop=True)
    if write_files:
      results_df.to_csv(os.path.join(folder_out, prefix+'_madhyperesults.csv'), index=False)
    print(f"Number of pairs: {results_df.shape[0]}")
    return(results_df)

def correlation_process(prefix, folder_out, bigmas, bigmbs, min_wells=2, filter_before_top3=False, write_files = False, chunk_size = 500):
    #print("start load:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    #bigmas = mx.array(np.loadtxt(os.path.join(folder_out,prefix+'_bigmas.tsv'), delimiter='\t', dtype=np.float32))
    #bigmbs = mx.array(np.loadtxt(os.path.join(folder_out, prefix+'_bigmbs.tsv'), delimiter='\t', dtype=np.float32))
    #mdh = mx.array(np.loadtxt(prefix+'_mdh.tsv', delimiter='\t', dtype=np.int32))
    rowinds_bigmas=mx.arange(bigmas.shape[0])
    rowinds_bigmbs=mx.arange(bigmbs.shape[0])
    
    #now we need to downsize to min_wells
    non_zero_counts_bigmas = mx.sum(bigmas > 0, axis=1)
    non_zero_counts_bigmbs = mx.sum(bigmbs > 0, axis=1)
    # Find the rows that have more than min_wells non-zero elements
    # valid_rows_bigmas = mx.nonzero(non_zero_counts_bigmas > min_wells)[0]
    # valid_rows_bigmbs = mx.nonzero(non_zero_counts_bigmbs > min_wells)[0]
    valid_rows_bigmas = [i for i, m in enumerate((non_zero_counts_bigmas > min_wells).tolist()) if m]
    valid_rows_bigmbs = [i for i, m in enumerate((non_zero_counts_bigmbs > min_wells).tolist()) if m]
    # # Filter the bigmas and bigmbs arrays
    bigmas = bigmas[valid_rows_bigmas]
    bigmbs = bigmbs[valid_rows_bigmbs]
    # # Also retain the corresponding indices
    rowinds_bigmas = rowinds_bigmas[valid_rows_bigmas]
    rowinds_bigmbs = rowinds_bigmbs[valid_rows_bigmbs]
    
    results = []
    chunk_size = int(chunk_size)
    n_wells = bigmas.shape[1]  # Assuming n_wells is the number of columns in bigmas
    print('total number of chunks', bigmas.shape[0]//chunk_size)
    total_chunks = (bigmas.shape[0] // chunk_size)+1
    bigmb_w1_scaled = bigmbs - mx.mean(bigmbs, axis=1, keepdims=True)
    bigmb_w1_scaled = (bigmb_w1_scaled / mx.linalg.norm(bigmb_w1_scaled,ord=2,axis=1, keepdims=True)).T
    b_total = mx.sum(bigmbs > 0, axis=1,keepdims=True)
    bigmbs=(bigmbs > 0).T.astype(mx.float32)
    print("start processing time for T-Shell:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for ch in range(0, bigmas.shape[0], chunk_size):
        #print('Processing chunk', ch)
        percent_complete = int((ch // chunk_size + 1) / total_chunks * 100)
        # Print progress only on 5%, 10%, etc.
        if percent_complete % 10 == 0 and percent_complete > 0:
            print(f'Progress: {ch} ({percent_complete}%)')
        chunk_end = min(ch + chunk_size, bigmas.shape[0])
        row_range = slice(ch, chunk_end)
        bigma_w1_scaled = bigmas[row_range] - mx.mean(bigmas[row_range], axis=1, keepdims=True) # mask for madhype
        bigma_w1_scaled = bigma_w1_scaled / mx.linalg.norm(bigma_w1_scaled,ord=2,axis=1, keepdims=True) # mask for madhype
        a_total = mx.sum(bigmas[row_range] > 0, axis=1,keepdims=True)
        pairwise_cors_method2 = mx.matmul(bigma_w1_scaled, bigmb_w1_scaled) #mask for madhype
        overlaps = mx.matmul((bigmas[row_range] > 0).astype(mx.float32), bigmbs) #optimized
        
        if filter_before_top3:
            ### remove correlations for pairs with low overlap or high loss fraction
            overlap_mask = overlaps <= 2 # require overlap of >=3 wells
            ## calculate loss fraction
            wij = overlaps 
            wa = a_total
            wb = b_total.T
            loss_a_frac =(wb-wij)/(wij+(wb-wij)+(wa-wij))
            loss_b_frac =(wa-wij)/(wij+(wb-wij)+(wa-wij))
            loss_frac_sum = loss_a_frac+loss_b_frac
            loss_frac_mask = loss_frac_sum >= 0.5 # require (loss_a_frac+loss_b_frac)<0.5
            combined_mask = mx.logical_or(overlap_mask, loss_frac_mask)
            pairwise_cors_method2 = mx.where(combined_mask, -1, pairwise_cors_method2) # set correlation to -1 if not enough overlap or loss fraction is too high
        pairs = top3_indices(pairwise_cors_method2)
        result = {
            'alpha_nuc': 1 + to_numpy(rowinds_bigmas[row_range][pairs[:, 0]]),
            'beta_nuc': 1 + to_numpy(rowinds_bigmbs[pairs[:, 1]]),
            'r': to_numpy(pairwise_cors_method2[pairs[:, 0], pairs[:, 1]]),
            'wij': to_numpy(overlaps[pairs[:, 0], pairs[:, 1]]),
            'wa': to_numpy(a_total[:,0][pairs[:, 0]]),
            'wb': to_numpy(b_total[:,0][pairs[:, 1]])
        }
        results.append(result)
        # result_df = pd.DataFrame(result)
        # print('number of rows in this chunk: ' + str(result_df.shape[0]))
        # tmp_df = pd.concat([pd.DataFrame(result) for result in results])
        # print('number of rows total: ' + str(tmp.shape[0]))
        
#result is a list of dictionaries, each dictionary contains the results for a chunk of rows. You can convert it to a pandas DataFrame like this: 
    print("end time for T-Shell:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

#make pandas dataframe from each element of results and concatenate them
    results_df = pd.concat([pd.DataFrame(result) for result in results])
    results_df = results_df.reset_index(drop=True) 
    if write_files:
      results_df.to_csv(os.path.join(folder_out, prefix+'_corresults.csv'), index=False)
    return(results_df)

