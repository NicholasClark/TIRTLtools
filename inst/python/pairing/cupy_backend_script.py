import cupy as mx
import numpy as np
import pandas as pd
from datetime import datetime
import sys

def madhyper_process(prefix):
    print("start load:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    bigmas = mx.array(np.loadtxt(prefix+'_bigmas.tsv', delimiter='\t', dtype=np.float32))
    bigmbs = mx.array(np.loadtxt(prefix+'_bigmbs.tsv', delimiter='\t', dtype=np.float32))
    mdh = mx.array(np.loadtxt(prefix+'_mdh.tsv', delimiter='\t', dtype=np.int32))
    rowinds_bigmas=mx.arange(bigmas.shape[0])
    rowinds_bigmbs=mx.arange(bigmbs.shape[0])
    results = []
    chunk_size =500  # Define your chunk size
    n_wells = bigmas.shape[1]  # Assuming n_wells is the number of columns in bigmas
    print('total number of chunks', bigmas.shape[0]//chunk_size)
    total_chunks = bigmas.shape[0] // chunk_size
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
        pairs = mx.argwhere(mask_condition)
        #pairs = mx.array(np.argwhere(np.array(mask_condition, copy=False)))
    
        result = { #this works in cupy
              'alpha_nuc': 1+(rowinds_bigmas[row_range][pairs[:, 0]]).get(),
              'beta_nuc': 1+(rowinds_bigmbs[pairs[:, 1]]).get(),
              'wij': (overlaps[pairs[:, 0], pairs[:, 1]]).get(),
              'wa': (a_total[:,0][pairs[:, 0]]).get(),
              'wb': (b_total[:,0][pairs[:, 1]]).get()
        }

        #result = {
         #   'alpha_nuc': 1+np.array(rowinds_bigmas[row_range][pairs[:, 0]]),
          #  'beta_nuc': 1+np.array(rowinds_bigmbs[pairs[:, 1]]),
    #       'r': np.array(pairwise_cors_method2[pairs[:, 0], pairs[:, 1]]),
    #     'pval': pairwise_cors_method2_ps[pairs[:, 0], pairs[:, 1]],
    #     'pval3': pairwise_corsn[pairs[:, 0], pairs[:, 1]],
         #   'wij': np.array(overlaps[pairs[:, 0], pairs[:, 1]]),
        #    'wa': np.array(a_total[:,0][pairs[:, 0]]),
        #    'wb': np.array(b_total[:,0][pairs[:, 1]])
        #}

        results.append(result)
#result is a list of dictionaries, each dictionary contains the results for a chunk of rows. You can convert it to a pandas DataFrame like this: 
    print("end time for MAD-HYPE:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

#make pandas dataframe from each element of results and concatenate them
    results_df = pd.concat([pd.DataFrame(result) for result in results])
    results_df.to_csv(prefix+'_madhyperesults.csv', index=False)
    print(f"Number of pairs from MAD-HYPE: {results_df.shape[0]}")

def correlation_process(prefix,min_wells=2,filter_before_top3=False):
    print("start load for T-Shell:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    bigmas = mx.array(np.loadtxt(prefix+'_bigmas.tsv', delimiter='\t', dtype=np.float32))
    bigmbs = mx.array(np.loadtxt(prefix+'_bigmbs.tsv', delimiter='\t', dtype=np.float32))
    #mdh = mx.array(np.loadtxt(prefix+'_mdh.tsv', delimiter='\t', dtype=np.int32))
    rowinds_bigmas=mx.arange(bigmas.shape[0])
    rowinds_bigmbs=mx.arange(bigmbs.shape[0])
    print('file read done')
     #now we need to downsize to min_wells. Following filter does not work in mlx. add numpy workaround?
    non_zero_counts_bigmas = mx.sum(bigmas > 0, axis=1)
    non_zero_counts_bigmbs = mx.sum(bigmbs > 0, axis=1)
     # Find the rows that have more than min_wells non-zero elements
    valid_rows_bigmas = mx.nonzero(non_zero_counts_bigmas > min_wells)[0]
    valid_rows_bigmbs = mx.nonzero(non_zero_counts_bigmbs > min_wells)[0]
    # # Filter the bigmas and bigmbs arrays
    bigmas = bigmas[valid_rows_bigmas]
    bigmbs = bigmbs[valid_rows_bigmbs]
    # # Also retain the corresponding indices
    rowinds_bigmas = rowinds_bigmas[valid_rows_bigmas]
    rowinds_bigmbs = rowinds_bigmbs[valid_rows_bigmbs]
    
    results = []
    chunk_size =500  # Define your chunk size
    n_wells = bigmas.shape[1]  # Assuming n_wells is the number of columns in bigmas
    print('total number of chunks', bigmas.shape[0]//chunk_size)
    total_chunks = bigmas.shape[0] // chunk_size
    bigmb_w1_scaled = bigmbs - mx.mean(bigmbs, axis=1, keepdims=True) #uncomment for cor
    bigmb_w1_scaled = (bigmb_w1_scaled / mx.linalg.norm(bigmb_w1_scaled,ord=2,axis=1, keepdims=True)).T #uncomment for cor
    b_total = mx.sum(bigmbs > 0, axis=1,keepdims=True)
    bigmbs=(bigmbs > 0).T.astype(mx.float32)
    print("start processing time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
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

        mask_condition=(-pairwise_cors_method2<=mx.partition(pairwise_cors_method2*(-1), 3,axis=1)[:,2:3])
 
        pairs = mx.argwhere(mask_condition)# this is for cupy!
        #pairs = mx.array(np.argwhere(np.array(mask_condition, copy=False))) #this should work with mlx

        result = { #this works in cupy
              'alpha_nuc': 1+(rowinds_bigmas[row_range][pairs[:, 0]]).get(),
              'beta_nuc': 1+(rowinds_bigmbs[pairs[:, 1]]).get(),
              'r': (pairwise_cors_method2[pairs[:, 0], pairs[:, 1]]).get(),
              'wij': (overlaps[pairs[:, 0], pairs[:, 1]]).get(),
              'wa': (a_total[:,0][pairs[:, 0]]).get(),
              'wb': (b_total[:,0][pairs[:, 1]]).get()
        }

        #result = { #this works in mlx
        #    'alpha_nuc': 1+np.array(rowinds_bigmas[row_range][pairs[:, 0]]),
        #    'beta_nuc': 1+np.array(rowinds_bigmbs[pairs[:, 1]]),
        #    'r': np.array(pairwise_cors_method2[pairs[:, 0], pairs[:, 1]]),
    #     'pval': pairwise_cors_method2_ps[pairs[:, 0], pairs[:, 1]],
    #     'pval3': pairwise_corsn[pairs[:, 0], pairs[:, 1]],
        #    'wij': np.array(overlaps[pairs[:, 0], pairs[:, 1]]),
        #    'wa': np.array(a_total[:,0][pairs[:, 0]]),
        #    'wb': np.array(b_total[:,0][pairs[:, 1]])
        #}
        results.append(result)
#result is a list of dictionaries, each dictionary contains the results for a chunk of rows. You can convert it to a pandas DataFrame like this: 
    print("end time for T-Shell:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

#make pandas dataframe from each element of results and concatenate them
    results_df = pd.concat([pd.DataFrame(result) for result in results])
    results_df.to_csv(prefix+'_corresults.csv', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <prefix>")
        sys.exit(1)
        
    prefix = sys.argv[1]    
    madhyper_process(prefix)
    correlation_process(prefix)
