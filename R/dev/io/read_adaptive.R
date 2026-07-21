
path = "~/git/newell_benchmarking/data/BM03_processed-data/adaptive/A-F5_TRB.tsv"
read_adaptive = function(path) {
  df = fread(path) %>% as_tibble() %>%
    select(templates, frequency, amino_acid, cdr3_amino_acid, ,frame_type, productive_frequency, cdr3_length, v_gene, j_gene, rearrangement, cdr3_rearrangement, everything())

}
