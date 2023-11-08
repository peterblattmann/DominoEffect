DominoEffect <- function(mutation_dataset, gene_data, snp_data,
                         min_n_muts = 5, MAF_thresh = 0.01,
                         flanking_region = c(200, 300),
                         poisson.thr = 0.01, percentage.thr = 0.15,
                         ratio.thr = 45, approach = "percentage", write_to_file = "NO",
                         ens_release = "https://feb2023.archive.ensembl.org"){
  
  results <- identify_hotspots(mutation_dataset, gene_data, snp_data,
                               min_n_muts, MAF_thresh, flanking_region,
                               poisson.thr, percentage.thr, ratio.thr, 
                               approach)
  
  hotspot_mutations <- map_to_func_elem(results, write_to_file, 
                                        ens_release = ens_release)
  message ("   ***FINISHED***")
  return(hotspot_mutations)
}
