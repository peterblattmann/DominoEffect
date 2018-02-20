identify_hotspots <- function(mutation_dataset, gene_data , snp_data,
                            min_n_muts = 5, MAF_thresh = 0.01,
                            flanking_region = c(200, 300),
                            poisson.thr = 0.01, percentage.thr = 0.15,
                            ratio.thr = 45, approach = "percentage"){
    
    ### Checking for argument requirements
    if(length(min_n_muts) != 1L){
        stop("min_n_muts should be an integer.")
    }
    if(!(length(MAF_thresh) == 1 & class(MAF_thresh) == "numeric" & 
         MAF_thresh < 1)){
        stop("MAF_thresh must representing a valid minor allele 
             threshold below 1.")
    }
    if(!(length(poisson.thr) == 1 & class(poisson.thr) == "numeric" & 
         poisson.thr < 1)){
        stop("poisson.thr must representing a valid probability 
             threshold below 1.")
    }
    if(!(length(percentage.thr) == 1 & class(percentage.thr) == "numeric" & 
         percentage.thr < 1)){
        stop("percentage.thr must representing a valid percentage 
             threshold below 1.")
    }
    if(!(length(ratio.thr) == 1 & class(ratio.thr) == "numeric" & 
         ratio.thr > 1)){
        stop("ratio.thr must representing a valid ratio 
             threshold above 1.")
    }
    if(!(approach %in% c("percentage", "ratio", "both"))){
        stop("Approach must be defined as either percentage, ratio, or both.")
    }
    
    
    MutationData.colnames <- c("Ensembl_gene_id","Protein_residue", "Original_aa",
                               "Mutated_aa","Patient_id","Genomic_coordinate",
                               "Original_base", "Mutated_base")
    gene_id <- MutationData.colnames[1]
    protein_res <- MutationData.colnames[2]
    original_aa <- MutationData.colnames[3]
    mutated_aa <- MutationData.colnames[4]
    patient_id <- MutationData.colnames[5]
    genomic_coord <- MutationData.colnames[6]
    
    
    if(sum(colnames(mutation_dataset) %in% MutationData.colnames) == 8){
        mutation_dataset <- unique(mutation_dataset[, MutationData.colnames])
    } else {
        stop("Mutation_dataset should contain the following columns: 
             Ensembl_gene_id, Protein_residue, Original_aa, Mutated_aa, 
             Patient_id, Genomic_coordinate, Original_base, Mutated_base.")
    }
    
    DominoData.colnames <- c("Ensembl_gene_id","Representative_tr", 
                             "cDNA_length", "Gene_name", "Uniprot_id")
    if(sum(colnames(gene_data) %in% DominoData.colnames) == 5){
        gene_data <- unique(gene_data[, DominoData.colnames])
    } else {
        stop("Gene_data should contain the following columns: Ensembl_gene_id, 
             Representative_tr, cDNA_length, Gene_name, Uniprot_id.")
    }

    if(!(length(flanking_region) == 2L | length(flanking_region) == 1L)){
        stop("Flanking region size needs to be one or two integers.")
    }

    # if gene_data is provided in TxDB format convert to data frame
    if(class(gene_data) == "TxDb"){
        gene_data <- import_txdb(gene_data)
    }
    
    # if SnpData is provided in VCF format convert to data frame
    if(class(snp_data) == "CollapsedVCF"){
        snp_data <- import_vcf(snp_data, MAF_thresh = MAF_thresh)
    }
    
    if(length(snp_data) > 0){
        SnpData.colnames <- c("Chr_name", "Position_on_chr", 
                              "Minor_allele_freq")
        
        if(sum(SnpData.colnames %in% colnames(snp_data)) != 3){
            stop("SNP table needs to contain the following columns: Chr_name, 
                 Position_on_chr, Minor_allele_freq.")
        }
        
        if(class(snp_data[,"Minor_allele_freq"]) != "numeric"){
            stop("Minor_allele_freq must be of type numeric.")
        }
        
        if(length(MAF_thresh) > 0){
            snp_data <- snp_data[snp_data[,"Minor_allele_freq"] > MAF_thresh, ]
        }
        
        snp_data <- snp_data[,SnpData.colnames]
        snp_pos <- paste(snp_data[,"Chr_name"], snp_data[,"Position_on_chr"], 
                         sep = ":")
        
        # filter out all snps
        muts_snps <- mutation_dataset[,genomic_coord] %in% snp_pos
        
        if(sum(muts_snps) > 0){
            mutation_dataset <- mutation_dataset[!muts_snps,]
            message("  ", sum(muts_snps), " mutations overlap with common 
                    population SNPs and were removed.")
        }
    } else {
        warning("No SNP-table provided.")
    }
    ###
    
    ns_mut_dataset <- mutation_dataset[
        as.character(mutation_dataset[,original_aa]) != 
            as.character(mutation_dataset[,mutated_aa]), ]
    ns_mut_dataset <- unique(ns_mut_dataset[,c(gene_id,protein_res,patient_id)])

    # count how many times mutation exists
    ns_mut_dataset_key <- paste(ns_mut_dataset[,gene_id],
                                     ns_mut_dataset[,protein_res], sep = "-")
    t <- table(ns_mut_dataset_key)
    muts_to_check <- names(t[t >= min_n_muts])

    message("  There are in total ", length(muts_to_check), 
            " protein residues with ",  min_n_muts, " or more mutation.")

    if(length(muts_to_check) == 0){
        stop("Exiting since there are no mutations to check.")
    }

    message("  ***Checking potential hotspots.***")
    results <- data.frame(stringsAsFactors = FALSE)
    field <- 1
    for (s in muts_to_check){
        g_poz <- unlist(strsplit(s, "-"))
        gene <- g_poz[1]
        mut_pos_numeric <- as.numeric(g_poz[2])
        tot_mut_hotsp <- as.numeric(t[s])
        
        # calculate amino acid length of protein
        length_cdna <- as.numeric(gene_data[gene_data$Ensembl_gene_id==gene, 
                                            "cDNA_length"])[1]
        if((length(length_cdna) == 0) || is.na(length_cdna)){
            warning("No cDNA length provided for ", gene, ".")
        } else {
            length_aa_in <- (length_cdna/3)-1
            length_aa <- round(length_aa_in)
            
            # retrieve gene information
            rep_tr <- as.character(gene_data[gene_data$Ensembl_gene_id==gene, 
                                             "Representative_tr"])
            
            gene_symbol <- unique(as.character(gene_data[gene_data$Ensembl_gene_id==gene, 
                                                          "Gene_name"]))
            assoc_unip <- unique(as.character(gene_data [gene_data$Ensembl_gene_id==gene, 
                                                         "Uniprot_id"]))
            
            if(length(gene_symbol) > 1){
                #warning("Several gene symbols for:, ", gene , paste(gene_symbol,collapse = ","))
                gene_symbol <- gene_symbol[1]
            }
            
            boundary <- calculate_boundary(mut_pos_numeric, length_aa, flanking_region)
            
            to_the_left1 <- boundary[["region_1"]][1]
            to_the_left2 <- boundary[["region_2"]][1]
            to_the_right1 <- boundary[["region_1"]][2]
            to_the_right2 <- boundary[["region_2"]][2]
            
            all_mut_1 = mutation_dataset[mutation_dataset[,gene_id] == gene & 
                                             mutation_dataset[,protein_res] >=  to_the_left1 & 
                                             mutation_dataset[,protein_res] <= to_the_right1, 
                                         c(gene_id,protein_res,patient_id)]
            all_mut_2 = mutation_dataset[mutation_dataset[,gene_id] == gene & 
                                             mutation_dataset[,protein_res] >= to_the_left2 & 
                                             mutation_dataset[,protein_res] <= to_the_right2, 
                                         c(gene_id,protein_res,patient_id)]
            
            tot_n_mut_1 = nrow(all_mut_1)
            tot_n_mut_2 = nrow(all_mut_2)
            
            ratio_region_1 <- tot_mut_hotsp/tot_n_mut_1
            ratio_region_2 <- tot_mut_hotsp/tot_n_mut_2
            
            exp_per_res_1 <- tot_n_mut_1/(to_the_right1 - to_the_left1 + 1)
            exp_per_res_2 <- tot_n_mut_2/(to_the_right2 - to_the_left2 + 1)
            
            reg_leng1 = to_the_right1 - to_the_left1 + 1
            p_excp_av1 = tot_n_mut_1/reg_leng1
            p_result1 = ppois(q = tot_mut_hotsp, lambda = p_excp_av1, 
                              lower.tail = FALSE)
            
            reg_leng2 = to_the_right2 - to_the_left2 + 1
            p_excp_av2 = tot_n_mut_2/reg_leng2
            p_result2 = ppois(q = tot_mut_hotsp, lambda = p_excp_av2, 
                              lower.tail = FALSE)
            
            # select only residues that meet poisson.thr, percentage.thr and ratio.thr
            if((p_result1 < poisson.thr) & (p_result2 < poisson.thr)){
                write.line <- FALSE
                
                if(approach == "percentage"){
                    if((ratio_region_1 > percentage.thr) & (ratio_region_2 > percentage.thr)){
                        write.line <- TRUE
                    }
                }
                else if(approach == "ratio"){
                    if(((tot_mut_hotsp/exp_per_res_1) > ratio.thr) & 
                       ((tot_mut_hotsp/exp_per_res_1) > ratio.thr)){
                        write.line <- TRUE
                    }
                }
                else if(approach == "both"){
                    if(((ratio_region_1 > percentage.thr) & 
                        (ratio_region_2 > percentage.thr)) &
                    (((tot_mut_hotsp/exp_per_res_1) > ratio.thr) & 
                     ((tot_mut_hotsp/exp_per_res_1) > ratio.thr))){
                        write.line <- TRUE
                    }
                }
                
                if(write.line){
                    region_1 <- paste(to_the_left1, to_the_right1, sep="-")
                    region_2 <- paste(to_the_left2, to_the_right2, sep="-")
                    
                    ratio_region_1 = sprintf("%1.2f%%", 100*ratio_region_1)
                    ratio_region_2 = sprintf("%1.2f%%", 100*ratio_region_2)
                    results[field, "Gene"] = gene_symbol
                    results[field, "Ensembl_id"] = gene
                    results[field, "Repres_tr"] = rep_tr
                    results[field, "Prot_position"] = mut_pos_numeric
                    results[field, "N_mut"] = tot_mut_hotsp
                    results[field, "Prot_length"] = length_aa
                    results[field, "Assoc_unip_ids"] = assoc_unip
                    results[field, "Perc_region_1"] = ratio_region_1
                    results[field, "Perc_region_2"] = ratio_region_2
                    results[field, "Poisson_p_value_1"] = p_result1
                    results[field, "Poisson_p_value_2"] = p_result2
                    field = field + 1
                }
            }
            
            if(nrow(results) > 0){
                results <-  results[order(results[,"Poisson_p_value_2"], decreasing = FALSE), ]
                results[,"Adj_p_value_region_1"] <- p.adjust(results[,"Poisson_p_value_2"], 
                                                             method = "BH")
                results[,"Adj_p_value_region_2"] <- p.adjust(results[,"Poisson_p_value_1"], 
                                                             method = "BH")
                
                results <- results[results$Adj_p_value_region_1 < poisson.thr & 
                                       results$Adj_p_value_region_2 < poisson.thr,]
            }
        }
    }
    #message ("  ***Identified hotspots are stored in the object: results.***")
    return(results)
}
