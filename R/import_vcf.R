import_vcf <- function(vcf, MAF_thresh = 0.01){
    if(class(vcf) != "CollapsedVCF"){
        stop("Input file needs to be of class CollapsedVCF")
    }
    
    vcf_table <- info(vcf)    ### za frequency
    coords <- rowRanges(vcf)
    
    # obtain coordiantes of variants above a MAF_thresh
    common_vars <- row.names(vcf_table[vcf_table$LDAF > MAF_thresh, ])
    coords_common_vars <- coords[common_vars, ]
    variant_start <- start(ranges(coords_common_vars))
    variant_end <- end(ranges(coords_common_vars))
    
    # select only SNPs (same start and end position)
    snps <- variant_start == variant_end
    variant_start <- variant_start[snps]
    common_vars <- common_vars[snps]
    
    coords_common_vars <- coords_common_vars[common_vars, ]
    chromosomes <- as.character (seqnames (coords_common_vars))
    freq <- vcf_table[common_vars,"LDAF"]
    
    SnpData <- as.data.frame(cbind(chromosomes, variant_start, freq))
    colnames(SnpData) <- c("Chr_name", "Position_on_chr", "Minor_allele_freq")
    
    return(SnpData)
}