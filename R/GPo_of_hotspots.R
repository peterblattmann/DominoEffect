GPo_of_hotspots <- function(hotspot_mutations){
    hotspot_mutations$AA_position <- paste(hotspot_mutations$Ensembl_id, 
                                           "-", hotspot_mutations$Prot_position, sep = "")
    TestData$AA_position <- paste(TestData$Ensembl_gene_id, "-", 
                                  TestData$Protein_residue, sep = "")
    cols_of_interest <- c("Genomic_coordinate", "Original_base", "Mutated_base", "AA_position")
    hotspots_genomic_info <- unique(TestData[TestData$AA_position %in% 
                                                 hotspot_mutations$AA_position, cols_of_interest])
    coordinates_sep <- data.frame(do.call('rbind', 
                                          strsplit(as.character(hotspots_genomic_info$Genomic_coordinate),
                                                   ':',fixed=TRUE)))
    prot_residue_info <- data.frame(do.call('rbind', 
                                            strsplit(as.character(hotspots_genomic_info$AA_position),
                                                     '-',fixed=TRUE)))
    chr_info = paste("chr", coordinates_sep$X1)
    position_on_chr = coordinates_sep$X2
    hotspots_GPo <- GPos(seqnames=chr_info, pos= position_on_chr, stitch=FALSE)
    
    GenomicRanges::mcols(hotspots_GPo)$REF_NT <- as.character(hotspots_genomic_info$Original_base)
    GenomicRanges::mcols(hotspots_GPo)$MUT_NT <- as.character(hotspots_genomic_info$Mutated_base)
    GenomicRanges::mcols(hotspots_GPo)$Gene_ID <- as.character(prot_residue_info$X1)
    GenomicRanges::mcols(hotspots_GPo)$Protein_pos <- as.character(prot_residue_info$X2)
    
    return(hotspots_GPo)
}
    

