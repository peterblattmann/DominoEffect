import_txdb <- function(txdb_object){
    if(!is(txdb_object, "TxDb")){
        stop("Input file needs to be of class TxDb")
    }
    genes <- keys(txdb_object, keytype = "GENEID")
    EnsGenes <- select(txdb_object, keys = genes, columns = c("GENEID", "TXNAME", "TXTYPE", "CDSID", "CDSSTART", "CDSEND"), keytype = "GENEID")
    EnsGenes <- subset(EnsGenes, TXTYPE == "protein_coding")
    EnsGenes$CDS_length <- abs(EnsGenes$CDSEND - EnsGenes$CDSSTART + 1)
    
    EnsGenes_Transcripts <- aggregate(EnsGenes[,"CDS_length"], by=list(EnsGenes$GENEID, EnsGenes$TXNAME), function(x)sum(x, na.rm=TRUE))
    
    maxTrx <- aggregate(EnsGenes_Transcripts[,"x"], by=list(EnsGenes_Transcripts$Group.1), max)
    EnsGenes_longest_transcript <- merge(EnsGenes_Transcripts, maxTrx, by=c("Group.1", "x"))
    
    EnsGenes_longest_transcript <- EnsGenes_longest_transcript[,c(1,3,2)]
    
    colnames(EnsGenes_longest_transcript) <- c("Ensembl_gene_id", "Representative_tr", "cDNA_length")
    EnsGenes_longest_transcript$Gene_name <- NA
    EnsGenes_longest_transcript$Uniprot_id <- NA
    
    return(EnsGenes_longest_transcript)
}
