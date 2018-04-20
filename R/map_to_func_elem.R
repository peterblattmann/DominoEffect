map_to_func_elem <- function(hotspot_results, write_to_file = "NO", ens_release = "73"){
    ## Checking for argument requirements
    if(!(write_to_file %in% c("YES", "NO"))){
        stop("write_to_file should either be YES or NO")
    }
    ###
    
    message ("  Obtaining sequences of proteins encoded by the representative Ensembl transcripts.")
    
    if (ens_release == 73) {
        #ens_database <- "sep2013.archive.ensembl.org"
        ens_database <- "dec2013.archive.ensembl.org"
    } else {
        ens_database <- ens_release
    }
    
    #####
    # to fail gracefully if Ensembl is down
    # code from http://bioconductor.org/developers/how-to/web-query/
    N.TRIES <- as.integer(5)
    stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))
    
    while (N.TRIES > 0L) {
        ensembl <- tryCatch(useMart(biomart="ENSEMBL_MART_ENSEMBL", host = ens_database, dataset = "hsapiens_gene_ensembl"), error=identity)
        if (!inherits(ensembl, "error"))
            break
        N.TRIES <- N.TRIES - 1L
    }
    
    if (N.TRIES == 0L) {
        stop("'useMart()' failed:",
             "\n  error: ", conditionMessage(ensembl))
    }
    #####
    # ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl") # --- for current version
    
    transc_ids = hotspot_results$Repres_tr
    prot_seq = getSequence (id = transc_ids, type= "ensembl_transcript_id", seqType = "peptide", mart = ensembl)
    
    gff.file <- list()
    uni.sequences <- c()
    
    message ("  Obtaining UniProtKB sequences and coordinates of protein functional elements.")
    genes_w_hotspots <- unique(hotspot_results$Ensembl_id)
    for(gene_ens in genes_w_hotspots){
        gene_unip <- hotspot_results [hotspot_results$Ensembl_id==gene_ens, "Assoc_unip_ids"]
        assoc_unip_ids <- as.character(gene_unip)

        if(length(grep(",", assoc_unip_ids))){
            assoc_unip_ids <- trimws(unlist(strsplit(assoc_unip_ids, ",")))
        }
            
        for (unip_acc in assoc_unip_ids){
            if (!(is.na(unip_acc)))
            {
            if(nchar(unip_acc) == 6){
                retrieve_unip <- paste("http://www.uniprot.org/uniprot/", unip_acc, ".gff", sep = "")
                gff.file[[unip_acc]] <- read.delim(retrieve_unip, skip = 2, header = FALSE)
                    
                retrieve_unip <- paste ("http://www.uniprot.org/uniprot/", unip_acc, ".fasta", sep = "")
                uni.sequences <- c(uni.sequences, readAAStringSet(retrieve_unip))
                    
            }else {
                warning("IDs for following uniprot accessions could not be retrieved:", unip_acc, " that corresponds to:", gene_ens, ".")
            }
            }
        }
    }
    
    gff.file <- lapply(gff.file, function(x)x[,c(1,2,3,4,5,9)])
    gff.file <- lapply(gff.file, function(x){
        colnames(x) <- c("ID", "Source", "Type", "Start", "End", "Description"); x
    })
    
    uni.sequences <- do.call(c, uni.sequences)
    
    message ("  Assessing sequence agreement and mapping hotspots to functional elements.")
    pb <- txtProgressBar(min = 0, max = nrow(hotspot_results))

    for(i in seq_len(nrow(hotspot_results))){
        setTxtProgressBar(pb, i)
        #print(i)
        uni_genes <- hotspot_results[i, "Assoc_unip_ids"]
        ensembl_gene <- hotspot_results[i, "Ensembl_id"]
        ensembl_transcript <- hotspot_results[i, "Repres_tr"]
        ensembl_mut_position <- hotspot_results[i, "Prot_position"]
        ens.seq <- prot_seq[prot_seq$ensembl_transcript_id==ensembl_transcript, "peptide"]
        ids.uniprot <- c()
        all_unip <- unlist(strsplit(uni_genes, ","))
        for (uni_gene in all_unip)
        {
            id.unip <- grep(uni_gene, names(uni.sequences))
            ids.uniprot <- c(ids.uniprot, id.unip)
        }
        #### nema if, samp for loop za sve uniprot ids i posli uniq od overlapa
        
        unip_func_elem <- c()
        
        if((is.na(ens.seq)) || (ens.seq == "Sequence unavailable")) {
            warning ("No ensembl sequence available for ensembl_gene: ", ensembl_gene, " where the representative transcript is ", ensembl_transcript)
        } else if (is.na (ids.uniprot[1])) {
            # warning ("No UniProtKB identifier available for ensembl_gene: ", ensembl_gene)
        }
        else {
            ens.seq <- gsub("\\*", "", ens.seq)
            ens.seq <- AAString(ens.seq)
            
            for (id.uniprot in ids.uniprot) {
                uni.seq <- as.character(uni.sequences[[id.uniprot]])
                uni.seq <- gsub ("\\*", "", uni.seq)
                uni.seq <- AAString(uni.seq)
                
                names_details <- names (uni.sequences[id.uniprot])
                to_unip_access <- unlist(strsplit(names_details, "\\|"))
                unip_access <- to_unip_access[2]
                
                if(length(ens.seq) & length(uni.seq))
                {
                    #print (ensembl_transcript)
                    #print (unip_access)
                    
                    if(identical(ens.seq, uni.seq)){
                        hotsp.unip <- ensembl_mut_position
                        status <- "ok_match"
                    } else {
                        hotsp.in.unip <- align_to_unip(ens.seq, uni.seq, ensembl_mut_position)
                        status <- hotsp.in.unip[[1]]
                        hotsp.unip <- hotsp.in.unip[[2]]
                    }
                    
                    if (status == "ok_match") {
                        uniprot.data <- gff.file[[unip_access]]
                        if(length(uniprot.data)){
                            type_to_exclude = c ("CHAIN", "SEQUENCE CONFLICT", "NATURAL VARIANT", "VAR_SEQ", "TOPO_DOM", "STRAND", "HELIX", "COILED", "COMPBIAS", "NON_TER", "TURN", "MUTAGENESIS", "ALTERNATIVE SEQUENCE", "MODIFIED RESIDUE")
                            uniprot.data$Type = toupper (uniprot.data$Type)
                            uniprot.data = uniprot.data [!(uniprot.data$Type %in% type_to_exclude), ]
                            overlaping.elements = uniprot.data[(uniprot.data$Start <= hotsp.unip) & (uniprot.data$End >= hotsp.unip),]
                        
                            for (el in seq_len(nrow(overlaping.elements))){
                                el_type <- overlaping.elements[el, "Type"]
                                el_desc_all <- overlaping.elements[el, "Description"]
                                el_desc_spl <- unlist(strsplit(as.character(el_desc_all), ";"))
                                all_notes <- grep ("Note", el_desc_spl)
                                for(w_note in all_notes){
                                    note <- el_desc_spl[w_note]
                                    note <- gsub ("Note=", "", note)
                                    el_type <- tolower (el_type)
                                    func_elem <- paste (el_type, note, sep = ": ")
                                    unip_func_elem <- c(unip_func_elem, func_elem)
                                }
                            }
                        }
                    }
                }
                unip_func_elem = unique (unip_func_elem)
                if(length(unip_func_elem)){
                    hotspot_results[hotspot_results$Ensembl_id == ensembl_gene & hotspot_results$Prot_position == ensembl_mut_position, "Protein_funcional_region"] <- paste(unip_func_elem, collapse ="; ")
                } else {
                    hotspot_results[hotspot_results$Ensembl_id == ensembl_gene & hotspot_results$Prot_position == ensembl_mut_position, "Protein_funcional_region"] <- "NA"
                }
            }
        }
    }
    close(pb)
    
    hotspot_results = hotspot_results[,c(1,2,4,5,14,3,7,8,9,12,13,6)]
    if (write_to_file == "YES"){
      message ("   Saving results in the file Protein_hotspot_residues.txt.")
      write.table (hotspot_results, file = "Protein_hotspot_residues.txt", sep = "\t", row.names = FALSE, col.names = colnames (hotspot_results),quote = FALSE)
    }
    return (hotspot_results)
}
