
align_to_unip <- function(ens.seq, uni.seq, ensembl_mut_position){
    window.size = 15
    flank_size <- as.integer((window.size-1)/2)
    
    # select sub-sequence from specified window size
    if(ensembl_mut_position > flank_size & ensembl_mut_position < nchar(ens.seq)- flank_size){
        ens.subseq <- subseq(ens.seq, ensembl_mut_position - flank_size, ensembl_mut_position + flank_size)
        residue_shift <- 7
    } else {
        #if mutation is too near to end or start use the specified window size from start or end
        if(ensembl_mut_position <= flank_size){
            ens.subseq <- subseq(ens.seq, 1, window.size)
            residue_shift <- ensembl_mut_position - 1
        } else {
            ens.subseq <- subseq(ens.seq, nchar(ens.seq)-window.size, nchar(ens.seq))
            residue_shift <- window.size - (nchar(ens.seq)-ensembl_mut_position) - 1
        }
    }
    
    ### ALIGNMENT FUNCTION CAN STILL BE IMPROVE TO INCLUDE ADDITIONAL NON-PERFECT MATCHES
    alignment.subseq <- pairwiseAlignment(subject = uni.seq, pattern = ens.subseq, 
                                          substitutionMatrix = "BLOSUM50", type = "local")
    
    n.match <- countPattern(pattern = ens.subseq, subject = uni.seq)
    seq.match <- matchPattern(pattern = ens.subseq, subject = uni.seq)
    
    aa.match <- "go"
    uniprot_mut_position <- "NA"
    if(!(n.match == 0)) {
        best.match = max(width(seq.match))
        if(best.match == 15) {
            count_i <- 1
            all_starts = start (seq.match)
            for (loc_start in all_starts)
            {
                pattern_width = width(seq.match[count_i])
                if ((aa.match == "go") & (pattern_width == best.match)) {
                    identity.subseq <- pid(alignment.subseq[count_i], type = "PID3")
                    if ((identity.subseq/100*best.match) >= 13) {
                        ### problem with gaps is when the begining doesn't start there
                        ## when we address that we can propery correct the residue shift
                        ## residue_shift <- residue_shift + gaps.ens - gaps.uni
                        gaps.ens <- grep("-", pattern(alignment.subseq[1]))
                        gaps.uni <- grep("-", subject(alignment.subseq[1]))
                        
                        if ((is.na (gaps.ens[1])) & (is.na (gaps.uni[1]))) {
                            uniprot_mut_position = loc_start + residue_shift
                            ens.aa <- subseq(ens.seq, ensembl_mut_position, ensembl_mut_position)
                            uni.aa <- subseq(uni.seq, uniprot_mut_position, uniprot_mut_position)
                            if(ens.aa == uni.aa) {
                                aa.match = "ok_match"
                            }
                        }
                    }
                }
                count_i = count_i + 1
            }
        }
    }
    
    return (list (aa.match, uniprot_mut_position))
}