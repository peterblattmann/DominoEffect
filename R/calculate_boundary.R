calculate_boundary <- function(mut_pos_numeric, length_aa, flanking_region){
    ### Checking for argument requirements
    if(length(mut_pos_numeric) != 1L){
        stop("mut_pos_numeric must be an integer.")
    }
    
    if(length(length_aa) != 1L){
        stop("length_aa must be an integer.")
    }
    
    if(length(flanking_region) != 2L){
        if(length(flanking_region) == 1L){
            flanking_region <- c(flanking_region, flanking_region)
        } else {stop("flanking_region must be one or two integers")}
    }
    ###
    
    # sort and divide flanking region by 2
    flanking_region <- sort(flanking_region)/2
    
    # calculate boundary positions including flanking base pairs
    min <- round(mut_pos_numeric - flanking_region)
    plus <- round(mut_pos_numeric + flanking_region)
    left_boundary <- right_boundary <- rep(NA, 2)

    for(i in seq_len(length(flanking_region))){
        # if region is in the middle of sequence it does not need to shift
        if ((min[i] >= 1) && (plus[i] <= length_aa)) {
            left_boundary[i] <- min[i]
            right_boundary[i] <- plus[i]
        } 
        
        #if region is at the beginning, start region at 1
        else if (min[i] < 1) {
            left_boundary[i] <- 1
            right_boundary[i] <- 1 + 2*flanking_region[i]
            if(right_boundary[i] > length_aa){right_boundary[i] <- length_aa}
        } 
        
        #if region is at the end, end region at maximum position
        else if (plus[i] > length_aa) {
            right_boundary[i] <- length_aa
            left_boundary[i] <-  length_aa - 2*flanking_region[i]  
            if(left_boundary[i] < 1){left_boundary[i] <- 1}
        }
        
    }
 
    return(list(region_1 = c(left_boundary[1], right_boundary[1]),
                region_2 = c(left_boundary[2], right_boundary[2])))
  }
