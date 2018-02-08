calculate_boundary <- function(mut_pos_numeric, length_aa, flanking_region){
    if(length(flanking_region) == 1){
        flanking_region = c(flanking_region, flanking_region)
    }
    if (!(flanking_region[1] < flanking_region[2])) {
        tempp = flanking_region[1]
        flanking_region[1] = flanking_region[2]
        flanking_region[2] = tempp
    }

    # divide flanking region by 2
    flanking_region <- flanking_region / 2

    # calculate boundary positions including flanking base pairs
    min1 <- round (mut_pos_numeric - flanking_region[1])
    plus1 <- round (mut_pos_numeric + flanking_region[1])
    min2 <- round (mut_pos_numeric - flanking_region[2])
    plus2 <- round (mut_pos_numeric + flanking_region[2])

    if ((min1 >= 1) && (plus1 <= length_aa)) {
        to_the_left1 <- min1
        to_the_right1 <- plus1
    } else if (min1 < 1) {
        to_the_left1 <- 1
        rest <- 2 * flanking_region[1] - (mut_pos_numeric - 1)
        to_the_right1 <- mut_pos_numeric + rest
        temp <- to_the_right1
        if (temp > length_aa) {
          to_the_right1 <- length_aa
        }
    } else if (plus1 > length_aa) {
        to_the_right1 <- length_aa
        rest <- 2 * flanking_region[1] - (length_aa - mut_pos_numeric)
        to_the_left1 <- mut_pos_numeric - rest
        temp <- to_the_left1
        if(temp < 1){
            to_the_left1 = 1
        }
    }

    if ((min2 >= 1) && (plus2 <= length_aa)) {
        to_the_left2 <- min2
        to_the_right2 <- plus2
    } else if (min2 < 1) {
        to_the_left2 <- 1
        rest <- 2 * flanking_region[2] - (mut_pos_numeric - 1)
        to_the_right2 <- mut_pos_numeric + rest
        temp <- to_the_right2
        if (temp > length_aa) {
            to_the_right2 <- length_aa
        }
    } else if (plus2 > length_aa) {
        to_the_right2 <- length_aa
        rest <- 2 * flanking_region[2] - (length_aa - mut_pos_numeric)
        to_the_left2 <- mut_pos_numeric - rest
        temp <- to_the_left2
        if(temp < 1){
            to_the_left2 = 1
        }
    }

    return(list(
      region_1 = c(to_the_left1, to_the_right1),
      region_2 = c(to_the_left2, to_the_right2)
    ))
  }
