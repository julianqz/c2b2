
# adapted from removeNtrackGaps from helpers.R for JI 2020

# Depending on `reading_frame` of start position of `sequence`,
# add any necessary padding (denoted using `reserved_char`)
# Input
# - sequence: does NOT need to have reading frame 1
# - reading_frame: of start position of `sequence`; integer {1,2,3}
# - reserved_char: reserved character to be used in padding; 
#                  must not be used in `sequence`
# Output
# - a new string containing the padded version of `sequence`
add_rf_padding = function(sequence, reading_frame, reserved_char) {
    require(seqinr)
    
    # expect reading_frame to be integer ranging [1, 3]
    stopifnot(reading_frame>=1 & reading_frame<=3)
    
    # expect reserved_char to not be in sequence
    stopifnot(!any(reserved_char %in% s2c(sequence)))
    
    # depending on reading frame at 1st position, add any necessary padding
    rf_padding = paste(rep(reserved_char, times=reading_frame-1), 
                       collapse="")
    sequence_rf1 = paste0(rf_padding, sequence)
    
    return(sequence_rf1)
}


# Split in-frame nt sequence into triplets
# IMPORTANT: assumes that `sequence` has reading frame 1
# Input:
# - an in-frame nt sequence
# Output:
# - a vector of triplets. Possible for there to be non-triplet overhang at RHS
get_triplets = function(sequence) {
    require(stringi)
    
    # break into triplets
    # ok if there's non-triplet overhang at the RHS end
    idx_triplet_first_pos = seq(from=1, to=nchar(sequence), by=3)
    
    triplets_gapped = stri_sub(str=sequence, 
                               from=idx_triplet_first_pos, 
                               to=(idx_triplet_first_pos+2))
    
    return(triplets_gapped)
}


# Identify IMGT gaps 
# Important: assumes reading frame 1
# Input: 
# - vec_triplets: a vector of triplets. Possible for there to be non-triplet overhang at RHS
# - def_imgt_gap: definition of an instance of IMGT gap; expected to be 3-char long
# Output:
# - a boolean vector. Triplets that are "..." have `TRUE` values.
identify_gaps = function(vec_triplets, def_imgt_gap="...") {
    
    stopifnot(nchar(def_imgt_gap)==3)
    
    # locate gaps
    bool_triplets_gapped = vec_triplets==def_imgt_gap
    
    return(bool_triplets_gapped)
}


#' Remove IMGT gaps in germline sequence and optionally also 
#' in observed sequence based on gaps present in germline sequence.
#' 
#' An improved version of remove_imgt_gaps() 
#' 
#' @param  germ               IMGT-gapped germline sequence.
#' @param  obsv               IMGT-gapped observed sequence(s); optional.
#' @param  reading_frame      Integer {1..3} specifying the reading frame of the first position.
#' @param  reserved_char      A reserved character not used in `germ` or `obsv` that is
#'                            used in any necessary padding.
#' @param  return_gap_boolean Whether to return a vector of boolean that indicates which of the 
#'                            post-padding triplets are IMGT gaps.
#'                       
#' @return   If return_gap_boolean=FALSE, a vector containing `germ_no_gaps` and `obsv_no_gaps`.
#'           If return_gap_boolean=TRUE, a list containing `germ_no_gaps`, `obsv_no_gaps`,
#'           and a boolean vector indicating gap positions.
#' 
#' @details  IMGT gaps are removed from germline. 
#'           Corresponding gap positions in observed are also removed.
#'           
#'           If `obsv` supplied, assumes:
#'           - `germ` and `obsv` have the same length
#'           - `germ` and `obsv` have the same reading frame
#'           
#'           Only in-frame "..." in `germ` are considered IMGT gaps and thus removed.
#'           
#'           If the intention is to just remove IMGT gaps (regardless of the 
#'           input sequence is germline or not), pass the input sequence to `germ`.
#'        
#' @examples
#' germ = "ABC...EDF..GTGH...JHK.....LOP....W"
#' obsv = "ABCXXXEDFXXGTGHXXXJHKXXXXXLOPXXXXW"
#' remove_imgt_gaps_2(germ, NULL)
#' remove_imgt_gaps_2(obsv, NULL)
#' remove_imgt_gaps_2(germ, obsv)
#' 
remove_imgt_gaps_2 = function(germ, obsv=NULL,
                              reading_frame=1, reserved_char="%",
                              return_gap_boolean=FALSE) {
    
    require(stringi)
    
    ### germ
    
    # add padding based on reading frame
    germ_rf1 = add_rf_padding(sequence=germ, 
                              reading_frame=reading_frame,
                              reserved_char=reserved_char)
    
    # identify gaps
    # wrt padded gapped seq w rf1
    germ_triplets_gapped = get_triplets(germ_rf1)
    germ_bool_triplets_gapped = identify_gaps(germ_triplets_gapped)
    
    if (any(germ_bool_triplets_gapped)) {
        # remove gaps
        germ_no_gaps_rf1 = paste(germ_triplets_gapped[!germ_bool_triplets_gapped], collapse="")
        
        # remove padding 
        germ_no_gaps = substring(text=germ_no_gaps_rf1, first=reading_frame)
        
        ### obsv
        
        if (!is.null(obsv)) {
            
            stopifnot(nchar(germ)==nchar(obsv))
            
            # add padding based on reading frame
            obsv_rf1 = add_rf_padding(sequence=obsv, 
                                      reading_frame=reading_frame,
                                      reserved_char=reserved_char)
            # wrt padded gapped seq w rf1
            obsv_triplets_gapped = get_triplets(obsv_rf1)
            # remove gaps according to where they are in germ
            obsv_no_gaps_rf1 = paste(obsv_triplets_gapped[!germ_bool_triplets_gapped], collapse="")
            # remove padding
            obsv_no_gaps = substring(text=obsv_no_gaps_rf1, first=reading_frame)
            
        } else {
            obsv_no_gaps = obsv
        }
        
    } else {
        germ_no_gaps = germ
        obsv_no_gaps = obsv
    }
    
    if (return_gap_boolean) {
        return_obj = vector(mode="list", length=3)
        names(return_obj) = c("germ_no_gaps", "obsv_no_gaps", "gap_boolean")
        return_obj[["germ_no_gaps"]] = germ_no_gaps
        return_obj[["obsv_no_gaps"]] = obsv_no_gaps
        return_obj[["gap_boolean"]] = germ_bool_triplets_gapped
    } else {
        return_obj = c(germ_no_gaps=germ_no_gaps, obsv_no_gaps=obsv_no_gaps)
    }
    
    return(return_obj)
}


# adapted from restoreGaps from helpers.R for JI 2020

# Add back IMGT gaps to a non-gapped sequence that has had its IMGT gaps removed
# Designed to work hand-in-hand with `remove_imgt_gaps_2`
#
# Input:
# - sequence_no_gap: string of non-gapped sequence
# - gap_boolean: boolean vector indicating which post-padding triplets are gaps
#                in the gapped sequence
# - reading_frame: reading frame of the starting position
# - reserved_character: character not used in `sequence_no_gap` to be used to 
#                       denote padding
# Output:
# - A string of gapped sequence (without padding)
# 
restore_gaps = function(sequence_no_gap, gap_boolean, 
                        reading_frame=1, reserved_char="%") {
    
    # add padding back
    sequence_no_gap_rf1 = add_rf_padding(sequence_no_gap, reading_frame, reserved_char)
    
    # split into triplets
    sequence_no_gap_rf1_triplets = get_triplets(sequence_no_gap_rf1)
    
    # fill in non-gapped positions
    stopifnot(sum(!gap_boolean)==length(sequence_no_gap_rf1_triplets))
    sequence_gapped_rf1_triplets = rep("...", length(gap_boolean))
    sequence_gapped_rf1_triplets[!gap_boolean] = sequence_no_gap_rf1_triplets
    
    # collapse
    sequence_gapped_rf1 = paste(sequence_gapped_rf1_triplets, collapse="")
    
    # remove padding
    sequence_gapped = substring(text=sequence_gapped_rf1, first=reading_frame)
    
    return(sequence_gapped)
}


#### tests ####

library(testthat)

### Should remove gaps as planned

test_that("Remove gaps from germline only; reading frame 1", {
    
    seq_gapped = "ABCD..EF......GHIJ..."
    reading_frame = 1
    
    # not returning boolean vector
    test_obj_1 = remove_imgt_gaps_2(germ=seq_gapped, obsv=NULL, 
                                    reading_frame=reading_frame, 
                                    return_gap_boolean=F)
    
    expect_equivalent(test_obj_1["germ_no_gaps"],
                      expected="ABCD..EF...GHIJ")
    
    # returning boolean vector
    test_obj_2 = remove_imgt_gaps_2(germ=seq_gapped, obsv=NULL, 
                                    reading_frame=reading_frame, 
                                    return_gap_boolean=T)
    
    expect_equivalent(test_obj_2[["gap_boolean"]],
                      expected=c(F,F,F,T,F,F,T))
})


test_that("Remove gaps from germline only; reading frame 2", {
    
    seq_gapped = "ABCD..EF......GHIJ..."
    reading_frame = 2
    
    # not returning boolean vector
    test_obj_1 = remove_imgt_gaps_2(germ=seq_gapped, obsv=NULL, 
                                    reading_frame=reading_frame, 
                                    return_gap_boolean=F)
    
    expect_equivalent(test_obj_1["germ_no_gaps"],
                      expected="ABCD..EFGHIJ...")
    
    # returning boolean vector
    test_obj_2 = remove_imgt_gaps_2(germ=seq_gapped, obsv=NULL, 
                                    reading_frame=reading_frame, 
                                    return_gap_boolean=T)
    
    expect_equivalent(test_obj_2[["gap_boolean"]],
                      expected=c(F,F,F,T,T,F,F,F))
})


test_that("Remove gaps from germline only; reading frame 3", {
    
    seq_gapped = "ABCD..EF......GHIJ..."
    reading_frame = 3
    
    # not returning boolean vector
    test_obj_1 = remove_imgt_gaps_2(germ=seq_gapped, obsv=NULL, 
                                  reading_frame=reading_frame, 
                                  return_gap_boolean=F)
    
    expect_equivalent(test_obj_1["germ_no_gaps"],
                      expected="ABCD..EF...GHIJ...")
    
    # returning boolean vector
    test_obj_2 = remove_imgt_gaps_2(germ=seq_gapped, obsv=NULL, 
                                    reading_frame=reading_frame, 
                                    return_gap_boolean=T)
    
    expect_equivalent(test_obj_2[["gap_boolean"]],
                      expected=c(F,F,F,F,T,F,F,F))
})


test_that("Remove gaps from germline and observed; reading frame 1", {
    
    germ_gapped = "ABCD..EF......GHIJ..."
    obsv_gapped = "123456789012345678901"
    reading_frame = 1
    
    # not returning boolean vector
    test_obj_1 = remove_imgt_gaps_2(germ=germ_gapped, obsv=obsv_gapped, 
                                    reading_frame=reading_frame, 
                                    return_gap_boolean=F)
    
    expect_equivalent(test_obj_1[["obsv_no_gaps"]],
                      expected="123456789345678")
    
    # returning boolean vector
    test_obj_2 = remove_imgt_gaps_2(germ=germ_gapped, obsv=obsv_gapped, 
                                    reading_frame=reading_frame, 
                                    return_gap_boolean=T)
    
    expect_equivalent(test_obj_2[["obsv_no_gaps"]],
                      expected="123456789345678")
    
    expect_equivalent(test_obj_2[["gap_boolean"]],
                      expected=c(F,F,F,T,F,F,T))
})

test_that("Remove gaps from germline and observed; reading frame 2", {
    
    germ_gapped = "ABCD..EF......GHIJ..."
    obsv_gapped = "123456789012345678901"
    reading_frame = 2
    
    # not returning boolean vector
    test_obj_1 = remove_imgt_gaps_2(germ=germ_gapped, obsv=obsv_gapped, 
                                    reading_frame=reading_frame, 
                                    return_gap_boolean=F)
    
    expect_equivalent(test_obj_1[["obsv_no_gaps"]],
                      expected="123456785678901")
    
    # returning boolean vector
    test_obj_2 = remove_imgt_gaps_2(germ=germ_gapped, obsv=obsv_gapped, 
                                    reading_frame=reading_frame, 
                                    return_gap_boolean=T)
    
    expect_equivalent(test_obj_2[["obsv_no_gaps"]],
                      expected="123456785678901")
    
    expect_equivalent(test_obj_2[["gap_boolean"]],
                      expected=c(F,F,F,T,T,F,F,F))
})

test_that("Remove gaps from germline and observed; reading frame 3", {
    
    germ_gapped = "ABCD..EF......GHIJ..."
    obsv_gapped = "123456789012345678901"
    reading_frame = 3
    
    # not returning boolean vector
    test_obj_1 = remove_imgt_gaps_2(germ=germ_gapped, obsv=obsv_gapped, 
                                    reading_frame=reading_frame, 
                                    return_gap_boolean=F)
    
    expect_equivalent(test_obj_1[["obsv_no_gaps"]],
                      expected="123456789045678901")
    
    # returning boolean vector
    test_obj_2 = remove_imgt_gaps_2(germ=germ_gapped, obsv=obsv_gapped, 
                                    reading_frame=reading_frame, 
                                    return_gap_boolean=T)
    
    expect_equivalent(test_obj_2[["obsv_no_gaps"]],
                      expected="123456789045678901")
    
    expect_equivalent(test_obj_2[["gap_boolean"]],
                      expected=c(F,F,F,F,T,F,F,F))
})

### should be able to recover gapped sequence after removing gaps

test_that("Remove and restore gaps; reading frame 1", {
    
    seq_gapped = "ABCD..EF......GHIJ..."
    reading_frame = 1
    
    # returning boolean vector
    test_obj_1a = remove_imgt_gaps_2(germ=seq_gapped, obsv=NULL, 
                                     reading_frame=reading_frame, 
                                     return_gap_boolean=T)
    
    test_obj_1b = restore_gaps(sequence_no_gap=test_obj_1a[["germ_no_gaps"]],
                               gap_boolean=test_obj_1a[["gap_boolean"]],
                               reading_frame=reading_frame)
    
    expect_equivalent(test_obj_1b,
                      expected=seq_gapped)
})

test_that("Remove and restore gaps; reading frame 2", {
    
    seq_gapped = "ABCD..EF......GHIJ..."
    reading_frame = 2
    
    # returning boolean vector
    test_obj_1a = remove_imgt_gaps_2(germ=seq_gapped, obsv=NULL, 
                                     reading_frame=reading_frame, 
                                     return_gap_boolean=T)
    
    test_obj_1b = restore_gaps(sequence_no_gap=test_obj_1a[["germ_no_gaps"]],
                               gap_boolean=test_obj_1a[["gap_boolean"]],
                               reading_frame=reading_frame)
    
    expect_equivalent(test_obj_1b,
                      expected=seq_gapped)
})


test_that("Remove and restore gaps; reading frame 3", {
    
    seq_gapped = "ABCD..EF......GHIJ..."
    reading_frame = 3
    
    # returning boolean vector
    test_obj_1a = remove_imgt_gaps_2(germ=seq_gapped, obsv=NULL, 
                                     reading_frame=reading_frame, 
                                     return_gap_boolean=T)
    
    test_obj_1b = restore_gaps(sequence_no_gap=test_obj_1a[["germ_no_gaps"]],
                               gap_boolean=test_obj_1a[["gap_boolean"]],
                               reading_frame=reading_frame)
    
    expect_equivalent(test_obj_1b,
                      expected=seq_gapped)
})

