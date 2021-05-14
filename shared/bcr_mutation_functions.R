#' Calculate AA replacement mutation for a single nucleotide sequence
#' 
#' @param  germ  Germline nucleotide sequence.
#' @param  obsv  Observed nucleotide sequence.
#' @param  upper_bound  Upper bound nt position (inclusive) of sequence to use.
#'                      Either `NULL` or a whole number.
#' 
#' @return   A vector, `c(count, freq, denom)`, where
#'           - count is the number of AA (non-synonymous) mutation(s)
#'           - freq is the frequency of AA (non-synonymous) mutation(s)
#'           - denom is the total number of AA position(s) considered
#'           - if denom is 0, count and freq will be NA; otherwise, 
#'             count and freq will be numeric (incl 0)
#' 
calc_aa_mutation = function(germ, obsv, upper_bound=NULL) {
    
    suppressPackageStartupMessages(require(alakazam))
    suppressPackageStartupMessages(require(seqinr))
    suppressPackageStartupMessages(require(stringi))
    
    germ_len = stri_length(germ)
    obsv_len = stri_length(obsv)
    
    if (germ_len<3 | obsv_len<3) {
        # NA if no valid position exists for calculation
        return(c(count=NA, freq=NA, denom=0))
    }
    
    if (is.null(upper_bound)) {
        # if NULL, use up to common_len
        
        # ok if not a multiple of 3
        # alakazam::translateDNA will ignore any non-triplet overhang
        common_len = min(germ_len, obsv_len)
        
        if (common_len<3) {
            # NA if no valid position exists for calculation
            return(c(count=NA, freq=NA, denom=0))
        }
        
        # use up to common_len
        germ_2 = translateDNA(seq=stri_sub(germ, from=1, to=common_len), trim=F)
        obsv_2 = translateDNA(seq=stri_sub(obsv, from=1, to=common_len), trim=F)
        
    } else{
        # if not NULL, expect a whole number
        
        if (upper_bound<3) {
            # germ_2/obsv_2 would be NA below
            # NA if no valid position exists for calculation
            return(c(count=NA, freq=NA, denom=0))
        }
        
        # use up to IMGT V 1..upper_bound
        germ_2 = translateDNA(seq=stri_sub(germ, from=1, to=min(germ_len, upper_bound)), trim=F)
        obsv_2 = translateDNA(seq=stri_sub(obsv, from=1, to=min(obsv_len, upper_bound)), trim=F)
        
        germ_2 = stri_sub(germ_2, from=1, to=min(stri_length(germ_2), stri_length(obsv_2)))
        obsv_2 = stri_sub(obsv_2, from=1, to=min(stri_length(germ_2), stri_length(obsv_2)))
    }
    
    stopifnot(stri_length(germ_2)==stri_length(obsv_2))
    
    germ_2_c = s2c(germ_2)
    obsv_2_c = s2c(obsv_2)
    
    # only consider positions where neither germline nor observed is X
    bool = (germ_2_c!="X") & (obsv_2_c!="X")
    
    if (any(bool)) {
        denom = sum(bool)
        germ_2_c_sub = germ_2_c[bool]
        obsv_2_c_sub = obsv_2_c[bool]
        count = sum(germ_2_c_sub!=obsv_2_c_sub)
        freq = count/denom
        return(c(count=count, freq=freq, denom=denom))
    } else {
        # NA if no valid position exists for calculation
        return(c(count=NA, freq=NA, denom=0))
    }
    
}


#' Parse output of `shazam::calcObservedMutations` with `returnRaw=TRUE`
#' 
#' @param   obj  Output from a call to `shazam::calcObservedMutations` with 
#'               `returnRaw=TRUE`. Either a `data.frame` or `NA`.
#' 
#' @return  A vector `c(nuc_denom, nuc_R, nuc_S, nuc_RS, nuc_RS_freq)`.
#' 
parse_shazam_mutation = function(obj) {

    # obj is a list, with components "pos" and "nonN"
    # nonN is a whole number
    # pos is either a data.frame, or NA
    #     if data.frame, pos has columns "position", "R" or "r", "S" or "s", "region"
    
    # do not use is.na(obj$pos) -- will return multiple booleans if obj$pos is data.frame
    if (is.data.frame(obj$pos)) {
        # depending on shazam version, colnames could be either {R,S} or {r,s}
        if ("R" %in% colnames(obj$pos)) {
            col_r = "R"
            col_s = "S"
        } else {
            col_r = "r"
            col_s = "s"
        }
        nuc_R = sum(obj$pos[[col_r]])
        nuc_S = sum(obj$pos[[col_s]])
        nuc_RS = nuc_R+nuc_S
    } else {
        nuc_R = 0
        nuc_S = 0
        nuc_RS = 0
    }
    
    nuc_denom = obj$nonN
    # otherwise would mess up the name in the returned vector
    names(nuc_denom)=NULL
    
    if (nuc_denom!=0) {
        nuc_RS_freq = nuc_RS / nuc_denom
    } else {
        nuc_RS_freq = NA
    }
    
    return( c( nuc_denom=nuc_denom, nuc_R=nuc_R, nuc_S=nuc_S, nuc_RS=nuc_RS, 
               nuc_RS_freq=nuc_RS_freq ) )
}

