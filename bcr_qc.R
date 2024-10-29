# Julian Q. Zhou
# https://github.com/julianqz

# functions to perform QC on BCRs

#' Check chain consistency of V, D, J, C gene annotations of a sequence
#' 
#' @param    vg          V gene annotation(s) of a sequence. 
#' @param    dg          D gene annotation(s) of a sequence
#' @param    jg          J gene annotation(s) of a sequence
#' @param    cg          C gene annotation(s) of a sequence. 
#' @param    loci        One of "IG" or "TR".
#' @param    separator   Character that separates multiple annotations. 
#'                       Default to `,` (comma).
#' @param    verbose     Boolean. If TRUE, prints the inconsistency when check 
#'                       is failed.
#' 
#' @returns  TRUE or FALSE.
#' 
#' @details  Checks if all non-NULL, non-NA, non-empty input annotations have 
#'           the same chain type, e.g. all "IGH", "IGK", "IGL", "TRA", etc.

inspect_chain_consistency = function(vg, dg, jg, cg, 
                                     loci=c("IG", "TR"), 
                                     separator=",",
                                     verbose=F){
    require(stringi)
    
    # check value for `loci` is valid 
    stopifnot(loci %in% c("IG", "TR"))
    
    # if any is NULL, would be dropped by c()
    chains = c(vg, dg, jg, cg)
    # at least one should be non-NULL
    stopifnot(length(chains)>=1)
    
    # at least one annotation must be non-NA and non-empty 
    # TRUE if non-NA and non-empty
    bool_ok = !is.na(chains) & chains!="" 
    stopifnot( sum(bool_ok)>= 1 )
    
    # NA and empty annotations ignored
    chains = chains[bool_ok]
    
    # assumed format:
    # IG[HKL]*
    # TR[ABGD]*
    
    # convert all to uppercase
    chains = toupper(chains)
    
    # unpack multiple annotations (if any)
    # first split by ","
    # then get IG/TR*
    # then extract first 3 chars (IG[HKL], TR[ABGD])
    chains_unpacked = sapply(chains, function(s){
        s_split = strsplit(s, split=separator)[[1]]
        s_split_extr = sapply(s_split, function(ss) {
            stri_extract_first(str=ss,
                               regex="[IT][GR][HKLABGD][[:alnum:]]?")
        }, USE.NAMES=F)
        return( substr(s_split_extr, 1, 3) )
    }, simplify=F)
    
    # unique chain types present
    uniq_chain = unique(unlist(chains_unpacked))
    
    # sanity check
    # at least 1
    stopifnot( length(uniq_chain)>=1 )
    # all starts with IG* or TR*
    stopifnot( all(sapply(uniq_chain, grepl, pattern=paste0("^", loci))) )
    
    # consistent vs. inconsistent
    if (length(uniq_chain)==1) {
        return(T)
    } else {
        if (verbose) { cat("Failed chain consistency check:", uniq_chain, "\n") }
        return(F)
    }
}


#' Perform sequence-level QC for BCR sequences
#'
#' @param   db                          data.frame
#' @param   chain_type                  One of "IG" or "TR.
#' @param   col_v_call                  Column name for V call. Required.
#' @param   col_j_call                  Column name for J call. Required.
#' @param   col_d_call                  Column name for D call. Can be `NA`.
#' @param   col_c_call                  Column name for C call. Can be `NA`.
#' @param   check_valid_vj              Boolean. Whether to check that both 
#'                                      V and J gene annotations are valid & non-empty.
#' @param   check_chain_consistency     Boolean. Whether to check that gene
#'                                      annotations have chain consistency.
#' @param   check_N                     Boolean. Whether to check the # or % of 
#'                                      N's in `col_N`.
#' @param   max_N                       Max # or % of N's in `col_N` allowed from
#'                                      position 1 thru `last_pos_N`.
#' @param   col_N                       Column name(s) in which to perform `check_N`.
#' @param   last_pos_N                  Last position(s) in `col_N` thru which
#'                                      to perform `check_N`. Length should 
#'                                      match that of `col_N`.
#' @param   as_perc_N                   Boolean. If `TRUE`, check the % of N's.
#'                                      If `FALSE`, check the number of N's.                                    
#' @param   check_nonATGC               Boolean. Whether to check the number of 
#'                                      non-ATGC positions in `col_obsv` using
#'                                      `col_germ` as reference.
#' @param   max_nonATGC                 Max number of non-ATGC positions allowed
#'                                      in `col_obsv` from position 1 thru 
#'                                      `last_pos_nonATGC`.
#' @param   last_pos_nonATGC            The last position thru `col_obsv` to perform
#'                                      `check_nonATGC`.
#' @param   as_perc_nonATGC             Boolean. If `TRUE`, check the % of non-ATGC
#'                                      positions. If `FALSE`, check the number.                                     
#' @param   check_none_empty            Boolean. Whether to check for `[Nn]one`
#'                                      and empty (`""`) values.
#' @param   col_none_empty              Column name(s) in which to perform 
#'                                      `check_none_empty`. Note that such column(s)
#'                                      should be of class `character`.                                      
#' @param   check_NA                    Boolean. Whether to check for `NA`.
#' @param   col_NA                      Column name(s) in which to perform `check_NA`.
#' @param   check_len_mod3              Boolean. Whether to check that lengths
#'                                      are a multiple of 3.
#' @param   col_len_mod3                Column name(s) in which to perform 
#'                                      `check_len_mod3`. Note that column(s) containing
#'                                      strings is/are expected. Lengths will be calculated
#'                                      on the fly.
#' 
#' @returns A bool vector of length `nrow(db)` indicating whether each row
#'          passed all checks of choice.
#'
#' @details - `check_valid_vj`: both V and J start with "IG" or "TR"
#'
#'          - `check_chain_consistency`: all of V, D, J, and C calls, 
#'            where supplied, are "IGH", "IGK", "IGL", or "TRA", etc.
#'            Ok to one or more, but not all, of `col_[vdjc]_call` as `NA`.
#'            Also see `inspect_chain_consistency`.
#'
#'          - `check_N` checks if the # or % of positions that are N in 
#'            `col_N` from position 1 thru `last_pos_N` is `<=` (at most)
#'            `max_N`. This check is skipped for a row which has
#'            `""`, `"[Nn]one"`, or `NA` for `col_germ`.
#'
#'          - `check_nonATGC` checks if the number of positions that are
#'            non-ATGC in `col_obsv`, between position 1 thru `last_pos_nonATGC`,
#'            at positions that are ATGC in `col_germ`, is `<=` (at most)
#'            `max_nonATGC`. This check is skipped for a row which has
#'            `""`, `"[Nn]one"`, or `NA` for `col_germ`.
#'            
#'          - When `as_perc_N` or `as_perc_nonATGC` is `TRUE`, `max_N` or 
#'            `max_nonATGC` should be in the range of [0, 100]. A value of 10 
#'            would mean 10%.
#'
#'          - `check_none_empty` checks if `col_none_empty` is `[Nn]one` or `""`.
#'            This check is for column(s) of class `character`.
#'            This check is skipped for a row which has `NA`. 
#'
#'          - `check_NA` checks if `col_NA` is `NA`.
#'
#'          - `check_len_mod3` check if `nchar(col_len_mod3) %% 3 == 0`.
#'            This check is skipped for a row which has
#'            `""`, `"[Nn]one"`, or `NA` for `col_germ`.
#'
#'            
#'          After all checks of choice are performed, an `AND` operation is 
#'          performed across the check results for each row.
                                                                                                                                                                                                                                                                                                                                                  
perform_qc_seq = function(db, chain_type=c("IG", "TR"),
                          col_v_call, col_d_call, col_j_call, col_c_call,
                          check_valid_vj=F, 
                          check_chain_consistency=F, 
                          check_N=F, max_N, col_N, last_pos_N, as_perc_N,
                          check_nonATGC=F, col_obsv, col_germ,
                          max_nonATGC, last_pos_nonATGC, as_perc_nonATGC,
                          check_none_empty=F, col_none_empty,
                          check_NA=F, col_NA,
                          check_len_mod3=F, col_len_mod3) {
    
    require(stringi)
    require(seqinr)
    
    if (check_N) {
        stopifnot( length(col_N) == length(last_pos_N) )
        if (as_perc_N) {
            stopifnot( max_N>=0 & max_N<=100 )
        } else {
            stopifnot( max_N>=0 )
        }
    }
    
    if (check_nonATGC) {
        if (as_perc_nonATGC) {
            stopifnot( max_nonATGC>=0 & max_nonATGC<=100 )
        } else {
            stopifnot( max_nonATGC>=0 )
        }
    } 
    
    nseqs = nrow(db)
    cat("\nPerforming sequence-level QC on", nseqs, "sequences\n")
    
    if (check_valid_vj) {
        # IgBLAST V/D/J calls: IG[HKL][VDJ]... or TR[ABDG][VDJ]...
        pattern_valid_1 = paste0("^", chain_type)

        # IMGT High/V-QUEST calls: "Musmus IGHV1-81*01 F"
        # space before IG/TR
        # \s: space
        pattern_valid_2 = paste0("\\s", chain_type)
        
        # TRUE means valid 
        bool_valid_v = grepl(pattern=pattern_valid_1, x=db[[col_v_call]]) | grepl(pattern=pattern_valid_2, x=db[[col_v_call]])
        bool_valid_j = grepl(pattern=pattern_valid_1, x=db[[col_j_call]]) | grepl(pattern=pattern_valid_2, x=db[[col_j_call]])
        
        # missing V
        table(db[[col_v_call]][!bool_valid_v], useNA="ifany")
        
        # missing J
        table(db[[col_j_call]][!bool_valid_j], useNA="ifany")
        
        # valid v and j
        bool_valid_vj = bool_valid_v & bool_valid_j
        
        # count
        cat("\ncheck_valid_vj:\n")
        cat("cols:", col_v_call, col_j_call, "\n")
        print(table(bool_valid_vj, useNA="ifany"))
        cat("\n")
        
    } else {
        bool_valid_vj = rep(T, nseqs)
    }
    
    
    if (check_chain_consistency) {
        
        # if col_[vdjc]_call does not exist, db[[ ]][i] evaluates to NULL, 
        #    which can be handled by inspect_chain_consistency
        
        # whether columns are supplied and exist
        bool_v = !is.na(col_v_call) && col_v_call %in% colnames(db)
        bool_d = !is.na(col_d_call) && col_d_call %in% colnames(db)
        bool_j = !is.na(col_j_call) && col_j_call %in% colnames(db)
        bool_c = !is.na(col_c_call) && col_c_call %in% colnames(db)
        
        # TRUE means consistent
        bool_chain = sapply(1:nrow(db), function(i){
            # DO NOT use `NULL` inside ifelse; will cause error; use `NA` instead
            inspect_chain_consistency(vg=ifelse(bool_v, db[[col_v_call]][i], NA),
                                      dg=ifelse(bool_d, db[[col_d_call]][i], NA),
                                      jg=ifelse(bool_j, db[[col_j_call]][i], NA),
                                      cg=ifelse(bool_c, db[[col_c_call]][i], NA),
                                      loci=chain_type,
                                      verbose=F)
        })
        
        # count
        cat("\ncheck_chain_consistency:\n")
        cat("cols:", col_v_call, col_d_call, col_j_call, col_c_call, "\n")
        print(table(bool_chain, useNA="ifany"))
        cat("\n")
        
    } else {
        bool_chain = rep(T, nseqs)
    }
    
    
    if (check_N) {
        
        # the next stopifnot will fail if column does not exist in db
        # remove such cols first
        col_N = col_N[col_N %in% colnames(db)]
        stopifnot(length(col_N)>=1)
        
        # check that all columns to be checked are characters
        stopifnot(all( sapply(col_N, 
                              function(s){ class(db[[s]]) }) == "character" ))
        
        # matrix, even when col_N is of length 1
        # each col is a col in col_pN
        # each row is a seq in db
        # TRUE means #/% of N's in pos 1 thru last_pos_N <= max_N
        mtx_N = do.call(cbind, 
                             sapply(1:length(col_N), function(i){
                                 s = col_N[i]
                                 p = last_pos_N[i]
                                 cur_s = db[[s]]
                                 
                                 # convert to uppercase so no need to deal with cases
                                 cur_s = toupper(cur_s)
                                 
                                 # skip if NA, "", "[Nn]one"
                                 bool_skip = is.na(cur_s) | cur_s=="" | cur_s=="NONE"
                                 
                                 # truncate to pos 1 to last_pos_N
                                 idx_truncate = which( nchar(cur_s) > p )
                                 if (length(idx_truncate)>0) {
                                     cur_s[idx_truncate] = sapply(idx_truncate, function(i_db){
                                         return(substr(cur_s[i_db], 1, p))
                                     }, USE.NAMES=F)
                                 }
                                 
                                 # count number of occurrences of N characters
                                 cur_s_count = stri_count_fixed(str=cur_s,
                                                                pattern="N")
                                 # calc %
                                 if (as_perc_N) {
                                     cur_s_count = cur_s_count / nchar(cur_s) * 100   
                                 }
                                 
                                 return( (cur_s_count <= max_N) | bool_skip )
                             }, simplify=F))
        # sum across columns
        # if all TRUE, rowSums should be the same as number of cols
        bool_N = rowSums(mtx_N, na.rm=T) == ncol(mtx_N)
        
        # count
        cat("\ncheck_N ( max_N=", max_N, 
            ifelse(as_perc_N, "%", ""), "; <=):\n")
        cat("cols:", col_N, "\n")
        print(table(bool_N, useNA="ifany"))
        cat("\n")
        
    } else {
        bool_N = rep(T, nseqs)
    }
    
    
    if (check_nonATGC) {
        
        # check that col_obsv and col_germ exist
        stopifnot(all(c(col_obsv, col_germ) %in% colnames(db)))
        
        bool_skip = is.na(db[[col_germ]]) | db[[col_germ]]=="" | tolower(db[[col_germ]])=="none"
        
        # convert to uppercase so no need to deal with cases
        # truncate to pos 1 to last_pos_nonATGC
        germ_upper_trunc = sapply(db[[col_germ]], function(s){
            toupper( substr(s, 1, min(nchar(s), last_pos_nonATGC) ) )
        }, USE.NAMES=F)
        # same for observed
        obsv_upper_trunc = sapply(db[[col_obsv]], function(s){
            toupper( substr(s, 1, min(nchar(s), last_pos_nonATGC) ) )
        }, USE.NAMES=F)
        
        # positions in germline from pos 1 thru last_pos_nonATGC
        # that are non-ATGC (e.g. ".")
        # this should be a list
        germ_idx = sapply(germ_upper_trunc, function(s){
            return(which(!s2c(s) %in% c("A","T","G","C")))
        }, simplify=F)
        
        # count the number of non-ATGCs in obsv, excl non-ATGC positions in germ
        obsv_count = sapply(1:nrow(db), function(i){
            cur_idx_excl = germ_idx[[i]]
            
            if (length(cur_idx_excl)>0) {
                # exclude positions that are non-ATGC in germ from obsv
                cur_obsv = c2s( s2c(obsv_upper_trunc[i])[-cur_idx_excl] )
            } else {
                cur_obsv = obsv_upper_trunc[i]
            }
            
            cur_obsv_count = stri_count_regex(str=cur_obsv, pattern="[^ATGC]")
            
            if (as_perc_nonATGC) {
                cur_obsv_count = cur_obsv_count / nchar(cur_obsv) *100
            }
            
            names(cur_obsv_count) = NULL
            
            return(cur_obsv_count)
        })
        
        # TRUE means # of non-ATGCs in pos 1 thru last_pos_nonATGC <= max_nonATGC
        # skip check if germline missing (and set to TRUE for that row)
        bool_nonATGC = (obsv_count <= max_nonATGC) | bool_skip
        
        # count
        cat("\ncheck_nonATGC ( max_nonATGC=", max_nonATGC, 
            ifelse(as_perc_nonATGC, "%", ""), "; <=):\n")
        cat("observed col:", col_obsv, "; germline col:", col_germ, "\n")
        print(table(bool_nonATGC, useNA="ifany"))
        cat("\n")
         
        
    } else {
        bool_nonATGC = rep(T, nseqs)
    }
    
    
    if (check_none_empty) {
        
        # the next stopifnot will fail if column does not exist in db
        # remove such cols first
        col_none_empty = col_none_empty[col_none_empty %in% colnames(db)]
        stopifnot(length(col_none_empty)>=1)
        
        # check that all columns to be checked are characters
        stopifnot(all( sapply(col_none_empty, 
                              function(s){ class(db[[s]]) }) == "character" ))
        
        # matrix, even when col_none_empty is of length 1
        # each col is a col in col_none_empty
        # each row is a seq in db
        # TRUE means not none AND not empty (row with NA skipped)
        mtx_none_empty = do.call(cbind, 
                                 sapply(col_none_empty, function(s){
                                     # convert to lowercase so no need to deal with cases
                                     return( (tolower(db[[s]]) != "none" & db[[s]] != "") | is.na(db[[s]]) )
                                 }, simplify=F))
        # sum across columns
        # if all TRUE, rowSums should be the same as number of cols
        bool_none_empty = rowSums(mtx_none_empty, na.rm=T) == ncol(mtx_none_empty)
        
        # count
        cat("\ncheck_none_empty:\n")
        cat("cols:", col_none_empty, "\n")
        print(table(bool_none_empty, useNA="ifany"))
        cat("\n")
        
        
    } else {
        bool_none_empty = rep(T, nseqs)
    }
    
    
    if (check_NA) {
        
        # remove non-existing columns
        col_NA = col_NA[col_NA %in% colnames(db)]
        stopifnot(length(col_NA)>=1)
        
        # matrix, even when col_NA is of length 1
        # each col is a col in col_NA
        # each row is a seq in db
        # TRUE means not NA
        mtx_NA = do.call(cbind, 
                         sapply(col_NA, function(s){
                             return( !is.na(db[[s]]) )
                         }, simplify=F))
        # sum across columns
        # if all TRUE, rowSums should be the same as number of cols
        bool_NA = rowSums(mtx_NA, na.rm=T) == ncol(mtx_NA)
        
        # count
        cat("\ncheck_NA:\n")
        cat("cols:", col_NA, "\n")
        print(table(bool_NA, useNA="ifany"))
        cat("\n")
        
    } else {
        bool_NA = rep(T, nseqs)
    }
    
    
    if (check_len_mod3) {
        
        # the next stopifnot will fail if column does not exist in db
        # remove such cols first
        col_len_mod3 = col_len_mod3[col_len_mod3 %in% colnames(db)]
        stopifnot(length(col_len_mod3)>=1)
        
        # check that all columns to be checked are characters
        stopifnot(all( sapply(col_len_mod3, 
                              function(s){ class(db[[s]]) }) == "character" ))
        
        # matrix, even when col_len_mod3 is of length 1
        # each col is a col in col_len_mod3
        # each row is a seq in db
        # TRUE means length is a multipel of 3 (row with NA/empty/[Nn]one skipped)
        mtx_len_mod3 = do.call(cbind, 
                               sapply(col_len_mod3, function(s){
                                   bool_skip = is.na(db[[s]]) | db[[s]]=="" | tolower(db[[s]])=="none"
                                   return( (nchar(db[[s]]) %% 3 == 0) | bool_skip )
                               }, simplify=F))
        # sum across columns
        # if all TRUE, rowSums should be the same as number of cols
        bool_len_mod3 = rowSums(mtx_len_mod3, na.rm=T) == ncol(mtx_len_mod3)
        
        # count
        cat("\ncheck_len_mod3:\n")
        cat("cols:", col_len_mod3, "\n")
        print(table(bool_len_mod3, useNA="ifany"))
        cat("\n")
        
    } else {
        bool_len_mod3 = rep(T, nseqs)
    }
    
    
    # combine
    stopifnot(!any(is.na(bool_valid_vj)))
    stopifnot(!any(is.na(bool_chain)))
    stopifnot(!any(is.na(bool_N)))
    stopifnot(!any(is.na(bool_nonATGC)))
    stopifnot(!any(is.na(bool_none_empty)))
    stopifnot(!any(is.na(bool_NA)))
    stopifnot(!any(is.na(bool_len_mod3)))

    bool = bool_valid_vj & bool_chain & bool_N & bool_nonATGC & bool_none_empty & bool_NA & bool_len_mod3
    
    # count
    cat("\nAfter seq-level QC, number of seqs:\n")
    print(table(bool, useNA="ifany"))
    cat("\n")
    
    return(bool)
}


#' Perform cell-level QC for BCR sequences
#' 
#' @param   db                          data.frame
#' @param   chain_type                  One of "IG" or "TR.
#' @param   col_locus                   Column name for locus. Expected values
#'                                      for `chain_type="IG"` are `{IGH, IGK, IGL}`.
#' @param   col_cell                    Column name for cell ID.
#' @param   col_umi                     Column name for UMI count.                                     
#' @param   check_locus                 Boolean. Whether to perform check on the 
#'                                      consistency between V call and locus annotation.
#' @param   col_v_call                  Column name for V call.
#' @param   check_num_HL                Boolean. Whether to perform check on the 
#'                                      number of heavy and light chains per cell.
#' @param   logic_num_HL                The logic to be applied to the check on the
#'                                      number of heavy and light chains per cell.
#'                                      One of `1H_1L`, `1H_min1L`, 
#'                                      `1H_min1L_or_min1H_1L`, or `1H`,
#'
#' @returns A bool vector of length `nrow(db)` indicating whether each row
#'          passed all checks of choice.
#'          
#' @details While each row in the `db` supplied as input to the function is 
#'          expected to represent a sequence, the cell-level QC is applied 
#'          on a cell-by-cell basis, and considers all sequences linked to a cell
#'          during the check for that cell. The result of the check is then 
#'          propagated to all sequences linked to that cell.    
#'          
#'          `1H_1L`: pass if a cell has exactly 1 heavy chain and exactly 1 light chain
#'          `1H_min1L`: pass if a cell has exactly 1 heavy chain and at least 1 light chain
#'          `1H_min1L_or_min1H_1L`: pass if 
#'          - either a cell has exactly 1 heavy chain and at least 1 light chain
#'          - or a cell has at least 1 heavy chain and exactly 1 light chain
#'          `1H`: pass if a cell has exactly 1 heavy chain    
#'          
#'          In the cases of `1H_min1L` and `1H_min1L_or_min1H_1L`, within each 
#'          cell, for the chain type with more than 1 sequence, the most abundant
#'          sequence (in terms of UMI count) is kept while the rest discarded. 
#'          When tied, the sequence that appears earlier in the db is kept.
#'          In other words, exactly 1 heavy and 1 light per cell is kept ultimately.
#'                                                                                                                                                                          
perform_qc_cell = function(db, chain_type=c("IG", "TR"), 
                           col_locus, col_cell, col_umi,
                           check_locus, col_v_call, 
                           check_num_HL, 
                           logic_num_HL=c("1H_1L", "1H",
                                          "1H_min1L", "1H_min1L_or_min1H_1L")
                           ) {
    
    stopifnot( all(c(col_locus, col_cell) %in% colnames(db)) )
    if (check_locus) { stopifnot( col_v_call %in% colnames(db) ) }
    
    uniq_cells = unique(db[[col_cell]])
    cat("\nNumber of unique cells before cell-level QC:", 
        length(uniq_cells), "\n")
    
    idx_uniq_cells = match(uniq_cells, db[[col_cell]])
    stopifnot( all.equal(uniq_cells, db[[col_cell]][idx_uniq_cells]) )
    
    # validate locus
    if (check_locus) {
        # based on V call
        vec_chain = substr(db[[col_v_call]], 1, 3)
        
        #* currently only supports chain_type="IG"
        #  expect values {IGH, IGK, IGL}
        stopifnot( all.equal(vec_chain, db[[col_locus]]) )
    }
    
    # number of heavy and light chain(s) per cell
    if (check_num_HL) {
        #* currently only supports chain_type="IG"
        chain_count_mtx = matrix(NA, nrow=length(uniq_cells), ncol=2)
        colnames(chain_count_mtx) = c("heavy", "light")
        rownames(chain_count_mtx) = uniq_cells
        
        # row: IG[HKL]
        # col: cell ID
        tab_cell_chain = table(db[[col_locus]], db[[col_cell]], useNA="ifany")
        
        # wrt tab_cell_chain
        idx_tab = match(uniq_cells, colnames(tab_cell_chain))
        stopifnot(all.equal(uniq_cells, colnames(tab_cell_chain)[idx_tab]))
        
        chain_count_mtx[, "heavy"] = tab_cell_chain["IGH", idx_tab]
        chain_count_mtx[, "light"] = colSums(tab_cell_chain[c("IGL", "IGK"), idx_tab])
        
        # sanity check
        stopifnot( all.equal( rowSums(chain_count_mtx), 
                              colSums(tab_cell_chain)[idx_tab], 
                              check.attributes=F ) )
        
        cat("\nNumber of heavy chains per cell:\n")
        print(table(chain_count_mtx[, "heavy"], useNA="ifany"))
        cat("\n")
        
        cat("\nNumber of light chains per cell:\n")
        print(table(chain_count_mtx[, "light"], useNA="ifany"))
        cat("\n")
        
        cat("\nNumber of heavy and light chains per cell:\n")
        print(table(chain_count_mtx[, "heavy"], chain_count_mtx[, "light"]))
        cat("\n")
        
        cat("\nConfig:", logic_num_HL, "\n")
        
        if (logic_num_HL=="1H_1L") {
            # exactly 1 heavy, exactly 1 light
            bool_num_HL = chain_count_mtx[, "heavy"]==1 & chain_count_mtx[, "light"]==1
        } else if (logic_num_HL=="1H_min1L") {
            # exactly 1 heavy, at least 1 light
            bool_num_HL = chain_count_mtx[, "heavy"]==1 & chain_count_mtx[, "light"]>=1
        } else if (logic_num_HL=="1H") {
            # excatly 1 heavy, regardless of the number of light chain(s)
            # (could be 0 light chain)
            bool_num_HL = chain_count_mtx[, "heavy"]==1
        } else if (logic_num_HL=="1H_min1L_or_min1H_1L") {
            # exactly 1 heavy, at least 1 light
            # OR
            # at least 1 heavy, exactly 1 light
            bool_1H_min1L = chain_count_mtx[, "heavy"]==1 & chain_count_mtx[, "light"]>=1
            bool_min1H_1L = chain_count_mtx[, "heavy"]>=1 & chain_count_mtx[, "light"]==1
            bool_num_HL = bool_1H_min1L | bool_min1H_1L
        } else {
            warning("Unrecognized option for `logic`. `check_num_HL` skipped.\n")
            bool_num_HL = rep(T, nrow(chain_count_mtx))
        }
        
        cat("\ncheck_num_HL, number of cells:\n")
        cat("col:", col_cell, "\n")
        print(table(bool_num_HL), useNA="ifany")
        
        # map back to db
        uniq_cells_pass = uniq_cells[bool_num_HL]
        
        bool_num_HL_db = db[[col_cell]] %in% uniq_cells_pass
        
        # filter down to exactly 1 H and exactly 1 L
        if (logic_num_HL %in% c("1H_min1L", "1H_min1L_or_min1H_1L")) {
            # if there are cells with more than 1 H or L
            if ( sum(bool_num_HL_db) != length(uniq_cells_pass)*2 ) {
                
                # initialize
                bool_keep_seq = rep(T, nrow(db))
                
                # cells with >1 heavy, or {either >1 heavy or >1 light}
                
                chain_count_mtx_2 = chain_count_mtx[bool_num_HL, ]
                # sanity check: every cell in this table should have at least 1H
                # and/or at least 1L
                stopifnot(all( rowSums(chain_count_mtx_2>=1)==2 ))
                
                if (logic_num_HL=="1H_min1L") {
                    cells_min1 = rownames(chain_count_mtx_2)[chain_count_mtx_2[, "light"]>1]
                } else {
                    # in each cell, how many chains have count >1?
                    chain_count_mtx_2_rowsum = rowSums(chain_count_mtx_2>1)
                    # expect at most 1 chain has count >1
                    stopifnot(all(chain_count_mtx_2_rowsum<=1))
                    # cells with chain with count >1
                    cells_min1 = rownames(chain_count_mtx_2)[chain_count_mtx_2_rowsum==1]
                }
                
                # For each cell, set the boolean in bool_kep_seq for the non-majority
                # heavy/light (depending on which chain has >1) to F
                # If tied, keep the seq that appears earlier in the db (which.max) 
                for (cur_cell in cells_min1) {
                    # wrt db
                    idx_cell = which(db[[col_cell]]==cur_cell)
                    
                    # which has >1? heavy or light
                    cur_cell_locus_tab = table(db[[col_locus]][idx_cell])
                    # wrt cur_cell_locus_tab
                    i_tab_h = which(names(cur_cell_locus_tab)=="IGH")
                    
                    if (cur_cell_locus_tab[i_tab_h]>1) {
                        # if >1 heavy, must be only 1 light
                        stopifnot(sum(cur_cell_locus_tab[-i_tab_h])==1)
                        # keep most abundant heavy; disregard remaining heavy
                        # wrt db
                        idx_db_h = which(db[[col_cell]]==cur_cell & db[[col_locus]]=="IGH")
                        # wrt idx_db_h
                        idx_db_h_max = which.max(db[[col_umi]][idx_db_h])
                        bool_keep_seq[idx_db_h[-idx_db_h_max]] = F
                    } else {
                        # if 1 heavy, must be >1 light
                        stopifnot(sum(cur_cell_locus_tab[-i_tab_h])>1)
                        # keep most abundant light; disregard remaining light
                        # wrt db
                        idx_db_l = which(db[[col_cell]]==cur_cell & db[[col_locus]]!="IGH")
                        # wrt idx_db_l
                        idx_db_l_max = which.max(db[[col_umi]][idx_db_l])
                        bool_keep_seq[idx_db_l[-idx_db_l_max]] = F
                    }
                }
                # sanity check
                # there should be F's in bool_keep_seq (can't be still all T)
                stopifnot(!all(bool_keep_seq))
                
                bool_num_HL_db = bool_num_HL_db & bool_keep_seq
                
                # more sanity checks
                # dimension should still match nrow(db)
                stopifnot(length(bool_num_HL_db) == nrow(db))
                # number of cells passed should remain unchanged
                stopifnot( length(uniq_cells_pass) == 
                           length(unique(db[[col_cell]][bool_num_HL_db])) )
                # after filtering there should be exactly 1 H and exactly 1 L per cell
                stopifnot( sum(bool_num_HL_db) == length(uniq_cells_pass)*2 )
                stopifnot( sum(db[[col_locus]][bool_num_HL_db]=="IGH") == 
                               sum(db[[col_locus]][bool_num_HL_db]!="IGH") )
                
            } 
        }

    } else {
        bool_num_HL_db = rep(T, nrow(db))
    }
    
    bool = bool_num_HL_db
    
    # count
    cat("\nNumber of unique cells after cell-level QC:", 
        length(unique(db[[col_cell]][bool])), "\n")
    
    cat("\nAfter cell-level QC, number of seqs:\n")
    print(table(bool, useNA="ifany"))
    cat("\n")
    
    return(bool)
}


#' BCR sequence-level and/or cell-level QC
#' 
#' @param   db_name     Name of tab-separated file with headers that contains
#'                      input data
#' @param   seq_level   Boolean. Whether to perform sequence-level QC.
#' @param   cell_level  Boolean. Whether to perform cell-level QC.
#' @param   sequential  Boolean. Whether to perform sequence-level QC first 
#'                      before performing cell-level QC, as opposed to performing
#'                      both in parallel and then performing an `&` operation. 
#'                      Only applicable if both `seq_level` and `cell_level` are 
#'                      `TRUE`.
#' @param   outname     Stem of output filename. Prefix to 
#'                      `_qc.tsv`.
#' @param   outdir      Path to output directory.
#' @param   ...         All other parameters are passed to helper functions.        
#' 
#' @returns Writes a `[outname]_qc-pass/fail.tsv` to `outdir`.
#' 
#' @details When both `seq_level` and `cell_level` are `TRUE`, `sequential` being
#'          `FALSE` could give slightly different results than `TRUE`. For example,
#'          a cell may be linked to 2 light chains pre-QC. Suppose one of the two light
#'          chains gets filtered by sequence-level QC. Suppose that `logic_num_HL`
#'          is set to `1H_1L`. If cell-level QC is applied sequentially after
#'          one of the two light chains is filtered, this cell is considered by
#'          cell-level QC to be linked to 1 light chain only, and would pass QC.
#'          In contrast, if cell-level QC is applied in parallel with sequence-level
#'          QC, then this cell is considered by cell-level QC to be linked to 2
#'          light chains and would therefore fail cell-level QC. After the `&`
#'          operation, all sequences from this cell would then fail QC. 

perform_qc = function(db_name, seq_level=T, cell_level=F, sequential=F,
                      outname, outdir,
                      chain_type,
                      col_v_call, col_d_call, col_j_call, col_c_call,
                      check_valid_vj=F, 
                      check_chain_consistency=F, 
                      check_N=F, max_N, col_N, last_pos_N, as_perc_N,
                      check_nonATGC=F, col_obsv, col_germ,
                      max_nonATGC, last_pos_nonATGC, as_perc_nonATGC,
                      check_none_empty=F, col_none_empty,
                      check_NA=F, col_NA,
                      check_len_mod3=F, col_len_mod3,
                      col_locus, col_cell, col_umi,
                      check_locus,
                      check_num_HL, logic_num_HL) {
    
    db = read.table(db_name, header=T, sep="\t", stringsAsFactors=F)
    
    if (seq_level) {
        bool_seq = perform_qc_seq(db, chain_type,
                                  col_v_call, col_d_call, col_j_call, col_c_call,
                                  check_valid_vj, 
                                  check_chain_consistency, 
                                  check_N, max_N, col_N, last_pos_N, as_perc_N,
                                  check_nonATGC, col_obsv, col_germ,
                                  max_nonATGC, last_pos_nonATGC, as_perc_nonATGC,
                                  check_none_empty, col_none_empty,
                                  check_NA, col_NA,
                                  check_len_mod3, col_len_mod3)
    } else {
        bool_seq = rep(T, nrow(db))
    }
    
    if (cell_level) {
        
        if (seq_level & sequential) {
            cat("\nPerforming seq-level and cell-level QCs sequentially\n")
            
            # wrt db
            idx_use_for_cell_qc = which(bool_seq)
            
            if (length(idx_use_for_cell_qc)==0) {
                stop("No data left after sequence-level QC. Halted before cell-level QC.")
            }
            
        } else {
            # use all
            idx_use_for_cell_qc = 1:nrow(db)
        }
        
        # - If seq_level & sequential, bool_cell_intermediate's length could be shorter than 
        #   bool_seq, because db[idx_use_for_cell_qc, ] nrow could have changed due to
        #   subsetting by idx_use_for_cell_qc
    
        # - Otherwise, bool_cell_intermediate length should be the same as bool_seq,
        #   because db[idx_use_for_cell_qc, ] nrow would have stayed unchanged due to
        #   idx_use_for_cell_qc being 1 through nrow(db)
        
        # wrt db[idx_use_for_cell_qc, ]
        bool_cell_intermediate = perform_qc_cell(db[idx_use_for_cell_qc, ], chain_type,
                                                 col_locus, col_cell, col_umi,
                                                 check_locus, col_v_call,
                                                 check_num_HL, logic_num_HL)

        # map back to full db
        # any seq not included in idx_use_for_cell_qc will be FALSE
        # this ensures db, bool_seq, and bool_cell all have the same length/nrow

        bool_cell = rep(F, nrow(db))
        bool_cell[idx_use_for_cell_qc] = bool_cell_intermediate

    } else {
         bool_cell = rep(T, nrow(db))
    }
    
    stopifnot(length(bool_seq)==length(bool_cell))

    # cases
    # seq_level  cell_level
    # F          F
    # T          F
    # F          T
    # T          T           sequential
    # T          T           parallel
    if (seq_level & cell_level & sequential) {
        # T T sequential
        bool_final = bool_cell
    } else {
        # all other cases
        bool_final = bool_seq & bool_cell
    }

    stopifnot(length(bool_final)==nrow(db))
    
    setwd(outdir)
    
    if (any(bool_final)) {
        f = paste0(outname, "_qc-pass.tsv")
        write.table(db[bool_final, ], file=f, quote=F, sep="\t", 
                    row.names=F, col.names=T)
    }
    
    if (any(!bool_final)) {
        f = paste0(outname, "_qc-fail.tsv")
        write.table(db[!bool_final, ], file=f, quote=F, sep="\t", 
                    row.names=F, col.names=T)
    }
    
}


#' Post-QC split
#' 
#' Split by heavy/light and by productive/non-productive
#' 
#' @param   db_name     Name of tab-separated file with headers that contains
#'                      input data
#' @param   col_v_call  Column name for V gene annotation. 
#' @param   col_prod    Column name for productive/non-productive.
#' @param   val_prod    Value in `col_prod` indicating productive.
#' @param   outname     Stem of output filename. Prefix to 
#'                      `_[heavy|light]_[pr|npr].tsv`.
#' @param   outdir      Path to output directory.
#' 
#' @returns Writes `[outname]_[heavy|light]_[pr|npr].tsv` to `outdir`.
#' 
#' @details Currently only supports chain_type="IG" ("TR" not yet supported).
#'          Relies on "IG[HKL]" in V gene annotation to identify the locus.

split_db = function(db_name, col_v_call, col_prod, val_prod, outname, outdir) {
    
    db = read.table(db_name, header=T, sep="\t", stringsAsFactors=F)
    
    #bool_heavy = tolower(substr(db[[col_v_call]], 1, 3))=="igh"
    # to work with IMGT High/V-QUEST output which adds species name in front of V gene annotation
    bool_heavy = grepl(pattern="igh", x=db[[col_v_call]])

    bool_pr = db[[col_prod]]==val_prod
    
    db_heavy_pr = db[bool_heavy & bool_pr, ]
    db_heavy_npr = db[bool_heavy & !bool_pr, ]
    db_light_pr = db[!bool_heavy & bool_pr, ]
    db_light_npr = db[!bool_heavy & !bool_pr, ]
    
    setwd(outdir)
    
    if (nrow(db_heavy_pr)>0) {
        f = paste0(outname, "_heavy_pr.tsv")
        write.table(db_heavy_pr, file=f, quote=F, sep="\t", row.names=F, col.names=T)
        cat("\n# seqs for heavy & productive:", nrow(db_heavy_pr), "\n")
    } else {
        cat("\nNo data for heavy & productive.\n")
    }
    
    if (nrow(db_heavy_npr)>0) {
        f = paste0(outname, "_heavy_npr.tsv")
        write.table(db_heavy_npr, file=f, quote=F, sep="\t", row.names=F, col.names=T)
        cat("\n# seqs for heavy & non-productive:", nrow(db_heavy_npr), "\n")
    } else {
        cat("\nNo data for heavy & non-productive.\n")
    }
    
    if (nrow(db_light_pr)>0) {
        f = paste0(outname, "_light_pr.tsv")
        write.table(db_light_pr, file=f, quote=F, sep="\t", row.names=F, col.names=T)
        cat("\n# seqs for light & productive:", nrow(db_light_pr), "\n")
    } else {
        cat("\nNo data for light & productive.\n")
    }
    
    if (nrow(db_light_npr)>0) {
        f = paste0(outname, "_light_npr.tsv")
        write.table(db_light_npr, file=f, quote=F, sep="\t", row.names=F, col.names=T)
        cat("\n# seqs for light & non-productive:", nrow(db_light_npr), "\n")
    } else {
        cat("\nNo data for light & non-productive.\n")
    }
    
}
