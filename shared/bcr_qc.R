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
    
    # unpack multiple annotations (if any)
    # first split by ","
    # then extract first 3 chars
    chains_unpacked = sapply(chains, function(s){
        return( substr( strsplit(s, split=separator)[[1]], 1, 3) )
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
#' @param   check_perc_N                Boolean. Whether to check the % of 
#'                                      N's in `col_perc_N`.
#' @param   max_perc_N                  Max % of N's in `col_perc_N` allowed from
#'                                      position 1 thru `last_pos_N`.
#' @param   col_perc_N                  Column name(s) in which to perform 
#'                                      `check_perc_N`.
#' @param   last_pos_N                  Last position(s) in `col_perc_N` thru which
#'                                      to perform `check_perc_N`. Length should 
#'                                      match that of `col_perc_N`.
#' @param   check_num_nonATGCN          Boolean. Whether to check the number of 
#'                                      non-ATGCN positions in `col_obsv` using
#'                                      `col_germ` as reference.
#' @param   max_num_nonATGCN            Max number of non-ATGCN positions allowed
#'                                      in `col_obsv` from position 1 thru 
#'                                      `last_pos_nonATGCN`.
#' @param   last_pos_nonATGCN           The last position thru `col_obsv` to perform
#'                                      `check_num_nonATGCN`.
#' @param   check_none_empty            Boolean. Whether to check for `[Nn]one`
#'                                      and empty (`""`) values.
#' @param   col_none_empty              Column name(s) in which to perform 
#'                                      `check_none_empty`.                                      
#' @param   check_NA                    Boolean. Whether to check for `NA`.
#' @param   col_NA                      Column name(s) in which to perform `check_NA`.
#' @param   check_len_mod3              Boolean. Whether to check that lengths
#'                                      are a multiple of 3.
#' @param   col_len_mod3                Column name(s) in which to perform 
#'                                      `check_len_mod3`. Note that column(s) containing
#'                                      strings is/are expected. Lengths will be calculated
#'                                      on the fly.
#' 
#' @returns A new `db` subset to rows passing all checks of choice.
#'
#' @details - `check_valid_vj`: both V and J start with "IG" or "TR"
#'
#'          - `check_chain_consistency`: all of V, D, J, and C calls, 
#'            where supplied, are "IGH", "IGK", "IGL", or "TRA", etc.
#'            Ok to one or more, but not all, of `col_[vdjc]_call` as `NA`.
#'            Also see `inspect_chain_consistency`.
#'
#'          - `check_perc_N` checks if the % of positions that are N in 
#'            `col_perc_N` from position 1 thru `last_pos_N` is `<` (strictly)
#'            `max_perc_N`. This check is skipped for a row which has
#'            `""`, `"[Nn]one"`, or `NA` for `col_germ`.
#'
#'          - `check_num_nonATGCN` checks if the number of positions that are
#'            non-ATGCN in `col_obsv`, between position 1 thru `last_pos_nonATGCN`,
#'            at positions that are ATGCN in `col_germ`, is `<` (strictly)
#'            `max_num_nonATGCN`. This check is skipped for a row which has
#'            `""`, `"[Nn]one"`, or `NA` for `col_germ`.
#'
#'          - `check_none_empty` checks if `col_none_empty` is `[Nn]one` or `""`.
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
                          check_perc_N=F, max_perc_N, col_perc_N, last_pos_N,
                          check_num_nonATGCN=F, col_obsv, col_germ,
                          max_num_nonATGCN, last_pos_nonATGCN,
                          check_none_empty=F, col_none_empty,
                          check_NA=F, col_NA,
                          check_len_mod3=F, col_len_mod3) {
    
    require(stringi)
    require(seqinr)
    
    if (check_perc_N) {
        stopifnot( length(col_perc_N) == length(last_pos_N) )
    }
    
    nseqs = nrow(db)
    cat("\nPerforming sequence-level QC on", nseqs, "sequences\n")
    
    if (check_valid_vj) {
        # IG[HKL][VDJ]... or TR[ABDG][VDJ]...
        pattern_valid = paste0("^", chain_type)
        
        # TRUE means valid 
        bool_valid_v = grepl(pattern=pattern_valid, x=db[[col_v_call]])
        bool_valid_j = grepl(pattern=pattern_valid, x=db[[col_j_call]])
        
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
            inspect_chain_consistency(vg=ifelse(bool_v, db[[col_v_call]][i], NULL),
                                      dg=ifelse(bool_d, db[[col_d_call]][i], NULL),
                                      jg=ifelse(bool_j, db[[col_j_call]][i], NULL),
                                      cg=ifelse(bool_c, db[[col_c_call]][i], NULL),
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
    
    
    if (check_perc_N) {
        
        # check that all columns to be checked are characters
        stopifnot(all( sapply(col_perc_N, 
                              function(s){ class(db[[s]]) }) == "character" ))
        
        # matrix, even when col_perc_N is of length 1
        # each col is a col in col_perc_N
        # each row is a seq in db
        # TRUE means % of N's in pos 1 thru last_pos_N < max_perc_N
        mtx_perc_N = do.call(cbind, 
                             sapply(1:length(col_perc_N), function(i){
                                 s = col_perc_N[i]
                                 p = last_pos_N[i]
                                 
                                 # skip if NA, "", "[Nn]one"
                                 bool_skip = is.na(db[[s]]) | db[[s]]=="" | tolower(db[[s]])=="none"
                                 
                                 # convert to uppercase so no need to deal with cases
                                 # truncate to pos 1 to last_pos_N
                                 cur_s = toupper(substr(db[[s]], 1, p))
                                 # count number of occurrences of N characters
                                 cur_s_count = stri_count_fixed(str=cur_s,
                                                                pattern="N")
                                 # calc %
                                 cur_s_perc = cur_s_count / p * 100
                                 
                                 return( (cur_s_perc < max_perc_N) | bool_skip )
                             }, simplify=F))
        # sum across columns
        # if all TRUE, rowSums should be the same as number of cols
        bool_perc_N = rowSums(mtx_perc_N, na.rm=T) == ncol(mtx_perc_N)
        
        # count
        cat("\ncheck_perc_N ( max_perc_N=", max_perc_N, "; strict <):\n")
        cat("cols:", col_perc_N, "\n")
        print(table(bool_perc_N, useNA="ifany"))
        cat("\n")
        
    } else {
        bool_perc_N = rep(T, nseqs)
    }
    
    
    if (check_num_nonATGCN) {
        
        bool_skip = is.na(db[[col_germ]]) | db[[col_germ]]=="" | tolower(db[[col_germ]])=="none"
        
        # convert to uppercase so no need to deal with cases
        # truncate to pos 1 to last_pos_nonATGCN
        germ_upper_trunc = sapply(db[[col_germ]], function(s){
            toupper( substr(s, 1, min(nchar(s), last_pos_nonATGCN) ) )
        })
        # same for observed
        obsv_upper_trunc = sapply(db[[col_obsv]], function(s){
            toupper( substr(s, 1, min(nchar(s), last_pos_nonATGCN) ) )
        })
        
        # positions in germline from pos 1 thru last_pos_nonATGCN 
        # that are non-ATGC (e.g. ".")
        # this should be a list
        germ_idx = sapply(germ_upper_trunc, function(s){
            return(which(!s2c(s) %in% c("A","T","G","C","N")))
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
            
            cur_obsv_count = stri_count_regex(str=cur_obsv, pattern="[^ATGCN]")
            
            return(cur_obsv_count)
        })
        
        # TRUE means # of non-ATGCs in pos 1 thru last_pos_nonATGCN < max_num_nonATGCN
        # skip check if germline missing (and set to TRUE for that row)
        bool_num_nonATGCN = (obsv_count < max_num_nonATGCN) | bool_skip
          
        # count
        cat("\ncheck_num_nonATGCN ( max_num_nonATGCN=", max_num_nonATGCN, "; strict <):\n")
        cat("observed col:", col_obsv, "; germline col:", col_germ, "\n")
        print(table(bool_num_nonATGCN, useNA="ifany"))
        cat("\n")
         
        
    } else {
        bool_num_nonATGCN = rep(T, nseqs)
    }
    
    
    if (check_none_empty) {

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
    stopifnot(!any(is.na(bool_perc_N)))
    stopifnot(!any(is.na(bool_num_nonATGCN)))
    stopifnot(!any(is.na(bool_none_empty)))
    stopifnot(!any(is.na(bool_NA)))
    stopifnot(!any(is.na(bool_len_mod3)))

    bool = bool_valid_vj & bool_chain & bool_perc_N & bool_num_nonATGCN & bool_none_empty & bool_NA & bool_len_mod3
    
    # count
    cat("\nCombined:\n")
    print(table(bool, useNA="ifany"))
    cat("\n")
    
    # subset
    return(db[bool, ])
}


# TODO
# perform_qc_cell = function() { }


#' BCR QC & post-QC split
#' 
#' @param   db          data.frame
#' @param   seq_level   Boolean. Whether to perform sequence-level QC.
#' @param   cell_level  Boolean. Whether to perform cell-level QC.
#' @param   split       Boolean. Whether to split by heavy/light and by 
#'                      productive/non-productive.
#' @param   col_prod    Column name for productive/non-productive.
#' @param   val_prod    Value in `col_prod` indicating productive.
#' @param   outname     Stem of output filename. Prefix to 
#'                      `_qc_[heavy|light]_[pr|npr].tsv`.
#' @param   outdir      Path to output directory.
#' @param   ...         All other parameters are passed to helper functions.        

perform_qc = function(db, seq_level=T, cell_level=F, 
                      split=T, col_prod, val_prod, outname, outdir,
                      col_v_call, col_d_call, col_j_call, col_c_call,
                      check_valid_vj=F, 
                      check_chain_consistency=F, 
                      check_perc_N=F, max_perc_N, col_perc_N, last_pos_N,
                      check_num_nonATGCN=F, col_obsv, col_germ,
                      max_num_nonATGCN, last_pos_nonATGCN,
                      check_none_empty=F, col_none_empty,
                      check_NA=F, col_NA,
                      check_len_mod3=F, col_len_mod3) {
    
    # currently only supporing IG, not TR
    chain_type = "IG"
    
    if (seq_level) {
        db = perform_qc_seq(db, 
                            chain_type=chain_type,
                            col_v_call, col_d_call, col_j_call, col_c_call,
                            check_valid_vj, 
                            check_chain_consistency, 
                            check_perc_N, max_perc_N, col_perc_N, last_pos_N,
                            check_num_nonATGCN, col_obsv, col_germ,
                            max_num_nonATGCN, last_pos_nonATGCN,
                            check_none_empty, col_none_empty,
                            check_NA, col_NA,
                            check_len_mod3, col_len_mod3)
    }
    
    if (cell_level) {
        # TODO
        # db = perform_qc_cell(db)
    }
    
    
    setwd(outdir)
    
    if (split) {
        bool_heavy = tolower(substr(db[[col_v_call]], 1, 3))=="igh"
        bool_pr = db[[col_prod]]==val_prod
        
        db_heavy_pr = db[bool_heavy & bool_pr, ]
        db_heavy_npr = db[bool_heavy & !bool_pr, ]
        db_light_pr = db[!bool_heavy & bool_pr, ]
        db_light_npr = db[!bool_heavy & !bool_pr, ]
        
        if (nrow(db_heavy_pr)>0) {
            f = paste0(outname, "_qc_heavy_pr.tsv")
            write.table(db_heavy_pr, file=f, quote=F, sep="\t", row.names=F, col.names=T)
        } else {
            cat("\nNo data left for heavy & productive.\n")
        }
        
        if (nrow(db_heavy_npr)>0) {
            f = paste0(outname, "_qc_heavy_npr.tsv")
            write.table(db_heavy_npr, file=f, quote=F, sep="\t", row.names=F, col.names=T)
        } else {
            cat("\nNo data left for heavy & non-productive.\n")
        }
        
        if (nrow(db_light_pr)>0) {
            f = paste0(outname, "_qc_light_pr.tsv")
            write.table(db_light_pr, file=f, quote=F, sep="\t", row.names=F, col.names=T)
        } else {
            cat("\nNo data left for light & productive.\n")
        }
        
        if (nrow(db_light_npr)>0) {
            f = paste0(outname, "_qc_light_npr.tsv")
            write.table(db_light_npr, file=f, quote=F, sep="\t", row.names=F, col.names=T)
        } else {
            cat("\nNo data left for light & non-productive.\n")
        }
        
    } else {
        f = paste0(outname, "_qc.tsv")
        write.table(db, file=f, quote=F, sep="\t", row.names=F, col.names=T)
    }
    
}

