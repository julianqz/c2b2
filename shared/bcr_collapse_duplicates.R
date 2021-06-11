# Functions to perform the task of collapsing duplicate IMGT-aligned sequences

# notes to self on how collapseDuplicates behave
RUN=F
if (RUN) {
    
    # 1) behavior around "" in text_fields
    # - before collapsing, "" in c_call[1:3] dropped 
    # - after collapsing,  "" not in c_call for the collapse of [1:3]
    
    # 2) behavior around NA in text_fields
    # - before collapsing, NA in c_call[5:6] dropped 
    # - after collapsing,  NA not in c_call for the collapse of [5:6]
    
    # 3) behavior around sequences with varying lengths
    # - before collapsing, [8] and [9] differ in lengths; [8] is a substring of [9]
    # - [8] and [9] are not collapsed; a warning about length difference is displayed
    
    # [1:3] collapsed
    # [4:6] collapsed
    # [7] not collapsed 
    # [8:9] not collapsed due to difference in length
    
    tmp = data.frame(sequence_id=LETTERS[1:9],
                     sequence_alignment=c("CCCCTGGG","CCCCTGGN","CCCCTGNN",
                                          "ATGCATGC","ATGCATGC","ATGCATGC",
                                          "NAACTGGN",
                                          "TTTTAAAG","TTTTAAAGC"),
                     c_call=c("IGHM", "IGHG", "",
                              "IGHD", NA, NA,
                              "IGHA", 
                              "IGHE", "IGHE"),
                     sample_id=c("S1", "S1", "S2", 
                                 "S4", "S4", "S4",
                                 "S2",
                                 "S3", "S3"),
                     duplicate_count=1:9,
                     stringsAsFactors=FALSE)
    
    alakazam::collapseDuplicates(tmp, 
                                 text_fields=c("c_call", "sample_id"), 
                                 num_fields="duplicate_count", 
                                 add_count=TRUE, verbose=TRUE)
}


#' Call alakazam::collapseDuplicate on a by-clone basis with enhanced options
#' 
#' @param   db               A `data.frame`.
#' @param   col_clone        Column containing clone ID.
#' @param   col_seq          Column in which collapse is to be performed. 
#'                           Passed to `seq` in `collapseDuplicates`.
#' @param   col_id           Column containing sequence ID. 
#'                           Passed to `id` in `collapseDuplicates`.
#' @param   col_text_fields  Passed to `text_fields` in `collapseDuplicates`. Optional.
#' @param   col_num_fields   Passed to `num_fields` in `collapseDuplicates`. Optional.
#' @param   col_seq_fields   Passed to `seq_fields` in `collapseDuplicates`. Optional.
#' @param   col_preserve     Column containing info on which sequences to preserve
#'                           without collapsing. Optional. If supplied, must be
#'                           a single value.
#' @param   val_preserve_vec Value(s) in `col_preserve` that will trigger preservation.
#'                           Optional. If supplied, can be one or more values.
#' @param   col_distinct_vec Column(s) which will be used to keep sequence-level 
#'                           duplicates separate/uncollapsed. Optional. 
#'                           If supplied, can be one or more values.
#' 
#' @return  An updated db with collapse performed and a new column, `collapse_count`.
#' 
#' @details Be sure to check out the helpful doc of `collapseDuplicates`!
#' 
run_collapse_duplicates = function(db, nproc=1,
                                   col_clone, col_seq, col_id,
                                   col_text_fields=NULL, 
                                   col_num_fields=NULL, 
                                   col_seq_fields=NULL,
                                   col_preserve=NULL, val_preserve_vec=NULL,
                                   col_distinct_vec=NULL) {
    
    suppressPackageStartupMessages(require(stringi))
    suppressPackageStartupMessages(require(alakazam))
    suppressPackageStartupMessages(require(doParallel))
    suppressPackageStartupMessages(require(foreach))
    
    #### make_unique ####
    # located within run_collapse_duplicates so it can be exported to cluster
    
    #' Make entries unique
    #' 
    #' @param  vec  A vector. May or may not contain one or more distinct duplicates.
    #' 
    #' @return  A vector. If applicable, the return vector will be updated 
    #'          with all entries now being unique.
    #'         
    #' @details Suffix _1/2/3/... are added to duplicate entries. 
    #' 
    make_unique = function(vec) {
        # make a table
        tab = table(vec)
        
        # any dup?
        # - no duplicate, return unchanged
        # - otherwise, make unique and return updated vec
        if (any(tab>1)) {
            
            # all distinct duplicates
            dups = names(tab[tab>1])
            
            # for each distinct duplicate
            for (i in 1:length(dups)) {
                cur_dup = dups[i]
                cur_idx = which(vec==cur_dup)
                # add suffix to all entries
                vec[cur_idx] = paste0(cur_dup, "_", 1:length(cur_idx))
            }
            
            # sanity check: there should no more duplicate
            stopifnot( !base::anyDuplicated(vec) )
        } 
        return(vec)
    }
    
    
    #### checks ####
    
    stopifnot(nproc>=1)
    
    if (!is.null(col_preserve) | !is.null(col_distinct_vec)) {
        stopifnot( all( c(col_preserve, col_distinct_vec) %in% colnames(db) ) )
    }
    
    if (!is.null(col_preserve)) {
        # expect exactly 1 column
        stopifnot( length(col_preserve)==1 )
        # expect at least 1 non-NA value in val_preserve_vec
        stopifnot( !is.null(val_preserve_vec) )
        stopifnot( sum(!is.na(val_preserve_vec))>=1 )
    }
    
    n_clones_bf = length(unique(db[[col_clone]]))
    cat("\nBefore collapsing - # clones:",
        n_clones_bf, "\n")
    cat("\nBefore collapsing - # seqs:",
        nrow(db), "\n")
    
    #### preserve part of db ####
    # no collapsing
    
    if (!is.null(col_preserve)) {
        
        cat("\nBefore collapsing -", col_preserve, ":\n")
        print(table(db[[col_preserve]]))
        
        lst_bool_val = vector(mode="list", length=length(val_preserve_vec))
        lst_db_val = vector(mode="list", length=length(val_preserve_vec))
        names(lst_bool_val) = as.character(val_preserve_vec)
        names(lst_db_val) = as.character(val_preserve_vec)
        
        for (val in val_preserve_vec) {
            
            bool_val = db[[col_preserve]]==val
            lst_bool_val[[as.character(val)]] = bool_val
            
            if (any(bool_val)) {
                db_val = db[bool_val, ]
                
                cat("\nBefore collapsing - # seqs with",
                    col_preserve, "=", val, ":", 
                    nrow(db_val), "(preserved)\n")
            } else {
                db_val = NULL
            }
            
            lst_db_val[[as.character(val)]] = db_val
            rm(db_val, bool_val)
        }
        
        mtx_bool_val = do.call(cbind, lst_bool_val)
        # OR across bool_val
        bool_preserve = rowSums(mtx_bool_val)>=1
        
        db_preserve = do.call(rbind, lst_db_val)
        
        # ok if all preserved with no more seq left to collapse
        # if so, nrow(db_to_collapse) will be 0 and handled subsequently
        db_to_collapse = db[!bool_preserve, ]
        if (all(bool_preserve)) {
            stopifnot( nrow(db_to_collapse)==0 )
        }

        rm(mtx_bool_val, lst_bool_val, lst_db_val, bool_preserve)
        
    } else {
        cat("\nNothing to preserve.\n")
        db_to_collapse = db
        db_preserve = NULL
    }
    rm(db)
    
    #### collapse ####
    
    if (nrow(db_to_collapse)>0) {
        
        # collapse
        cat("\nCollapsing...\n")
        
        if (!is.null(col_distinct_vec)) {
            
            cat("\nKeeping sequence-level duplicates separate by:", 
                col_distinct_vec, "\n")
            
            # A hack to prevent duplicates derived from different `col_distinct_vec`
            # (e.g. subset, timepoint, etc.) from being collapsed
            
            # Create an artificial column by combining sequence and 
            # values of `col_distinct_vec` columns
            
            # Because collapseDuplicates requires that input seqs have the same length,
            # use one-letter code to represent values in `col_distinct_vec` columns
            # (since using full names could introduce varying lengths)
            
            lst_encoded = vector(mode="list", length=length(col_distinct_vec))
            names(lst_encoded) = col_distinct_vec
            
            for (col_d in col_distinct_vec) {
                # unique values
                col_d_uniq_vals = sort(unique(db_to_collapse[[col_d]]))
                # only support <=25 unique values
                # N is excluded from LETTERS because it is part of `ignore`
                CODES = LETTERS[-which(LETTERS=="N")]
                if (length(col_d_uniq_vals)>length(CODES)) {
                    stop("Currently only supprting <=", length(CODES), 
                         " unique values in `", col_d, "`\n")
                }
                # one-letter code
                col_d_code = CODES[1:length(col_d_uniq_vals)]
                names(col_d_code) = as.character(col_d_uniq_vals)
                # encode
                lst_encoded[[col_d]] = col_d_code[db_to_collapse[[col_d]]]
            }
            
            # each row is a sequence
            # cols: sequence, col_distinct_vec in encoded form
            db_to_collapse_encoded = data.frame(cbind(db_to_collapse[[col_seq]], 
                                                      do.call(cbind, lst_encoded)),
                                                stringsAsFactors=F, row.names=NULL)
            # concat
            col_tmp = "CD_TMP"
            col_tmp_separator = "@"
            db_to_collapse[[col_tmp]] = apply(db_to_collapse_encoded, 1, 
                                              stri_join, collapse=col_tmp_separator)
            
            stopifnot(!any(is.na(db_to_collapse[[col_tmp]])))
            
            rm(db_to_collapse_encoded, lst_encoded)
            
        } else {
            db_to_collapse[[col_tmp]] = db_to_collapse[[col_seq]]
        }
        
        
        ### collapse on a by-clone basis
        
        clones = sort(unique(db_to_collapse[[col_clone]]))
        
        ## set up cluster
        
        # Create cluster of nproc size and export namespaces
        # If user wants to parallellize this function and specifies nproc > 1, then
        # initialize and register slave R processes/clusters & 
        # export all necesseary environment variables, functions and packages.
        if (nproc==1) {
            # If needed to run on a single core/cpu then, registerDoSEQ
            # Without doing this, foreach will give warning (though will still run)
            registerDoSEQ()
        } else if (nproc>1) {
            cluster = parallel::makeCluster(nproc, type="PSOCK")
            registerDoParallel(cluster)
            
            # export to cluster
            export_functions <- list("clones", 
                                     "db_to_collapse",
                                     "col_clone",
                                     "col_seq",
                                     "col_id",
                                     "col_tmp",
                                     "col_tmp_separator",
                                     "col_text_fields",
                                     "col_num_fields",
                                     "col_seq_fields",
                                     "make_unique",
                                     "collapseDuplicates"
            )
            parallel::clusterExport(cluster, export_functions, envir=environment())
        }
        
        
        ## loop thru clones
        lst_collapse = foreach(i=1:length(clones)) %dopar% {
            
            cur_db = db_to_collapse[db_to_collapse[[col_clone]]==clones[i], ]
            
            # collapseDuplicates() requires that all entries in the ID field be unique
            # this may not always be the case
            # e.g. a seq from sample 5 and a seq from sample 6 (both samples from the same subject)
            #      could happen to have the same ID (UMI=GGTTATGATAAAGGGTT)
            #      => append suffix so these become unique: GGTTATGATAAAGGGTT_[12]
            if (length(unique(cur_db[[col_id]])) < nrow(cur_db)) {
                cur_db[[col_id]] = make_unique(cur_db[[col_id]])
            }
            
            # text_fields: Character vector of textual columns to collapse
            
            # num_fields:  Vector of numeric columns to collapse
            #              The numeric annotations of duplicate sequences will be summed
            
            # seq_fields:  Vector of nucleotide sequence columns to collapse
            #              The sequence with the fewest number of non-informative characters will be retained
            #              Note, this is distinct from the seq parameter which is used to determine duplicates.
            
            cur_db_collapse = collapseDuplicates(data=cur_db,
                                                 id=col_id,
                                                 seq=col_tmp,
                                                 text_fields=col_text_fields,
                                                 num_fields=col_num_fields,
                                                 seq_fields=col_seq_fields,
                                                 add_count=T, # adds $collapse_count
                                                 ignore=c("N","-",".","?"), # default
                                                 sep=",", # default
                                                 verbose=F) # default
            
            # extract col_seq back from $TMP
            cur_db_collapse[[col_seq]] = sapply(cur_db_collapse[[col_tmp]], 
                                                function(s){
                                                    strsplit(s, col_tmp_separator)[[1]][1]
                                                }, USE.NAMES=F)
            
            # remove $TMP (no longer needed)
            cur_db_collapse[[col_tmp]] = NULL
            
            rm(cur_db)
            
            return(cur_db_collapse)
        }
        
        ## stop the cluster
        if (nproc>1) { parallel::stopCluster(cluster) }
        
        rm(db_to_collapse)
        
        db_collapse = do.call(rbind, lst_collapse)
        rm(lst_collapse)
        
    } else {
        # nothing to collapse
        cat("\nNo seq left to collapse. Skipped collapsing.\n")
        db_collapse = NULL
    }

    # can't both be NULL
    stopifnot( ! (is.null(db_collapse) & is.null(db_preserve)) )
    
    # add missing col in db_preserve (added by collapseDuplicates to db_collapse)
    if (!is.null(db_preserve)) {
        db_preserve[["collapse_count"]] = 1
    }
    
    if ( (!is.null(db_collapse)) & (!is.null(db_preserve)) ) {
        stopifnot(all.equal(colnames(db_preserve),
                            colnames(db_collapse)))
    }
    
    # ok if one of the two is NULL
    db_new = rbind(db_preserve, db_collapse)
    rm(db_preserve, db_collapse)
    
    n_clones_af = length(unique(db_new[[col_clone]]))
    cat("\nAfter collapsing - # clones",
        n_clones_af, "\n")
    # total number of clones should stay the same
    stopifnot( n_clones_bf == n_clones_af )
    
    cat("\nAfter collapsing - # seqs:",
        nrow(db_new), "\n")
    
    if (!is.null(col_preserve)) {
        cat("\nAfter collapsing -", col_preserve, ":\n")
        print(table(db_new[[col_preserve]]))
    }
    
    return(db_new)
}


#' Verify the collapse performed by `run_collapse_duplicates`.
#' 
#' @param  db                A collapsed `data.frame` outputted by 
#'                           `run_collapse_duplicates`.
#' @param  N                 A positive integer. The top N biggest clones to 
#'                           perform verification on.
#' @param  col_clone         Same as that for `run_collapse_duplicates`.
#' @param  col_seq           Same as that for `run_collapse_duplicates`.
#' @param  col_preserve      Same as that for `run_collapse_duplicates`.
#' @param  val_preserve_vec  Same as that for `run_collapse_duplicates`.
#' @param  col_distinct_vec  Same as that for `run_collapse_duplicates`.
#' 
#' @return Prints a console message indicating the verification outcome.
#' 
verify_collapse_duplicates = function(db, N,
                                      col_clone, col_seq,
                                      col_preserve=NULL, val_preserve_vec=NULL,
                                      col_distinct_vec=NULL) {
    
    suppressPackageStartupMessages(require(stringi))
    
    # helper func
    count_uniq = function(db, col_seq, col_distinct_vec=NULL) {
        
        if (!is.null(col_distinct_vec)) {
            cols = c(col_seq, col_distinct_vec)
            vec = apply(db[, cols], 1, stri_join, collapse="@")
        } else {
            vec = db[[col_seq]]
        }
        
        n_uniq = length(unique(vec))
        
        return(n_uniq)
    }
    
    cat("\nPerforming verification...\n")
    
    # remove rows to preserve
    if (!is.null(col_preserve)) {
        for (val in val_preserve_vec) {
            db = db[db[[col_preserve]]!=val, ]
        }
    }
    
    if (nrow(db)>0) {
        
        # check top N biggest clones
        tab = sort(table(db[[col_clone]]), decreasing=T)
        # `min` in case there're fewer than N clones
        clones_ck = names(tab)[1:min(N, length(tab))]
        
        bool = rep(F, length(clones_ck))
        
        for (i in 1:length(clones_ck)) {
            
            cur_db = db[db[[col_clone]]==clones_ck[i], ]
            
            n_uniq = count_uniq(db=cur_db, col_seq=col_seq,
                                col_distinct_vec=col_distinct_vec)
            
            bool[i] = nrow(cur_db)==n_uniq
            
            rm(cur_db, n_uniq)
        }
        
        if (all(bool)) {
            cat("\nAll verification PASSED.\n")
        } else {
            cat("\nVerification FAILED for:", clones_ck[!bool], "\n")
        }
        
    } else {
        cat("\nNo seq left to collapse after preserving. Verification skipped.\n")
    }
    
    rm(db)
}

