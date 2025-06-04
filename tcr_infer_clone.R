# Julian Q. Zhou
# https://github.com/julianqz


#' Infer TCR clonal relationships
#'
#' @param  db                A `data.frame`.
#' @param  v_call            Column name for V gene annotation.
#' @param  j_call            Column name for J gene annotation.
#' @param  junc_len          Column name for junction length.
#' @param  junc              Column name for junction. This can be either nt or aa.
#' @param  cell_id           Column name for cell ID.
#' @param  locus             Column name for locus. Expect values to be one of
#'                           `{TRB, TRA, TRD, TRG}`.
#' @param  single_cell_mode  Boolean. If `TRUE`, run in single-cell mode.
#' @param  use_only_vdj      Boolean. Only applies in single-cell mode. 
#'                           If `TRUE`, infer clonal relationships based on only
#'                           VDJ sequences, ignore VJ sequences. See details below.
#' @param  clone_id          Column name to be created for clone ID.
#'
#' @returns  A list containing `db_passed` and `db_failed`, both `data.frame`.
#'           `db_passed` contains an added column `clone_id`, and its number of rows 
#'           may be fewer than the input `db`. See details below.
#'           `db_failed` contains, if any, the rows from the input `db` that were removed.
#'
#' @details  In essence, TCR clonal relationships are inferred by first grouping 
#'           cells/sequences by VJL (V gene, J gene, junction length), and then, 
#'           within each VJL partition, by identical junction (nt or aa).
#'           
#'           If `single_cell_mode=FALSE`, bulk mode is run and all sequences in `db`
#'           are used for inferring clonal relationships, ignoring `cell_id` and `locus`.
#'           
#'           If `single_cell_mode=TRUE`, single-cell mode is run. 
#'           
#'           - First, cells which do not have both VDJ and VJ sequences are removed.
#'           This is why the output `db_passed` could potentially have fewer rows than
#'           the input `db`. It is assumed that each of the cells who do have both VDJ 
#'           and VJ sequences has 1 VDJ and 1 VJ sequence.
#'           
#'           - Next, clonal inference is performed. The inference can be based on both
#'           VDJ and VJ sequences, or VDJ sequences alone, depending on the value of 
#'           `use_only_vdj`. In the case that inference is performed based on VDJ 
#'           sequences only, the resultant clone IDs are still propagated to all 
#'           the corresponding VJ sequences. 
#'           
define_tcr_clone = function(db, v_call, j_call, junc_len, junc, cell_id, locus, 
                            single_cell_mode, use_only_vdj, clone_id) {
    
    require(alakazam)
    
    vec_loci_vdj = c("TRB", "TRD")
    vec_loci_vj = c("TRA", "TRG")
    col_vjl_group_alakazam = "vj_group" # from alakazam::groupGenes
    col_vjl_group = "vjl_group" # rename
    
    cat("\nColumn for junction =", junc, "\n")
    cat("\nSingle-cell mode =", single_cell_mode, "\n")
    
    cat("\nInput db has", nrow(db), "seqs\n")
    
    ### filtering of db, if necessary
    
    if (single_cell_mode) {
        
        stopifnot(!is.null(cell_id))
        
        # check locus values
        stopifnot(all(db[[locus]] %in% c(vec_loci_vdj, vec_loci_vj)))
        
        vec_uniq_cell_id = unique(db[[cell_id]])
        
        # cells with VDJ seqs
        vec_cell_id_vdj = db[[cell_id]][db[[locus]] %in% vec_loci_vdj]
        
        # cells with VJ seqs
        vec_cell_id_vj = db[[cell_id]][db[[locus]] %in% vec_loci_vj]
        
        bool_uniq_cell_id = ((vec_uniq_cell_id %in% vec_cell_id_vdj) & 
                             (vec_uniq_cell_id %in% vec_cell_id_vj))
        
        if (any(!bool_uniq_cell_id)) {
            
            vec_uniq_cell_id_rmv = vec_uniq_cell_id[!bool_uniq_cell_id]
            idx_db_rmv = which(db[[cell_id]] %in% vec_uniq_cell_id_rmv)
            
            cat("\nRemoved due to not having both VDJ and VJ seqs per cell:",
                length(vec_uniq_cell_id_rmv), "cells (",
                length(idx_db_rmv), "seqs)\n")
            
            db_failed = db[idx_db_rmv, ]
            
            db = db[-idx_db_rmv, ]
            
        } else {
            db_failed = NULL
        }
    } else {
        db_failed = NULL
    }
    
    
    ### VJL grouping
    
    # ?groupGenes
    # To invoke single-cell mode the cell_id argument must be specified and 
    # the locus column must be correct. 
    # Otherwise, groupGenes will be run with bulk sequencing assumptions, using all input sequences 
    # regardless of the values in the locus column.
    
    cat("\nPerforming VJL grouping...\n")
    
    if (single_cell_mode) {
        
        cat("\nuse_only_vdj =", use_only_vdj, "\n")
        
        # use_only_vdj=FALSE: infer relationships based on both VDJ and VJ seqs
        # use_only_vdj=TRUE:  infer relationships based on only VDJ seqs (despite being data being single-cell)
        # Either way, $vj_group values are propagated to both VDJ and VJ seqs
        
        db_vjl = groupGenes(data=db,
                            v_call=v_call,
                            j_call=j_call,
                            junc_len=junc_len,
                            cell_id=cell_id,
                            locus=locus,
                            only_heavy=use_only_vdj,
                            first=F)
        
    } else {
        # bulk mode
        # once cell_id=NULL, locus and only_heavy both ignored
        
        db_vjl = groupGenes(data=db,
                            v_call=v_call,
                            j_call=j_call,
                            junc_len=junc_len,
                            cell_id=NULL,
                            first=F)
    }
    
    # groupGenes adds a $vj_group column (note: not $vjl_group)
    # rename it to $vjl_group
    idx_col_vjl = which(colnames(db_vjl)==col_vjl_group_alakazam)
    stopifnot(length(idx_col_vjl)==1)
    colnames(db_vjl)[idx_col_vjl] = col_vjl_group
    
    ### prep format for helper func, if necessary
    
    vec_col_db_formatted_basic = c(cell_id, col_vjl_group, junc)
    
    if (single_cell_mode) {
        
        if (use_only_vdj) {
            # sc mode but only using vdj
            # keep only vdj seqs
            db_vjl_formatted = db_vjl[db_vjl[[locus]] %in% vec_loci_vdj, vec_col_db_formatted_basic]
                
        } else {
            # convert to wide format
            
            # Note that this could be different from vec_uniq_cell_id from above,
            # due to possible removal of cells from db
            vec_uniq_cell_id_2 = unique(db_vjl[[cell_id]])
            
            col_junc_vdj = paste0(junc, "_vdj")
            col_junc_vj = paste0(junc, "_vj")
            
            # vec_col_db_formatted_basic has junc instead of col_junc_vdj and col_junc_vj
            vec_col_db_formatted = c(cell_id, col_vjl_group, col_junc_vdj, col_junc_vj)
            
            db_vjl_formatted = data.frame(matrix(NA, nrow=length(vec_uniq_cell_id_2), 
                                                 ncol=length(vec_col_db_formatted)))
            colnames(db_vjl_formatted) = vec_col_db_formatted
            
            db_vjl_formatted[[cell_id]] = vec_uniq_cell_id_2
            
            # temp df's
            db_vjl_tmp_vdj = db_vjl[db_vjl[[locus]] %in% vec_loci_vdj, vec_col_db_formatted_basic]
            db_vjl_tmp_vj = db_vjl[db_vjl[[locus]] %in% vec_loci_vj, vec_col_db_formatted_basic]
            
            # wrt db_vjl_tmp_vdj
            idx_cell_id_vdj = match(vec_uniq_cell_id_2, db_vjl_tmp_vdj[[cell_id]])
            stopifnot(!any(is.na(idx_cell_id_vdj)))
            stopifnot(all.equal( vec_uniq_cell_id_2, db_vjl_tmp_vdj[[cell_id]][idx_cell_id_vdj] ))
            
            db_vjl_formatted[[col_vjl_group]] = db_vjl_tmp_vdj[[col_vjl_group]][idx_cell_id_vdj]
            db_vjl_formatted[[col_junc_vdj]] = db_vjl_tmp_vdj[[junc]][idx_cell_id_vdj]
            
            # wrt db_vjl_tmp_vj
            idx_cell_id_vj = match(vec_uniq_cell_id_2, db_vjl_tmp_vj[[cell_id]])
            stopifnot(!any(is.na(idx_cell_id_vj)))
            stopifnot(all.equal( vec_uniq_cell_id_2, db_vjl_tmp_vj[[cell_id]][idx_cell_id_vj] ))
            
            stopifnot(all.equal( db_vjl_formatted[[col_vjl_group]],
                                 db_vjl_tmp_vj[[col_vjl_group]][idx_cell_id_vj] ))
            
            db_vjl_formatted[[col_junc_vj]] = db_vjl_tmp_vj[[junc]][idx_cell_id_vj]
            
            stopifnot(!any(is.na(db_vjl_formatted)))
        }
        
    } else {
        # bulk mode; no formatting needed
        db_vjl_formatted = db_vjl[, vec_col_db_formatted_basic]
    }
    
    
    ### infer clonal relationships
    
    if (single_cell_mode) {
        
        if (use_only_vdj) {
            db_clonal = define_tcr_clone_helper(db=db_vjl_formatted, 
                                                col_vjl=col_vjl_group, 
                                                col_junc_1=junc, col_junc_2=NULL, 
                                                col_clone_id=clone_id)
        } else {
            db_clonal = define_tcr_clone_helper(db=db_vjl_formatted, 
                                                col_vjl=col_vjl_group, 
                                                col_junc_1=col_junc_vdj, col_junc_2=col_junc_vj, 
                                                col_clone_id=clone_id)
        }
        
    } else {
        db_clonal = define_tcr_clone_helper(db=db_vjl_formatted, 
                                            col_vjl=col_vjl_group, 
                                            col_junc_1=junc, col_junc_2=NULL, 
                                            col_clone_id=clone_id)
    }
    
    
    ### propagate clone IDs, if necessary
    
    cat("\nPropagating clone IDs...\n")
    db[[col_vjl_group]] = NA
    db[[clone_id]] = NA
    
    if (single_cell_mode) {
        # wrt db_clonal
        idx_cell_id_db_clonal = match(db[[cell_id]], db_clonal[[cell_id]])
        stopifnot(!any(is.na(idx_cell_id_db_clonal)))
        stopifnot(all.equal( db[[cell_id]], 
                             db_clonal[[cell_id]][idx_cell_id_db_clonal] ))
        
        db[[col_vjl_group]] = db_clonal[[col_vjl_group]][idx_cell_id_db_clonal]
        db[[clone_id]] = db_clonal[[clone_id]][idx_cell_id_db_clonal]
        
    } else {
        db[[col_vjl_group]] = db_clonal[[col_vjl_group]]
        db[[clone_id]] = db_clonal[[clone_id]]
    }
    stopifnot(!any(is.na(db[[col_vjl_group]])))
    stopifnot(!any(is.na(db[[clone_id]])))
    
    cat("\nOutput db has", nrow(db), "seqs, representing",
        length(unique(db[[clone_id]])), "clones\n")
    
    return(list(db_passed=db, db_failed=db_failed))
}


#' Infer TCR clonal relationships based on VJL grouping and identical junctions
#' 
#' @param  db            A `data.frame` containing `col_vjl`, and at least one of
#'                       `col_junc_1` and `col_junc_2`. See below for details.
#' @param  col_vjl       Column name for VJL grouping.
#'                       Defaults to `"vjl_group"`.
#' @param  col_junc_1    Column name for junction 1 (presumably those of VDJ seqs).
#'                       Defaults to `NULL`.
#' @param  col_junc_2    Column name for junction 2 (presumably those of VJ seqs).
#'                       Defaults to `NULL`.
#' @param  col_clone_id  Column name for clone ID (to be added).
#'                       Defaults to `"clone_id"`.
#' 
#' @returns A `data.frame` with an added column named `col_clone_id`.
#' 
#' @details This function assums that prior VJL grouping (based on V gene, J gene, 
#'          and junction length) has been performed and stored in the `col_vjl` column.
#'          
#'          Within each VJL group, the function performs further grouping based on one 
#'          or both of `col_junc_1` and `col_junc_2` to define clonal relationships.
#' 
#'          At least one of `col_junc_1` and `col_junc_2` must be specified. 
#'          If both specified, both junction sequences will be used simultaneously
#'          for clonal grouping.
#'          
#'          Notes on input `db`:
#'          
#'          1) If inferring clonal relationships using both VDJ (TRB/D) and VJ (TRA/G)
#'          sequences, `db` is expected to of "wide" format, meaning that each row
#'          should correspond to a cell, with both junctions from the VDJ and the VJ
#'          sequences stored in the same row.
#'          
#'          2) `db` should only contain sequences actually used in clonal inference.
#'          Specifically, if data has both VDJ and VJ seqs, but only VDJ seqs are 
#'          used for inferring clonal relationships, then `db` should only include 
#'          those VDJ seqs. Propagation of the inferred clonal relationships to VJ seqs 
#'          not used for the inference process is outside the scope of this function.
#'          
define_tcr_clone_helper = function(db, col_vjl="vjl_group", 
                                   col_junc_1=NULL, col_junc_2=NULL, 
                                   col_clone_id="clone_id") {
    
    require(dplyr)
    
    ### column checks
    
    # at least one of col_junc_ must be specified
    bool_null_junc_1 = is.null(col_junc_1)
    bool_null_junc_2 = is.null(col_junc_2)
    stopifnot(!bool_null_junc_1 | !bool_null_junc_2)
    
    if (!bool_null_junc_1) { 
        stopifnot(col_junc_1 %in% colnames(db)) 
        stopifnot(!any(is.na(db[[col_junc_1]])))
        stopifnot(!any(db[[col_junc_1]]==""))
    }
    
    if (!bool_null_junc_2) { 
        stopifnot(col_junc_2 %in% colnames(db)) 
        stopifnot(!any(is.na(db[[col_junc_2]])))
        stopifnot(!any(db[[col_junc_2]]==""))
    }
    
    stopifnot(col_vjl %in% colnames(db))
    stopifnot(!any(is.na(db[[col_vjl]])))
    stopifnot(!any(db[[col_vjl]]==""))
    
    if (col_clone_id %in% colnames(db)) {
        warning(col_clone_id, " already exists as a column; it will be overwritten.")
    }
    
    
    ### infer clonal relationships
    
    if (!bool_null_junc_1 & !bool_null_junc_2) {
        cat("\nInferring TCR clones based on", col_vjl, col_junc_1, col_junc_2, "\n")
        db_grouped = group_by(.data=db, 
                              !!rlang::sym(col_vjl), 
                              !!rlang::sym(col_junc_1), !!rlang::sym(col_junc_2))
        
    } else if (!bool_null_junc_1 & bool_null_junc_2) {
        cat("\nInferring TCR clones based on", col_vjl, col_junc_1, "\n")
        db_grouped = group_by(.data=db, 
                              !!rlang::sym(col_vjl), 
                              !!rlang::sym(col_junc_1))
        
    } else if (bool_null_junc_1 & !bool_null_junc_2) {
        cat("Inferring TCR clones based on", col_vjl, col_junc_2, "\n")
        db_grouped = group_by(.data=db, 
                              !!rlang::sym(col_vjl), 
                              !!rlang::sym(col_junc_2))
        
    }
    
    #group_indices(db_grouped)
    #group_size(db_grouped)
    #n_groups(db_grouped)
    
    
    ### generate clone IDs
    
    tbl_clone_def = group_keys(db_grouped)
    
    tbl_clone_def[[col_clone_id]] = NA
    
    vec_uniq_vjl = unique(tbl_clone_def[[col_vjl]])
    
    for (cur_vjl in vec_uniq_vjl) {
        cur_vjl_idx = which(tbl_clone_def[[col_vjl]]==cur_vjl)
        tbl_clone_def[[col_clone_id]][cur_vjl_idx] = paste0(cur_vjl, "_", 1:length(cur_vjl_idx))
    }
    stopifnot(!any(is.na(tbl_clone_def[[col_clone_id]])))
    
    
    ### distribute clone IDs to rows
    lst_group_rows = group_rows(db_grouped)
    db[[col_clone_id]] = NA
    
    for (i_cl in 1:length(lst_group_rows)) {
        cur_row_idx = lst_group_rows[[i_cl]]
        db[[col_clone_id]][cur_row_idx] = tbl_clone_def[[col_clone_id]][i_cl]
    }
    stopifnot(!any(is.na(db[[col_clone_id]])))
    
    return(db)
}

