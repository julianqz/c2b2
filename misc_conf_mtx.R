# confusion matrix for gex-based cell identity annotations

#' compute a confusion matrix for cell identity annotations
#' 
#' @param  df_1           A `data.frame`.
#' @param  df_2           A `data.frame`.
#' @param  col_cell_id_1  Column name for column in `df_1` containing cell IDs.
#' @param  col_cell_id_2  Column name for column in `df_2` containing cell IDs.
#' @param  col_anno_1     Column name for column in `df_1` containing cell 
#'                        identity annotations.
#' @param  col_anno_2     Column name for column in `df_2` containing cell 
#'                        identity annotations.
#'
#' @return A confusion matrix of class `table`.
#' 
#' @details
#' It is assumed and checked that `col_anno_1` and `col_anno_2` do not contain
#' character values that are `"na"` or `"absent"`.
#' 
#' Any values that are `NA`, `""` (empty character), or `NaN` in `col_anno_1` 
#' and `col_anno_2` will be converted to the character value `"na"`.
#' 
#' Cells that are present in one `data.frame` but missing in the other will be
#' counted under the `"absent"` category.
#' 
get_conf_mtx = function(df_1, df_2, 
                        col_cell_id_1, col_cell_id_2,
                        col_anno_1, col_anno_2) {
    
    stopifnot(nrow(df_1)>0 & nrow(df_2)>0)
    
    # check that columns exist
    stopifnot(col_cell_id_1 %in% colnames(df_1))
    stopifnot(col_cell_id_2 %in% colnames(df_2))
    stopifnot(col_anno_1 %in% colnames(df_1))
    stopifnot(col_anno_2 %in% colnames(df_2))
    
    # check that all cell IDs unique
    stopifnot(!any(duplicated(df_1[[col_cell_id_1]])))
    stopifnot(!any(duplicated(df_2[[col_cell_id_2]])))
    
    # assumes that col_anno does not contain values that are "na" (character) or "absent"
    # check that
    vec_forbidden = c("na", "absent")
    stopifnot(!any(df_1[[col_anno_1]] %in% vec_forbidden, na.rm=T))
    stopifnot(!any(df_2[[col_anno_2]] %in% vec_forbidden, na.rm=T))
    
    # if there's any NA or "" or NaN in col_anno, convert to "na" (character)
    bool_na_1 = is.na(df_1[[col_anno_1]]) | df_1[[col_anno_1]]=="" | is.nan(df_1[[col_anno_1]])
    bool_na_2 = is.na(df_2[[col_anno_2]]) | df_2[[col_anno_2]]=="" | is.nan(df_2[[col_anno_2]])
    if (any(bool_na_1)) { df_1[[col_anno_1]][bool_na_1] = "na" }
    if (any(bool_na_2)) { df_2[[col_anno_2]][bool_na_2] = "na" }
    # no NA
    stopifnot(!any(is.na(df_1[[col_anno_1]])))
    stopifnot(!any(is.na(df_2[[col_anno_2]])))
    # no ""
    # NA=="" gives NA ==> use na.rm=T
    stopifnot(!any(df_1[[col_anno_1]]=="", na.rm=T))
    stopifnot(!any(df_2[[col_anno_2]]=="", na.rm=T))
    # no NaN
    # is.nan(NA) gives FALSE ==> no need for na.rm=T
    stopifnot(!any(is.nan(df_1[[col_anno_1]])))
    stopifnot(!any(is.nan(df_2[[col_anno_2]])))
    
    # cell id's
    vec_cell_id_union = union(df_1[[col_cell_id_1]], df_2[[col_cell_id_2]])
    # wrt df_1
    idx_cell_id_1 = match(vec_cell_id_union, df_1[[col_cell_id_1]])
    bool_not_na_idx_cell_id_1 = !is.na(idx_cell_id_1)
    stopifnot(all.equal( vec_cell_id_union[bool_not_na_idx_cell_id_1],
                         df_1[[col_cell_id_1]][ idx_cell_id_1[bool_not_na_idx_cell_id_1] ] ))
    # wrt df_2
    idx_cell_id_2 = match(vec_cell_id_union, df_2[[col_cell_id_2]])
    bool_not_na_idx_cell_id_2 = !is.na(idx_cell_id_2)
    stopifnot(all.equal( vec_cell_id_union[bool_not_na_idx_cell_id_2],
                         df_2[[col_cell_id_2]][ idx_cell_id_2[bool_not_na_idx_cell_id_2] ] ))
    
    # holder
    col_mtx_cell_id = "cell_id"
    col_mtx_anno_1 = "anno_1"
    col_mtx_anno_2 = "anno_2"
    df_union_cols = c(col_mtx_cell_id, col_mtx_anno_1, col_mtx_anno_2)
    df_union = data.frame(matrix("absent", 
                                 nrow=length(vec_cell_id_union),
                                 ncol=length(df_union_cols)))
    colnames(df_union) = df_union_cols
    
    # fill in
    df_union[[col_mtx_cell_id]] = vec_cell_id_union
    df_union[[col_mtx_anno_1]][bool_not_na_idx_cell_id_1] = df_1[[col_anno_1]][ idx_cell_id_1[bool_not_na_idx_cell_id_1] ]
    df_union[[col_mtx_anno_2]][bool_not_na_idx_cell_id_2] = df_2[[col_anno_2]][ idx_cell_id_2[bool_not_na_idx_cell_id_2] ]
    stopifnot( sum(df_union[[col_mtx_anno_1]]!="absent") == nrow(df_1) )
    stopifnot( sum(df_union[[col_mtx_anno_2]]!="absent") == nrow(df_2) )
    # "absent" in col_mtx_anno_[12] means cell not present
    # there shouldn't be any cell that's absent in both df_1 and df_2
    stopifnot( !any( rowSums(df_union[, c(col_mtx_anno_1, col_mtx_anno_2)]=="absent")==2 ) )
    # no NA 
    stopifnot(!any(is.na(df_union)))
    
    mtx = table(df_union[[col_mtx_anno_1]], df_union[[col_mtx_anno_2]], 
                dnn=list("df_1", "df_2"), useNA="no")
    # sanity check
    stopifnot( sum(mtx[rownames(mtx)!="absent", ]) == nrow(df_1) )
    stopifnot( sum(mtx[, colnames(mtx)!="absent"]) == nrow(df_2) )
    
    return(mtx)
}

