#' Summarize clonal compositions
#' 
#' @param   db             A `data.frame`.
#' @param   col_clone      Column name containing clone ID.
#' @param   col_vec        A `vector` containing column names.
#' @param   order_by_size  Whether and how to order by clone size. 
#'                         One of `none`, `increasing`, or `decreasing`.
#' 
#' @return  A `data.frame` containing summary info on the clones in `db[[col_clone]]`.
#' 
#' @details Columns in the returned `data.frame` include `col_clone`, `clone_size`,
#'          and one column each for each unique value in each column in `col_vec`.
#'          
#'          Values in `col_vec` must not start with non-alphabet characters as 
#'          these values would become colnames of the returned `data.frame`, 
#'          which would not allow non-alphabet characters as the first character.
#'          
#'          If any `col_vec` has values of `""` (empty) or `NA`, these values
#'          will be converted to `"[col_vec]_na"` and recorded as such in the 
#'          returned `data.frame`. 
#'          For example, `""` and `NA` in the column `isotype` will be counted in
#'          the `isotype_na` column.
#'
#' @examples 
#' summarize_clone(db, col_clone="clone_id",
#'                 col_vec=c("seq_type", "timepoint", "tissue", "sorting", "isotype"),
#'                 order_by_size="decreasing")
#' 
summarize_clone = function(db, col_clone, col_vec, 
                           order_by_size=c("none", "increasing", "decreasing")) {

    # all cols present in db
    stopifnot( all( c(col_clone, col_vec) %in% colnames(db) ) )
    stopifnot( order_by_size %in% c("none", "increasing", "decreasing") )
    
    # if there's any "" or NA in any column in col_vec
    # convert that value to "na"
    for (col in col_vec) {
        
        cur_bool_NA = is.na(db[[col]])
        cur_bool_empty = db[[col]]==""
        # using || does not appear to work
        cur_bool = cur_bool_NA | cur_bool_empty
        
        stopifnot( !any(is.na(cur_bool)) )
        stopifnot( sum(cur_bool) == sum(cur_bool_empty, na.rm=T)+sum(cur_bool_NA) )
        
        if (any(cur_bool)) {
            # append prefix to distinguish potential multiple 'na' 
            # arising from different columns
            cur_na = paste0(col, "_na")
            cat(col, ": padding", sum(cur_bool_empty, na.rm=T), "''",
                "and", sum(cur_bool_NA), "NA as '", 
                cur_na, "'.\n")
            db[[col]][cur_bool] = cur_na
        }
    }
    
    # no overlap amongst unique values from col_vec
    val_uniq_lst = sapply(col_vec, function(col){ sort(unique(db[[col]])) }, 
                          USE.NAMES=T, simplify=F)
    val_uniq_vec = unlist(val_uniq_lst)
    stopifnot( length(unique(val_uniq_vec)) == length(val_uniq_vec) )
    
    # initiate clone_info
    
    clones = unique(db[[col_clone]])
    n_clones = length(clones)
    
    # +2 to include col_clone and clone_size
    clone_info = data.frame(matrix(NA, ncol=(length(val_uniq_vec)+2),
                                   nrow=n_clones))
    col_size = "clone_size"
    colnames(clone_info) = c(col_clone, col_size, val_uniq_vec)
    
    # if any val_uniq_vec starts with non-alphabet character, would cause
    # discrepency in colnames(clone_info) and val_uniq_vec
    if (any(!grepl(pattern="^[[:alpha:]]+", x=val_uniq_vec))) {
        stop("All values in col_vec must start with an alphebtic character.")
    }
    
    # fill in clone_info
    
    for (i_cl in 1:n_clones) {
        cl = clones[i_cl]
        # wrt db
        idx = which(db[[col_clone]]==cl)
        
        clone_info[[col_clone]][i_cl] = cl
        clone_info[[col_size]][i_cl] = length(idx)
        
        for (col in col_vec) {
            for (val in val_uniq_lst[[col]]) {
                clone_info[[val]][i_cl] = sum( db[[col]][idx]==val )
            }
        }
    }
    
    # sanity checks
    
    # no NA anywhere
    stopifnot(!any(is.na(clone_info)))
    
    # the math should add up
    for (col in col_vec) {
        cur_val = val_uniq_lst[[col]]
        if (length(cur_val)==1) {
            stopifnot( all.equal( clone_info[[col_size]], clone_info[[cur_val]] ) )
        } else {
            stopifnot( all.equal( clone_info[[col_size]], 
                                  rowSums(clone_info[, cur_val]) ) )
        }
    }
    
    # ordering by clone size?
    if (order_by_size=="increasing") {
        clone_info = clone_info[order(clone_info[[col_size]], decreasing=F), ]
    } else if (order_by_size=="decreasing") {
        clone_info = clone_info[order(clone_info[[col_size]], decreasing=T), ]
    }
    
    return(clone_info)
}

