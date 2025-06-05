# Julian Q. Zhou
# https://github.com/julianqz

# miscellaneous operational tasks

#' Compute full-range table
#' 
#' @param   tab          A `table` returned by `table()` (or a named numeric 
#'                       vector with integers as names).
#' @param   range_lower  An integer specifying the lower bound of the full range.
#'                       Default is `NULL`.
#' @param   range_upper  An integer specifying the upper bound of the full range.
#'                       Default is `NULL`.
#' @param   range_by     An integer specifying the interval between values in the
#'                       full range. Default is `1`.
#' @param   NA_for_zero  Whether to set a value in the full range for which the 
#'                       the count is 0 to `NA`. Default is `TRUE`.
#' 
#' @return  A named numeric vector with integers as names. The names cover the
#'          specified full range.
#' 
#' @details If left as `NULL`, `range_lower` is set to the lower bound of `tab`,
#'          and `range_upper` is set to the upper bound of `tab`.

get_full_range_table = function(tab, range_lower=NULL, range_upper=NULL, 
                                range_by=1, NA_for_zero=TRUE) {
    
    tab_range = range(as.numeric(names(tab)), na.rm=T)
    if (is.null(range_lower)) { range_lower=tab_range[1] }
    if (is.null(range_upper)) { range_upper=tab_range[2] }
    stopifnot(range_lower<=tab_range[1])
    stopifnot(range_upper>=tab_range[2])
    
    tab_full_range = seq(from=range_lower, to=range_upper, by=range_by)
    
    # initialize
    if (NA_for_zero) {
        tab_full = rep(NA, length(tab_full_range))
    } else {
        tab_full = rep(0, length(tab_full_range))
    }
    names(tab_full) = as.character(tab_full_range)
    
    idx = match(names(tab), names(tab_full))
    stopifnot(!any(is.na(idx)))
    stopifnot(all.equal( names(tab_full)[idx], names(tab) ))
    
    tab_full[idx] = tab
    
    return(tab_full)
}


#' Convert to numeric values to a specified scale
#' 
#' @param vec_orig    A numeric vector of original values.
#' @param max_orig    If specified, use this as the upper bound of original values.
#'                    If `NULL`, the upper bound will be derived from `vec_orig`.
#' @param min_orig    If specified, use this as the lower bound of original values.                    
#'                    If `NULL`, the lower bound will be derived from `vec_orig`.
#' @param max_scale   Maximum in new scale.
#' @param min_scale   Minimum in new scale.
#' @param keep_NA     Whether to keep NAs, if any, in `vec_orig`. Boolean.
#' 
#' @return  A numeric vector of scaled values.
#' 
#' @details If `keep_NA` is `TRUE`, any `NA` in `vec_orig` will remain `NA`.
#'          If `keep_NA` is `FALSE`, any `NA` in `vec_orig` will become `min_scale`.
#' 
scale_values = function(vec_orig, max_orig=NULL, min_orig=NULL,
                        max_scale, min_scale, keep_NA=FALSE) {
    
    # if max_orig or min_orig is not specified, derive from vec_orig
    if (is.null(max_orig)) { max_orig = max(vec_orig, na.rm=T) }
    if (is.null(min_orig)) { min_orig = min(vec_orig, na.rm=T) }
    
    range_orig = max_orig - min_orig
    range_scale = max_scale - min_scale
    
    vec_new = sapply(vec_orig, function(v){
        if (is.na(v)) {
            if (keep_NA) {
                return(NA)
            } else {
                # NA in vec_orig assigned min scaled value
                return(min_scale)
            }
        } else {
            if (range_orig>0) {
                # if there's variation in vec_orig entries
                v_2 = (v-min_orig) / range_orig * range_scale + min_scale
                return(v_2)
            } else {
                # if all vec_orig entries are the same
                return(min_scale)
            }
        }
    }, USE.NAMES=F)
    
    # sanity checks
    # length should stay unchanged
    stopifnot(length(vec_orig)==length(vec_new))
    if (keep_NA) {
        # if keeping NAs, the number of NAs should stay unchanged
        stopifnot( sum(is.na(vec_orig)) == sum(is.na(vec_new)) )
    } else {
        # if setting NAs to min, there should be no more NAs
        stopifnot(!any(is.na(vec_new)))
    }
    
    return(vec_new)
}


#' Select rows with specific values in specific columns
#'
#' @param  df       A `data.frame`.
#' @param  col_vec  A `vector` containing column names.
#' @param  val_lst  A `list` containing column values.
#' @param  rev_vec  A `vector` containing boolean values.
#' 
#' @return  A `boolean` vector matching `nrow(df)`.
#' 
#' @details There must be a one-to-one match between the entries in `col_vec`,
#'          `val_lst`, and `rev_vec`.
#'          Each entry in `col_vec` must be a single string.
#'          Each entry in `val_lst` must be one or more values.
#'          Each entry in `rev_vec` must be a single boolean value.
#'          
#'          For a column in `col_vec`, with `FALSE` in `rev_vec`, `%in%` is tested
#'          against value(s) in `val_lst`. With `TRUE` in `rev_vec`, `! %in%` is
#'          tested. 
#' 
#' @exapmles
#' # clone_info[["mab"]]!=0 & clone_info[["bulk"]]!=0 & clone_info[["tgx"]]==0
#' df_select(clone_info, 
#'           col_vec=c("mab", "bulk", "tgx"),
#'           val_lst=list(c(0), c(0), c(0)),
#'           rev_vec=c(T, T, F))
#' 
df_select = function(df, col_vec, val_lst, rev_vec) {
   # checks
    stopifnot( length(col_vec) == length(val_lst) )
    stopifnot( length(col_vec) == length(rev_vec) )
    stopifnot( all( col_vec %in% colnames(df) ) )
    stopifnot( all( rev_vec %in% c(TRUE, FALSE) ) )
    
    bool_lst = sapply(1:length(col_vec), 
                      function(i) {
                          # a single column
                          cur_col = col_vec[i]
                          # one or more values in that column
                          cur_val = val_lst[[i]]
                          cur_bool = df[[cur_col]] %in% cur_val
                          # reverse?
                          if (rev_vec[i]) {
                              cur_bool = !cur_bool
                          }
                          return(cur_bool)
                      },
                      USE.NAMES=F, simplify=F)
    
    # rows match rows in df
    # cols match col_vec
    bool_mtx = do.call(cbind, bool_lst)
    
    # AND across col_vec
    bool = rowSums(bool_mtx)==length(col_vec)
    
    return(bool)
}


#' Read in and combine a list of files 
#' 
#' @param  vec_path   char vector contaning paths to files
#' @param  file_type  file type; currently supported value: "tsv"
#' 
#' @returns  A `data.frame`.
#' 
#' @details  Any file that does not exist is skipped.
#'           Assumes files all have the same colnames (otherwise rbind won't work)
#' 
load_all_files = function(vec_path, file_type) {
    
    n_path = length(vec_path)
    lst_db = vector(mode="list", length=n_path)
    
    for (i in 1:n_path) {
        
        cur_path = vec_path[i]
        if (file.exists(cur_path)) {
            
            if (file_type=="tsv") {
                cur_db = read.table(cur_path, sep="\t", header=T, stringsAsFactors=F)
            } else {
                warning("\nUnsupported file_type:", file_type, "\n")
                cur_db = NULL
            }
            
            lst_db[[i]] = cur_db
            rm(cur_db)
        } else {
            cat("\nNot found:", cur_path, "\n")
        }
    }
    
    db_all = do.call(rbind, lst_db)
    return(db_all)
}

