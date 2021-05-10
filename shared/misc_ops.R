# miscellaneous operational tasks

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

