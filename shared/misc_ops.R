# miscellaneous operational tasks

#' Select rows with specific values in specific columns
#'
#' @param  db       A `data.frame`.
#' @param  col_vec  A `vector` containing column names.
#' @param  val_lst  A `list` containing column values.
#' 
#' @return  A `boolean` vector matching `nrow(db)`.
#' 
#' @details There must be a one-to-one match between the entries in `col_vec` and
#'          the entries in `val_lst`. 
#'          Each entry in `col_vec` must be a single string.
#'          Each entry in `val_lst` must be one or more values.
#' 
db_select = function(db, col_vec, val_lst) {
   # checks
    stopifnot( length(col_vec) == length(val_lst) )
    stopifnot( all( col_vec %in% colnames(db) ) )
    
    bool_lst = sapply(1:length(col_vec), 
                      function(i) {
                          # a single column
                          cur_col = col_vec[i]
                          # one or more values in that column
                          cur_val = val_lst[[i]]
                          cur_bool = db[[cur_col]] %in% cur_val
                          return(cur_bool)
                      },
                      USE.NAMES=F, simplify=F)
    
    # rows match rows in db
    # cols match col_vec
    bool_mtx = do.call(cbind, bool_lst)
    
    # AND across col_vec
    bool = rowSums(bool_mtx)==length(col_vec)
    
    return(bool)
}

