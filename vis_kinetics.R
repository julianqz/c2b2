# Visualization of temporal kinetics of B cell clones

#' get the onset timepoint
#' 
#' @param  df          A `data.frame` containing sequence/cell counts at various
#'                     timepoints for each clone. Each row is a clone. Columns
#'                     represent timepoints.
#' @param  vec_col_tp  A `character` vector containing column names representing
#'                     timepoints. All entries must be in `df`.
#'
#' @return  A list containing two vectors.
#'          - `vec_onset`: A `character` vector with length matching `nrow(df)` 
#'          in which each item matches a value in `vec_col_tp` and indicates 
#'          the onset timepoint of a clone.
#'          - `vec_order_by_onset`: An integer vector containing reordering index
#'          for rows in `df` based on onset timepoint
#'         
#' @details The order of timepoints specified in `vec_col_tp` matters and is 
#'          assumed to represent earlier to later timepoints from left to right.  
#' 
get_onset = function(df, vec_col_tp) {
    
    stopifnot(all(vec_col_tp %in% colnames(df)))
    stopifnot(length(vec_col_tp)>1)
    bin_mtx = df[, vec_col_tp]>0
    
    # each list item corresponds to a timepoint and is a boolean vector
    # TRUE if a clone first detected at a timepoint
    lst_bool = vector(mode="list", length=length(vec_col_tp))
    names(lst_bool) = vec_col_tp
    for (i in 1:length(vec_col_tp)) {
        if (i==1) {
            # first detected at first timepoint in vec_col_tp
            cur_bool = bin_mtx[, i]
        } else if (i==2) {
            # first detected at second timepoint in vec_col_tp
            # not detected at first timepoint
            # rowSums of a one-column matrix (1-D vector) throws error
            cur_bool = bin_mtx[, i] & !bin_mtx[, 1]
        } else {
            # first detected at i-th timepoint in vec_col_tp, i>=3
            # not detected at all timepoints before i-th
            cur_bool = bin_mtx[, i] & rowSums(bin_mtx[, 1:(i-1)])==0
        }
        lst_bool[[i]] = cur_bool
    }
    
    # convert list to df
    # each row corresponds to a clone
    # each col corresponds to a timepoint
    df_bool = do.call(cbind, lst_bool)
    # sanity check
    # there should be only 1 TRUE per row (each clone has only 1 onset timepoint)
    stopifnot(all(rowSums(df_bool)==1))
    
    # column index of onset timepoints
    vec_index = apply(df_bool, 1, function(x){which(x)})
    # sanity check
    # no NA
    stopifnot(!any(is.na(vec_index)))
    
    # vector of onset timepoints
    # each item corresponds to a clone
    vec_onset = vec_col_tp[vec_index]
    
    vec_order_by_onset = order(vec_index)
    
    return(list(vec_onset=vec_onset,
                vec_order_by_onset=vec_order_by_onset))
}


#' prepare a data.frame for visualizing temporal kinetics of clones
#' 
#' @param  df           A `data.frame` containing sequence/cell counts at various
#'                      timepoints for each clone. Each row is a clone. Columns
#'                      represent timepoints. Additional columns, such as clone 
#'                      ID, are allowed and will be propogated as is.
#' @param  vec_col_tp   A `character` vector containing column names representing
#'                      timepoints. All entries must be in `df`. Order matters 
#'                      (see Details).
#' @param  single_last  Boolean. If `TRUE`, place clones detected only at a single
#'                      timepoint last in the returned `df`. Defaults to `TRUE`.
#'                                          
#' @return  A modified `data.frame` based on `df`. 
#' 
#'          A column named `onset` is added. This column indicates the first 
#'          timepoint at which a clone was first detected.
#'          
#'          A column named `single` is added. This column indicates whether a 
#'          clone was detected only at a single timepoint.
#' 
#'          The total number of rows remains unchanged. Rows/clones are reordered 
#'          in the order of their onset timepoints. That is, clones with earlier 
#'          onset timepoints appear before those with later onset timepoints in 
#'          the returned `data.frame`. 
#'          
#'          If `single_last` is `TRUE`, rows/clones are further reordered such 
#'          that clones detected only at a single timepoint appear after clones 
#'          detected at more than one timepoint. If this is performed, the order
#'          previously introduced by reordering based on onset timepoint is
#'          preserved within the single-timepoint clones.
#'                    
#'          After all reordering, columns representing percentages and cumulative 
#'          percentages are added. The percentage of a clone at a timepoint is 
#'          calculated out of the total number of sequences/cells across clones 
#'          at that timepoint. The cumulative percentages are calculated according 
#'          to the row order of the returned `data.frame`.
#' 
#' @details The order of timepoints specified in `vec_col_tp` matters and is 
#'          assumed to represent earlier to later timepoints from left to right.
#'          
prep_kinetics_df = function(df, vec_col_tp, single_last=TRUE) {
    
    # at least 2 timepoints
    stopifnot(length(vec_col_tp)>=2)
    # all in df
    stopifnot(all(vec_col_tp %in% colnames(df)))
    
    # colname "single" is reserved
    stopifnot(!"single" %in% colnames(df))
    # colname "onset" is reserved
    stopifnot(!"onset" %in% colnames(df))
    
    # all counts >=0
    stopifnot(all(df[, vec_col_tp]>=0))
    
    # calculate onset timepoint
    lst_onset = get_onset(df, vec_col_tp)
    df[["onset"]] = lst_onset[["vec_onset"]]
    
    # reorder by onset timepoint
    df = df[lst_onset[["vec_order_by_onset"]], ]
    
    # binary matrix, whether cell is 0 (TRUE) or not (FALSE)
    df_bin_0 = df[, vec_col_tp]==0
    # no clone with all 0 counts
    stopifnot(!any(rowSums(df_bin_0)==length(vec_col_tp)))
    
    # clones with only 1 non-zero timepoint
    bool_single_tp = rowSums(df_bin_0)==(length(vec_col_tp)-1)
    
    # there must be at least 1 clone with more than 1 timepoint
    # (otherwise not much point visualizing)
    stopifnot(sum(!bool_single_tp)>=1)
    
    df[["single"]] = bool_single_tp
    
    # place clones detected only at a single timepoint last
    if (single_last) {
        # TRUE's (single tp) come after FALSE's
        df = df[order(bool_single_tp), ]
    }

    # IMPORTANT: cumulative percentage calculation should be performed after
    # ordering of `df` rows has been finalized
    # (Any further reordering would require recalculation of cumulative %)
    
    # data.frame to hold calculated statistics
    
    # percentage
    df_add_cols_perc = paste0(vec_col_tp, "_perc")
    # cumulative percentage
    df_add_cols_perc_cumu = paste0(vec_col_tp, "_perc_cumu")
    
    df_add_cols = c(df_add_cols_perc, df_add_cols_perc_cumu)
    
    df_add = data.frame(matrix(NA, nrow=nrow(df), ncol=length(df_add_cols)))
    colnames(df_add) = df_add_cols
    
    for (i in 1:length(vec_col_tp)) {
        cur_tp = vec_col_tp[i]
        
        # percentage columns are not essential
        # but useful in case of debugging
        cur_tp_perc = df_add_cols_perc[i]
        df_add[, cur_tp_perc] = df[[cur_tp]]/sum(df[[cur_tp]])*100
        
        # cumulative percentage columns are essential
        cur_tp_perc_cumu = df_add_cols_perc_cumu[i]
        df_add[, cur_tp_perc_cumu] = cumsum(df[[cur_tp]])/sum(df[[cur_tp]])*100
    }
    
    df = cbind(df, df_add)
    
    return(df)
}


#' visualize temporal kinetics of clones based on output from `prep_kinetics_df`
#' 
#' @param  df                 A `data.frame` returned by `prep_kinetics_df`.
#' @param  vec_col_perc_cumu  A `character` vector containing colnames for columns
#'                            containing sequence/cell counts (in cumulative 
#'                            percentages) at various timepoints.
#' @param  ymax               Upper limit of y-axis to be drawn. Expressed in `%` 
#'                            (eg. a `100` would mean `100%`). Defaults to `NULL`
#'                            in which case the upper limit is inferred from data.                           
#' @param  col_border         Passed to `polygon()` for controlling the color of
#'                            the border drawn for each clone's polygon.
#' @param  lwd_border         Passed to `polygon()` for controlling the width of
#'                            the border drawn for each clone's polygon.                           
#' @param  vec_color          Vector specifying color for each clone. Expects
#'                            1-to-1 correspondence with row order in `df`.
#' @param  vec_label_tp       X-axis text labels to be displayed for timepoints.
#'                            Expects 1-to-1 correspondence with `vec_col_perc_cumu`.
#' @param  cex_scale_label    `cex` passed to `text()` for controlling the size of 
#'                             scale bar text label.
#' 
#' @details  It is expected that all entries in `vec_col_perc_cumu` exist as
#'           columns in `df`, as well as a boolean column called `single`.
#'           
#'           If `show_single` is `FALSE` (default), clones detected at only a 
#'           single timepoint are not shown. If `TRUE`, they are shown. 
#'           
#'           The order in which the clones are drawn, from bottom to top along the
#'           vertical axis, is the same as the row order of `df`.
#' 
plot_kinetics = function(df, vec_col_perc_cumu, ymax=NULL,
                         col_border=NA, lwd_border=1, vec_color,
                         vec_label_tp, cex_scale_label=1) {
    
    stopifnot(all(vec_col_perc_cumu %in% colnames(df)))
    stopifnot(length(vec_color)==nrow(df))
    stopifnot(length(vec_col_perc_cumu)==length(vec_label_tp))
    
    n_tp = length(vec_col_perc_cumu)
    
    # x-coords of polygons
    # clockwise
    poly_x = c(1:n_tp, # top of polygon, from L to R
               n_tp:1) # bottom of polygon, from R to L
    
    # if unspecified, calculate ymax from data
    if (is.null(ymax)) {
        ymax = max(df[, vec_col_perc_cumu])
    }
    
    # initialize
    ymax_extend = 1.01
    
    plot(x=1:n_tp, y=rep(NA, n_tp), 
         xlim=c(0.5, n_tp),
         ylim=c(0, ymax*ymax_extend), 
         type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
    
    # x-axis
    axis(side=1, at=1:n_tp, labels=vec_label_tp, tick=F)
    
    # scale bar?
    # 10 times smaller than the closest multiple of 5
    
    # calculate from data in case this is different from manually specified `ymax`
    ymax_actual = max(df[, vec_col_perc_cumu])
    scale_bar_len = ceiling(ymax_actual/5)*5/10
    scale_bar_y_mid = ymax/2
    scale_bar_y_lb = scale_bar_y_mid-scale_bar_len/2
    scale_bar_y_ub = scale_bar_y_mid+scale_bar_len/2
    
    segments(x0=0.85, x1=0.85, y0=scale_bar_y_lb, y1=scale_bar_y_ub)
    text(x=0.85, y=scale_bar_y_lb*0.9, cex=cex_scale_label,
         labels=as.character(paste0(scale_bar_len, "%")))
    
    # one polygon per clone
    for (i in 1:nrow(df)) {
        
        # i_col = ifelse(i%%2==0, "blue", "red") # for debugging
        i_col = vec_color[i]
        
        if (i==1) {
            # first clone
            # use sapply to get a vector
            # c(df[i, vec_col_perc_cumu], ...) returns a list 
            polygon(x=poly_x,
                    y=c(sapply(vec_col_perc_cumu, function(s){df[[s]][i]}, USE.NAMES=F), 
                        rep(0, n_tp)),
                    col=i_col, border=col_border, lwd=lwd_border)
            
        } else {
            # subsequence clones
            # rev() for bottom side because direction there is from R to L
            polygon(poly_x,
                    c(sapply(vec_col_perc_cumu, function(s){df[[s]][i]}, USE.NAMES=F),
                      sapply(rev(vec_col_perc_cumu), function(s){df[[s]][i-1]}, USE.NAMES=F)
                      ),
                    col=i_col, border=col_border, lwd=lwd_border)
            
        }
    }

}

