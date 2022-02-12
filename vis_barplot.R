# Julian Q. Zhou
# https://github.com/julianqz

#' Enhanced barplot
#' 
#' @param  tab         Input counts. Either as a 1-D vector, or a 2-D matrix.
#'                     See details. 
#' @param  vec_ylim    Y-axis range limits. Passed on to `barplot`.
#' @param  vec_labels  Text labels on top of the bars. Usually counts of sort.
#'                     Optional. Default is `NULL`.
#'                     To specify, supply a character vector whose length matches 
#'                     the length (if 1-D vector) or nrow (if 2-D matrix) of `tab`.
#' @param  cex.text    Size of text labels. Passed on to `cex` in `text`.
#'                     Default is `1`. Ignored if `vec_labels` is `NULL`.                     
#' @param  las         Orientation of labels. Passed on to `barplot`.
#'                     Default is `2`.
#' @param  ylab        Y-axis label. Passed on to `barplot`. Default is `""`.
#' @param  cex.lab     Size of axis labels. Passed on to `barplot`. Default is `1.25`.
#' @param  col_bar     Color of bars. Passed on to `barplot`.
#'                     If `tab` is a 1-D vector, entries in `col_bar` should match
#'                     entries in `tab`. If `tab` is a 2-D matrix, entries in 
#'                     `col_bar` should match the columns of `tab`.
#' @param  col_border  Color of border of bars. Passed on to `barplot`. 
#'                     Default is `NA` (no border).
#' 
#' @return  A barplot with proper placement of x-axis labels. 
#' 
#' @details If `tab` is a 1-D vector, the plot rendered is identical to that 
#'          rendered by `graphics::barplot(tab)` (other than aesthetic details). 
#'          
#'          If `tab` is a 2-D matrix, the plot rendered is identical to that
#'          rendered by `graphics::barplot(t(tab))` (other than aesthetic details).
#'          The rows of `tab` correspond to columns in the barplot. Rownames show
#'          up as x-axis labels.
#'          The columns of `tab` correspond to segments in each column, and 
#'          are colored according to `col_bar`.
#'          
#' @example tab=matrix(c(7.7, 8.1, 7, 72.1, 5.1, 12.6, 9.1, 2, 70.8, 5.5), nrow=2, byrow=T)
#'          colnames(tab) = paste0("Ig", c("A","D","E","G","M"))
#'          rownames(tab) = c("donor_1", "donor_2")
#'          barplot_2(tab, vec_ylim=c(0,105), las=1, col_bar=1:5, ylab="%")
#'          
barplot_2 = function(tab, vec_ylim, vec_labels=NULL, cex.text=1,
                     las=2, ylab="", cex.lab=1.25, 
                     col_bar, col_border=NA) {
    
    if (is.null(dim(tab))) {
        len = length(tab)
    } else {
        len = nrow(tab)
    }
    
    if (is.null(vec_labels)) {
        vec_labels = rep("", len)
    }
    stopifnot(length(vec_labels)==len)
    
    bar_widths = rep(1, len)
    bar_spaces = rep(0.2, len)
    bar_spaces_abs = bar_widths * bar_spaces
    
    bar_se_x = sapply(1:len, function(i){
        if (i==1) {
            # special treatment: sum(widths[1:(i-1)]) won't work properly for i=1
            bar_widths[i]/2 + bar_spaces_abs[i]
        } else {
            # previous bar widths + half current bar width + previous inter-bar spaces 
            sum(bar_widths[1:(i-1)]) + bar_widths[i]/2 + sum(bar_spaces_abs[1:i])
        }
    })
    
    barplot(t(tab), width=bar_widths, space=bar_spaces, 
            ylim=vec_ylim, col=col_bar, border=col_border, 
            las=las, ylab=ylab, cex.lab=cex.lab)
    text(x=bar_se_x, y=vec_ylim[2]*0.98, labels=vec_labels, cex=cex.text)
}

