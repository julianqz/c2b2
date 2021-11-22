#### percentage barplot ####

# TODO: add doc

barplot_2 = function(tab, vec_ylim, vec_labels=NULL, 
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
    text(x=bar_se_x, y=vec_ylim[2]*0.98, labels=vec_labels)
}
