# adapted from vioplot::vioplot

# medStyle: whether to draw median as point or line
# cexMed: if point, cex of point; if line, lwd of line
# boxWidth: width of box 
# whiskerCol: color of whisker

library(vioplot)

vioplot.julian = function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
          horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
          lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
          at, add = FALSE, wex = 1, drawRect = TRUE,
          cexMed=1, boxWidth=0.3, medStyle=c("point", "line"), whiskerCol="black")  #*
{
    datas <- list(x, ...)
    n <- length(datas)
    if (missing(at)) 
        at <- 1:n
    upper <- vector(mode = "numeric", length = n)
    lower <- vector(mode = "numeric", length = n)
    q1 <- vector(mode = "numeric", length = n)
    q3 <- vector(mode = "numeric", length = n)
    med <- vector(mode = "numeric", length = n)
    base <- vector(mode = "list", length = n)
    height <- vector(mode = "list", length = n)
    baserange <- c(Inf, -Inf)
    args <- list(display = "none")
    if (!(is.null(h))) 
        args <- c(args, h = h)
    for (i in 1:n) {
        data <- datas[[i]]
        data.min <- min(data)
        data.max <- max(data)
        q1[i] <- quantile(data, 0.25)
        q3[i] <- quantile(data, 0.75)
        med[i] <- median(data)
        iqd <- q3[i] - q1[i]
        upper[i] <- min(q3[i] + range * iqd, data.max)
        lower[i] <- max(q1[i] - range * iqd, data.min)
        est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                                   data.max))
        smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                         args))
        hscale <- 0.4/max(smout$estimate) * wex
        base[[i]] <- smout$eval.points
        height[[i]] <- smout$estimate * hscale
        t <- range(base[[i]])
        baserange[1] <- min(baserange[1], t[1])
        baserange[2] <- max(baserange[2], t[2])
    }
    if (!add) {
        xlim <- if (n == 1) 
            at + c(-0.5, 0.5)
        else range(at) + min(diff(at))/2 * c(-1, 1)
        if (is.null(ylim)) {
            ylim <- baserange
        }
    }
    if (is.null(names)) {
        label <- 1:n
    }
    else {
        label <- names
    }
    
    #boxwidth <- 0.05 * wex #*
    boxwidth = boxWidth #*
    
    if (!add) 
        plot.new()
    if (!horizontal) {
        if (!add) {
            plot.window(xlim = xlim, ylim = ylim)
            axis(2)
            axis(1, at = at, label = label)
        }
        box()
        for (i in 1:n) {
            polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
                    c(base[[i]], rev(base[[i]])), col = col, border = border,  
                    lty = lty, lwd = lwd) 
            if (drawRect) {
                lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
                      lty = lty, col=whiskerCol) #*
                rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
                     q3[i], col = rectCol, lwd=lwd) #*
                #* 
                if (medStyle=="point") {
                    points(at[i], med[i], pch = pchMed, col = colMed, cex=cexMed)   
                } else if (medStyle=="line") {
                    segments(x0=at[i] - boxwidth/2, x1=at[i] + boxwidth/2,
                             y0=med[i], y1=med[i], lwd=cexMed, lty=1, col=colMed)
                }
                #* ^
            }
        }
    }
    else {
        if (!add) {
            plot.window(xlim = ylim, ylim = xlim)
            axis(1)
            axis(2, at = at, label = label)
        }
        box()
        for (i in 1:n) {
            polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                                    rev(at[i] + height[[i]])), col = col, border = border, 
                    lty = lty, lwd = lwd)
            if (drawRect) {
                lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
                      lty = lty, col=whiskerCol) #*
                rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
                         boxwidth/2, col = rectCol)
                #* 
                if (medStyle=="point") {
                    points(at[i], med[i], pch = pchMed, col = colMed, cex=cexMed)   
                } else if (medStyle=="line") {
                    segments(x0=at[i] - boxwidth/2, x1=at[i] + boxwidth/2,
                             y0=med[i], y1=med[i], lwd=cexMed, lty=1, col=colMed)
                }
                #* ^
            }
        }
    }
    invisible(list(upper = upper, lower = lower, median = med, 
                   q1 = q1, q3 = q3))
}

# vecLst: a list of vectors (names of list entries go into x axis)
# yHoriz: numeric value for horizontal dashed line
# xLas: 2 = vertical, 1=horizontal
# yLas: 3 = vertical, 2=horizontal

# useBoxplot: use boxplot instead of violin plot
# showPoints: show raw points (not recommended to use with useBoxplot=T if there're lots of points)
# configs for showPoints=T:
# - pointsPch: VECTORS matching length of vecLst
# - pointsCol: VECTORS matching length of vecLst (recommen using transparency via scales::alpha)
# - pointsCex: single value, cex of points
# - pointsJitterAmount: single value, passed to `amount` in `jitter``
# - pointsJitterBool: if TRUE, jitter; if FALSE, do not jitter
vioplotLst = function(vecLst, vioYmax=NULL, vioYmin=0, xLas=2, yLas=3, 
                      xNames=NULL, xTick=TRUE, xlab="", ylab="", title="", 
                      countCex=0.7, xaxisCex=1, yaxisCex=1, mainCex=1.2,
                      yHoriz=NULL, yHorizCol="firebrick2",
                      vioBkgCol=NA, vioMedCol="hotpink", vioRecCol="darkgoldenrod2",
                      vioMedPch=18, vioMedStyle="point", vioBoxWidth=0.2, vioWhiskerCol="black",
                      useBoxplot=F, showPoints=F, 
                      pointsPch, pointsCol, pointsCex=0.3, 
                      pointsJitterAmount=0.4, pointsJitterBool=TRUE) {
    nVec = length(vecLst)
    
    if (is.null(vioYmax)) {
        vioYmax = max(unlist(vecLst), na.rm=T)*1.05
    }
    
    # base plot
    plot(x=0:(nVec+1), y=rep(0, nVec+2),
         col=NA, ylim=c(vioYmin, vioYmax), xlim=c(0.5, nVec+0.5), 
         xaxt="n", yaxt="n", xlab=xlab, ylab=ylab, main=title, 
         cex.lab=1.2, cex.axis=1.2, cex.main=mainCex)
    axis(side=1, at=1:nVec, labels=xNames, las=xLas, cex.axis=xaxisCex, tick=xTick)
    axis(side=2, las=yLas, cex.axis=yaxisCex)
    
    # individual plots
    for (i in 1:nVec) {
        
        if (!useBoxplot) {
            
            # do not attemp,t to plot if NULL or if not NULL but of length 0
            if (!is.null(vecLst[[i]]) && length(vecLst[[i]])>0) {
                # vioplot may break if all values are the same
                # e.g. vioplot(rep(0, 113))
                if ( length(unique(vecLst[[i]])) > 1) {
                    
                    # if adding raw points
                    # do this BEFORE drawing vioplot
                    if (showPoints) {
                        if (pointsJitterBool) {
                            pts_x = jitter(x=rep(i, length=length(vecLst[[i]])), amount=pointsJitterAmount)
                        } else {
                            pts_x = rep(i, length=length(vecLst[[i]]))
                        }
                        points(x=pts_x, y=vecLst[[i]], cex=pointsCex, col=pointsCol[i], pch=pointsPch[i])
                    }
                    
                    if (is.null(xNames)) {
                        xName_i = ""
                    } else {
                        xName_i = xNames[i]
                    }
                    
                    vioplot.julian(vecLst[[i]], names=xName_i, add=T, at=i,
                                   col=vioBkgCol, rectCol=vioRecCol, whiskerCol=vioWhiskerCol,
                                   colMed=vioMedCol, pchMed=vioMedPch, medStyle=vioMedStyle, cexMed=1.5, 
                                   drawRect=T, boxWidth=vioBoxWidth)
                    #vioplot(vecLst[[i]], add=T, at=i,
                    #        col=vioBkgCol, rectCol=vioRecCol,
                    #        colMed=vioMedCol, pchMed=vioMedPch,
                    #        drawRect=T)                    
                    
                } else {
                    points(x=i, y=unique(vecLst[[i]]), col=vioMedCol, pch=vioMedPch, cex=1.5)
                }
                text(x=i, y=vioYmax, labels=length(vecLst[[i]]), cex=countCex)
            } else {
                text(x=i, y=vioYmax, labels="0", cex=countCex)
            }
        } else {
            # boxplot
            # do not attemp,t to plot if NULL or if not NULL but of length 0
            if (!is.null(vecLst[[i]]) && length(vecLst[[i]])>0) {

                if (is.null(xNames)) {
                    xName_i = ""
                } else {
                    xName_i = xNames[i]
                }
            
                # if adding raw points
                # do this AFTER drawing vioplot
                # in general, doesn't turn out well if there're too many points even if drawing 
                # afterwards
                if (showPoints) {
                    if (pointsJitterBool) {
                        pts_x = jitter(x=rep(i, length=length(vecLst[[i]])), amount=pointsJitterAmount)
                    } else {
                        pts_x = rep(i, length=length(vecLst[[i]]))
                    }
                    points(x=pts_x, y=vecLst[[i]], cex=pointsCex, col=pointsCol[i], pch=pointsPch[i])
                }
                
                # add boxplot after points so that it doesn't get obscured by points
                # do not plot outliers if showPoints is T
                # otherwise will end up double-plotting
                boxplot(vecLst[[i]], add=T, at=i, col=vioBkgCol, yaxt="n",
                        names=xName_i, outline=!showPoints)  
                
                text(x=i, y=vioYmax, labels=length(vecLst[[i]]), cex=countCex)
            } else {
                text(x=i, y=vioYmax, labels="0", cex=countCex)
            }
        }

    }
    
    # threshold
    if (!is.null(yHoriz)) {
        abline(h=yHoriz, col=2, lty=2, lwd=1.5)
    }
}
