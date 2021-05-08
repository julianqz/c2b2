# overview of plot
# upper section:
# - histogram of within-group DTN 
# - density curve (optional)
# lower section:
# - histogram of cross-group DTN 

# plotFreq: if TRUE, plot frequency (counts); if FALSE, plot density 
# vec1/2, line1_x/y: numeric vector
# line1_x/y is optional
# vecMax: if NULL, will be derived from vec1 and vec2
# transparent: [0, 1]; when 1, color is opaque
# vec1/2TxtPos, vec1/2TxtHeight: (0, 1]
plotCrossHam = function(vec1, vec2, plotFreq=T, binSize=0.01,
                        vecMax=NULL,
                        line1_x=NULL, line1_y=NULL,
                        vec1Col="royalblue", vec2Col="springgreen", transparency=1,
                        plotTitle="", vec1Txt="", vec2Txt="", 
                        vec1TxtPos=0.7, vec2TxtPos=0.7, 
                        vec1TxtHeight=0.8, vec2TxtHeight=0.8, 
                        labX="Hamming Distance (Normalized by Length)",
                        labY="",
                        zeroLine=F, zeroLineLty=2, zeroLineCol="darkgoldenrod2", 
                        cexLab=1, cexAxis=1, cexMain=1, cexTxt=1) {
    
    # if line1_x/y is supplied, only supports plotFreq=F
    if (!is.null(line1_x)) {
        stopifnot( !is.null(line1_y) )
        stopifnot( !plotFreq )
    }
    
    # check if all NA
    vec1NAbool = all(is.na(vec1))
    vec2NAbool = all(is.na(vec2))
    
    if (vec1NAbool & vec2NAbool) {
        message("vec1 and vec2 contain only NAs. Nothing is plotted.")
    } else {
        
        # global max
        # expect line1_x range to comform to vec1 
        # (which is the case if line1 is density based on vec1)
        if (is.null(vecMax)) {
            vecMax = max(c(vec1, vec2), na.rm=T)
        } 
        
        # vec1
        # assumes no negative in vec1
        vec1Hist = hist(vec1, breaks=seq(from=0, to=vecMax+0.05, by=binSize), plot=F)
        
        # vec2 
        # assumes no negative in vec2
        vec2Hist = hist(vec2, breaks=seq(from=0, to=vecMax+0.05, by=binSize), plot=F)
        
        # set density to negative (for plotting purpose)
        vec2Hist$density = -vec2Hist$density
        vec2Hist$counts = -vec2Hist$counts
        
        if (plotFreq) {
            vec1HistMax = max(vec1Hist$counts)
            vec2HistMax = min(vec2Hist$counts) # min because it's been negated
        } else {
            #vec1HistMax = max(vec1Hist$density)
            vec1HistMax = max(c(vec1Hist$density, line1_y)) # ok if line1_y is NULL
            vec2HistMax = min(vec2Hist$density) # min because it's been negated
        }
        
        
        ### plot
        # ?plot.histogram
        
        # vec1
        plot(vec1Hist, freq=plotFreq, 
             col=scales::alpha(vec1Col, transparency), border=NA, # if border not NA, tiny bit shows up for density=0
             xlim=c(0, vecMax*1.05), #* 
             ylim=c(vec2HistMax, vec1HistMax), #*
             xlab=labX, ylab=labY, yaxt="n",
             main=plotTitle, 
             cex.lab=cexLab, cex.axis=cexAxis, cex.main=cexMain)
        
        # line1
        if ( !is.null(line1_x) & !is.null(line1_y) ) {
            points(x=line1_x, y=line1_y, type="l", col="purple", lty=2, lwd=2)
        }
        
        # y-axis ticks
        yticks.by.base = c(2, 5, 10, 20, 50, 100)
        yticks.by.plot = yticks.by.base[which.min(abs( (vec1HistMax-vec2HistMax) / 10 - yticks.by.base ))]
        
        if ( (!vec1NAbool) & (!vec2NAbool) ) {
            # plot y-axis for both groups
            yticks.plot = sort(union(seq(from=0, to=vec1HistMax, by=yticks.by.plot),
                                     seq(from=0, to=vec2HistMax, by=-yticks.by.plot)))
        } else if ( vec1NAbool & (!vec2NAbool) ) {
            # only plot y-axis for group 2
            yticks.plot = sort(seq(from=0, to=vec2HistMax, by=-yticks.by.plot))
        } else if ( (!vec1NAbool) & vec2NAbool ) {
            # only plot y-axis for group 1
            yticks.plot = sort(seq(from=0, to=vec1HistMax, by=yticks.by.plot))
        }
        
        axis(side=2, at=yticks.plot, labels=as.character(abs(yticks.plot)),
             cex.axis=cexAxis, cex.lab=cexLab)
        
        # vec2
        plot(vec2Hist, freq=plotFreq, add=T, col=scales::alpha(vec2Col, transparency), border=NA)
        
        # horizontal line at count 0
        if (zeroLine) {
            abline(h=0, lwd=1.5, lty=zeroLineLty, col=zeroLineCol)
        }
        
        # text label
        text(x=vec1TxtPos*(vecMax+0.05), y=vec1TxtHeight*vec1HistMax, pos=4, labels=vec1Txt, cex=cexTxt)
        text(x=vec2TxtPos*(vecMax+0.05), y=vec2TxtHeight*vec2HistMax, pos=4, labels=vec2Txt, cex=cexTxt)
    
    }
    
}
