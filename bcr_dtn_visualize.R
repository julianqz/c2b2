#!/opt/conda/bin/Rscript

# Julian Q. Zhou
# https://github.com/julianqz

# wrapper to visualize dist-to-nearest

# assumes:
# - findThreshold was called together with distToNearest
# - pathCSV points to a comma-separated file with the following headers
#   "subj", "path_dtn" where "path_dtn" contains a .RData file from findThreshold

# If there is only 1 subject in --pathCSV, even if --plotBetween is TRUE, 
# between-subject visualization will be skipped 

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--pathCSV", action="store", default=NA, 
                type="character", help="Path to CSV containing subject list and paths to input files."),
    make_option("--pathWork", action="store", default=NA, type="character", 
                help="path_work."),
    make_option("--pathBtw", action="store", default=NA, type="character", 
                help="path_btw. Path to .RData file containing $cross_dist_nearest."),
    make_option("--pathHelper", action="store", default=NA, type="character", 
                help="path_helper. Path to helper script containing plotCrossHam()."),
    make_option("--plotWithin", action="store", default=FALSE, type="logical", 
                help="Whether to visualize within-subject dtn."),
    make_option("--plotBetween", action="store", default=FALSE, type="logical", 
                help="Whether to visualize bewteen-subject dtn."),
    make_option("--colSubj", action="store", default=NA, 
                type="character", help="Column name containing subject info."),
    make_option("--heavyLight", action="store", default=FALSE, type="logical", 
                help="Whether partition was based on both heavy and light chains."),
    make_option("--nPlotsPerRow", action="store", default=1, 
                type="numeric", help="n_plots_per_row."),
    make_option("--manualThresh", action="store", type="character", 
                help="thresh_man_vec. Length must match num of subjects. Can be all NA. Eg: 'NA,NA,NA'.")
)
opt = parse_args(OptionParser(option_list=option_list))


#### global config ####

suppressPackageStartupMessages(library(shazam))
# plotCrossHam
source(opt$pathHelper)

# suffix for filenames
if (opt$heavyLight) {
    out_suffix = "_groupByHL"
} else {
    out_suffix = "_groupByHonly"
}

subj_info = read.table(opt$pathCSV, header=T, sep=",", stringsAsFactors=F)
n_subj = nrow(subj_info)
subjects = subj_info[["subj"]]

path_work = opt$pathWork
path_btw = opt$pathBtw

plot_within = opt$plotWithin
plot_btw = opt$plotBetween

col_subj = opt$colSubj

n_plots_per_row = opt$nPlotsPerRow
n_rows = ceiling(n_subj / n_plots_per_row)

# parse
# \s is space
# ? means preceding item is optional and will be matched at most once

# first parse into strings
thresh_man_str = strsplit(opt$manualThresh, "\\s?,\\s?")[[1]]
# check length
stopifnot( length(thresh_man_str) == n_subj )
# convert non-NA into integer
thresh_man_vec = rep(NA, length=n_subj)
for (i in 1:n_subj) {
    if (thresh_man_str[i]!="NA") {
        thresh_man_vec[i] = as.numeric(thresh_man_str[i])
    }
}
names(thresh_man_vec) = subjects
# for debugging
# thresh_man_vec: NA NA NA ; class: logical 
# thresh_man_vec: NA 0.23 NA ; class: numeric 
#cat("\nthresh_man_vec:", thresh_man_vec, "; class:", class(thresh_man_vec), "\n") 


# not provided as command line options for now
cex_lab_within = 1.25
cex_axis_within = 1.25
cex_main_within = 1.25

cex_lab_btw = 2
cex_axis_btw = 2
cex_main_btw = 2

lab_x = "Distance to Nearest Neighbor"

#col_dist_within = "dist_nearest"
col_dist_btw = "cross_dist_nearest"


#### common ####

if (plot_within | plot_btw) {
    
    cat("\nLoading... \n")
    
    # loop thru every subject, $dist_nearest, store in a list
    # distances from all subjects are needed to get standardized xlim
    # do this before plotting to avoid having to load each db a second time
    
    # instead of loading db, can load thresh_obj which also contains $dist_nearest
    dist_lst = vector(mode="list", length=n_subj)
    xden_lst = vector(mode="list", length=n_subj)
    yden_lst = vector(mode="list", length=n_subj)
    thresh_auto_vec = rep(NA, length=n_subj)
    names(dist_lst) = subjects
    names(xden_lst) = subjects
    names(yden_lst) = subjects
    names(thresh_auto_vec) = subjects
    
    for (i in 1:n_subj) {
        
        # DensityThreshold object:
        # - x: input distance vector with NA or infinite values removed
        # - bandwith: for density estimation
        # - xdens: vector of distance value for smoothed density estimate
        # - ydens: vector of density value for smoothed density estimate
        # - threshold: distance threshold
        
        subj = subjects[i]
        cat(subj, "\n")
        
        # thresh_obj
        fn = paste0(subj_info[["path_dtn"]][i],
                    "thresh_density", out_suffix, "_", subj, ".RData")
        load(fn)
        
        if (!is.null(thresh_obj)) {
            dist_lst[[subj]] = thresh_obj@x
            xden_lst[[subj]] = thresh_obj@xdens
            yden_lst[[subj]] = thresh_obj@ydens
            thresh_auto_vec[subj] = thresh_obj@threshold
            
        } else {
            # this should set `all(is.na(cur_dist))` to be FALSE,
            # triggering a blank plot being produced
            dist_lst[[subj]] = NA
        }

        rm(thresh_obj)
    }
}


#### within-subject ####

if (plot_within) {
    
    cat("\nVisualizing within-subject distToNearest... \n")
    
    # standardized xlim for all subjects
    global_dist_range = range(unlist(dist_lst), na.rm=T)
    
    # 1: histogram
    # 2: histogram + density curve
    # 3: histogram + density curve + auto thresh
    # 4: histogram + density curve + manual thresh
    
    setwd(path_work)
    
    for (config in 1:4) {
        
        cat("config", config, "\n")
        
        fn = paste0("dtn_within", out_suffix, "_config-", config, ".pdf")
        pdf(fn, width=4*n_plots_per_row, height=4*n_rows)
        par(mfrow=c(n_rows, n_plots_per_row), mar=c(5,5,2,2)) # BLTR
        
        for (subj in subjects) {
            cat(subj, "\n")
            
            cur_dist = dist_lst[[subj]]
            cur_xden = xden_lst[[subj]]
            cur_yden = yden_lst[[subj]]
            cur_thresh_auto = thresh_auto_vec[subj]
            cur_thresh_man = thresh_man_vec[subj]
            
            if (all(is.na(cur_dist))) {
                
                # blank plot if all NA
                plot(0,0,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
                title(paste0(subj, " (all NAs)"))
                
            } else {
            
                ### histogram 
                
                # don't plot just yet
                cur_hist = hist(x=cur_dist, 
                               breaks=seq(from=0, to=global_dist_range[2]+0.03, by=0.02),
                               plot=F)
                
                local_y_max = max(c(cur_yden,
                                    cur_hist$density[!is.infinite(cur_hist$density)]), 
                                  na.rm=T)
                
                # plot
                plot(cur_hist, freq=F, xlim=c(0, global_dist_range[2]+0.03),
                     col=scales::alpha("skyblue", 0.7), border="white",
                     xlab=lab_x, main=subj,
                     cex.axis=cex_axis_within, cex.lab=cex_lab_within, 
                     cex.main=cex_main_within)
                
                ### density curve
                
                if (config %in% 2:4) {
                    points(x=cur_xden, y=cur_yden, type="l", col="purple", lty=2, lwd=1.5)
                }
                
                ### auto threshold
                
                if (config==3) {
                    if (!is.na(cur_thresh_auto)) {
                        abline(v=cur_thresh_auto, lty=2, lwd=1.5, col="firebrick2")
                        text(x=cur_thresh_auto, y=local_y_max*0.7, col="firebrick2", pos=4,
                             labels=as.character(round(cur_thresh_auto, 3)))
                    } else {
                        cat("auto thresh is NA; skipped plotting abline\n")
                    }
                }
                
                ### manual threshold
                
                if (config==4) {
                    if (!is.na(cur_thresh_man)) {
                        abline(v=cur_thresh_man, lty=2, lwd=1.5, col="black")
                        text(x=cur_thresh_man, y=local_y_max*0.9, col="black", pos=4,
                             labels=as.character(round(cur_thresh_man, 3)))
                    } else {
                        cat("manual thresh is NA; skipped plotting abline\n")
                    }
                }

            }
        }
        
        dev.off()
    }

}


#### between-subject ####

if (plot_btw) {
    
    if (n_subj>1) {
        
        cat("\nVisualizing within- & between-subject distToNearest... \n")
        
        # btw-subj distances
        # db, with $cross_dist_nearest
        load(path_btw)
        
        # standardized xlim for all subjects
        global_dist_range = range(c(unlist(dist_lst), db[[col_dist_btw]]), na.rm=T)
        
        # 1: histogram
        # 2: histogram + density curve
        # 3: histogram + density curve + auto thresh
        # 4: histogram + density curve + manual thresh
        
        setwd(path_work)
        
        for (config in 1:4) {
            
            cat("config", config, "\n")
            
            fn = paste0("dtn_btw", out_suffix, "_config-", config, ".pdf")
            pdf(fn, width=5*n_plots_per_row, height=5*n_rows)
            par(mfrow=c(n_rows, n_plots_per_row), mar=c(5,5,2,2)) # BLTR
            
            for (subj in subjects) {
                cat(subj, "\n")
                
                cur_dist = dist_lst[[subj]]
                cur_xden = xden_lst[[subj]]
                cur_yden = yden_lst[[subj]]
                cur_thresh_auto = thresh_auto_vec[subj]
                cur_thresh_man = thresh_man_vec[subj]
                cur_btw = db[[col_dist_btw]][ db[[col_subj]]==subj ]
                
                if (all(is.na(cur_dist))) {
                    
                    # blank plot if all NA
                    plot(0,0,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
                    title(paste0(subj, " (all NAs)"))
                    
                } else {
                    
                    # local ylims will be calculated by plotCrossHam
                    # but need to re-calculate in case threshold is to be noted
                    local_y_max = max(c( hist(cur_dist, breaks=seq(from=0, to=global_dist_range[2]+0.03, by=0.02), plot=F)$density,
                                         hist(cur_btw, breaks=seq(from=0, to=global_dist_range[2]+0.03, by=0.02), plot=F)$density,
                                         cur_yden ), na.rm=T)
                    
                    ### histogram and density curve
                    
                    if (config==1) {
                        cur_xden_plot = NULL
                        cur_yden_plot = NULL
                    } else {
                        cur_xden_plot = cur_xden
                        cur_yden_plot = cur_yden
                    }
                    
                    plotCrossHam(vec1=cur_dist, vec2=cur_btw, vecMax=global_dist_range[2],
                                 line1_x=cur_xden_plot, line1_y=cur_yden_plot,
                                 plotFreq=F, binSize=0.02,
                                 plotTitle=subj, labX=lab_x, labY="Density",
                                 vec1Col="skyblue", vec2Col="springgreen", transparency=1,
                                 zeroLine=F, zeroLineLty=2, zeroLineCol="darkgoldenrod2", 
                                 cexLab=cex_lab_btw, cexAxis=cex_axis_btw, 
                                 cexMain=cex_main_btw, cexTxt=1)
                    
                    ### auto threshold
                    
                    if (config==3) {
                        if (!is.na(cur_thresh_auto)) {
                            abline(v=cur_thresh_auto, lty=2, lwd=1.5, col="firebrick2")
                            text(x=cur_thresh_auto, y=local_y_max*0.7, col="firebrick2", pos=4,
                                 labels=as.character(round(cur_thresh_auto, 3)))
                        } else {
                            cat("auto thresh is NA; skipped plotting abline\n")
                        }
                    }
                    
                    ### manual threshold
                    
                    if (config==4) {
                        if (!is.na(cur_thresh_man)) {
                            abline(v=cur_thresh_man, lty=2, lwd=1.5, col="black")
                            text(x=cur_thresh_man, y=local_y_max*0.9, col="black", pos=4,
                                 labels=as.character(round(cur_thresh_man, 3)))
                        } else {
                            cat("manual thresh is NA; skipped plotting abline\n")
                        }
                    }
                    
                }
            }
            
            dev.off()
        }    
        
        
        rm(db)
        
    } else {
        cat("\nOnly 1 subject found. Skipping between-subject.\n")
    }
    
}

