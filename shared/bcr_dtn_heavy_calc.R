#!/opt/conda/bin/Rscript

# wrapper to calculate dist-to-nearest using heavy chains only

# assumes:
# - pathCSV points to a comma-separated file with the following headers
#   "subj", "path_db" where "path_db" points to a .tsv file

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--pathCSV", action="store", default=NA, 
                type="character", help="Path to CSV containing subject list and paths to input files."),
    make_option("--pathWork", action="store", default=NA, type="character", 
                help="path_work."),
    make_option("--nproc", action="store", default=1, type="numeric", 
                help="nproc."),
    make_option("--calcWithin", action="store", default=FALSE, type="logical", 
                help="Whether to calculate within-subject dtn."),
    make_option("--calcBetween", action="store", default=FALSE, type="logical", 
                help="Whether to calculate bewteen-subject dtn."),
    make_option("--subsampleWithin", action="store", default=NULL, type="numeric", 
                help="Within-subject subsampling. Do not specify via command line if NULL."),
    make_option("--subsampleBetween", action="store", default=NULL, type="numeric", 
                help="Between-subject subsampling. Do not specify via command line if NULL."),
    make_option("--colSubj", action="store", default=NA, 
                type="character", help="Column name containing subject info."),
    make_option("--colSeqID", action="store", default="sequence_id", 
                type="character", help="Column name containing sequence ID."),
    make_option("--colSeq", action="store", default="junction", 
                type="character", help="col_seq."),
    make_option("--colV", action="store", default="v_call", 
                type="character", help="col_v."),
    make_option("--colJ", action="store", default="j_call", 
                type="character", help="col_j.")
)
opt = parse_args(OptionParser(option_list=option_list))

subj_info = read.table(opt$pathCSV, header=T, sep=",", stringsAsFactors=F)

suppressPackageStartupMessages(library(shazam))

setwd(opt$pathWork)
sinkName = paste0("computingEnv_dtn_heavy_", Sys.Date(), "-", 
                  format(Sys.time(), "%H%M%S"), '.txt')
sink(sinkName)
cat("calcWithin:", opt$calcWithin, "\n")
# if NULL, will appear as "subsampleWithin: " (i.e. blank)
cat("subsampleWithin:", opt$subsampleWithin, "\n")
cat("calcBetween:", opt$calcBetween, "\n")
# if NULL, will appear as "subsampleBetween: " (i.e. blank)
cat("subsampleBetween:", opt$subsampleBetween, "\n")
sessionInfo()
sink()


#### within-subject ####

if (opt$calcWithin) {
    
    cat("\nCalculating within-subject distToNearest... \n")
    
    sink_name = "thresh_density_all.txt"
    
    for (i in 1:nrow(subj_info)) {
        
        db = read.table(subj_info[["path_db"]][i],
                        header=T, sep="\t", stringsAsFactors=F)
        
        subj = subj_info[["subj"]][i]
        
        cat("\n", subj, "; nrow(db):", nrow(db), "\n")
        
        fn = paste0("dtn_", subj,
                    ifelse(is.null(opt$subsampleWithin), "", 
                           paste0("_subsample-", opt$subsampleWithin)),
                    ".RData")
        
        # adds $dist_nearest and $vjl_group columns
        db = distToNearest(db, 
                           sequenceColumn=opt$colSeq,
                           vCallColumn=opt$colV,
                           jCallColumn=opt$colJ,
                           subsample=opt$subsampleWithin,
                           model="ham", normalize="len", 
                           first=F, VJthenLen=F, keepVJLgroup=T,
                           nproc=opt$nproc, progress=T)
        
        save(db, file=fn)
        
        # density estimate for threshold
        # if already subsampled in previous step, that will carry over 
        # (no need to re-subsample)
        thresh_obj = findThreshold(distances=db[["dist_nearest"]],
                                   method="density", progress=T)
        
        fn = paste0("thresh_density_", subj, ".RData")
        save(thresh_obj, file=fn)
        
        # print threshold
        if (i==1) {
            # initiate
            sink(file=sink_name, append=F)
            cat(subj, sep="", "\n")
            cat(thresh_obj@threshold, sep="", "\n")
            sink()
        } else {
            # append
            sink(file=sink_name, append=T)
            cat(subj, sep="", "\n")
            cat(thresh_obj@threshold, sep="", "\n")
            sink()
        }
        
        rm(db, thresh_obj)
    }
}


#### between-subject ####

if (opt$calcBetween) {
    
    # concat all subjects
    
    cols_keep = c(opt$colSeqID, opt$colSeq, opt$colV, opt$colJ)
    
    for (i in 1:nrow(subj_info)) {
        
        db_tmp = read.table(subj_info[["path_db"]][i],
                            header=T, sep="\t", stringsAsFactors=F)
        
        stopifnot(all(cols_keep %in% colnames(db_tmp)))
        db_tmp = db_tmp[, cols_keep]
        db_tmp[[opt$colSubj]] = subj_info[["subj"]][i]
        
        if (i==1) {
            # initiate db
            db = db_tmp
        } else {
            # append new rows to db
            db = rbind(db, db_tmp)
        }
        
        rm(db_tmp)
    }
    
    cat("\nBreakdown by subject:\n")
    print(table(db[[opt$colSubj]]))
    
    cat("\nCalculating between-subject distToNearest... \n")

    fn = paste0("dtn_btwSubj",
                ifelse(is.null(opt$subsampleBetween), "", 
                       paste0("_subsample-", opt$subsampleBetween)),
                ".RData")
    
    #* v1.0.2 stable release requires a bug fix from commit 47bffe0
    #* this bug fix is packed into julianqz/wu_cimm:main_0.1.1
    #* but if not using that docker image, v1.0.2 will fail next line
    db = distToNearest(db, 
                       sequenceColumn=opt$colSeq,
                       vCallColumn=opt$colV,
                       jCallColumn=opt$colJ,
                       cross=opt$colSubj,
                       subsample=opt$subsampleBetween,
                       model="ham", normalize="len", 
                       first=F, VJthenLen=F, keepVJLgroup=T,
                       nproc=opt$nproc, progress=T)
    
    save(db, file=fn)
    
    rm(db, fn)
}

