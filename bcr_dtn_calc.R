#!/opt/conda/bin/Rscript

# Julian Q. Zhou
# https://github.com/julianqz

# wrapper to calculate dist-to-nearest for heavy chains after partitioning
# - either based on heavy chains only
# - or based on both heavy and light chains

# assumes:
# - pathCSV points to a comma-separated file with the following headers
#   "subj", "path_db_heavy/light" where "path_db_*" points to a .tsv file
# - note that heavy and light chains are expected to be in separate files

# If the data is single-cell heavy:light paired, but partitioning based on
# heavy chains only is desired, simply set `--heavyLight` to `FALSE` and 
# light chains will be ignored.

# If there's only 1 subject in --pathCSV, even if --calcBetween is TRUE,
# between-subject calculation will be skipped.

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
    make_option("--calcThreshold", action="store", default=TRUE, type="logical", 
                help="Whether to run findThreshold."),
    make_option("--colSubj", action="store", default=NA, 
                type="character", help="Column name containing subject info."),
    make_option("--colSeqID", action="store", default="sequence_id", 
                type="character", help="Column name containing sequence ID."),
    make_option("--colSeq", action="store", default="cdr3", 
                type="character", help="sequenceColumn."),
    make_option("--colV", action="store", default="v_call", 
                type="character", help="vCallColumn."),
    make_option("--colJ", action="store", default="j_call", 
                type="character", help="jCallColumn."),
    make_option("--heavyLight", action="store", default=FALSE, type="logical", 
                help="Whether to partition using both heavy and light chains."),
    make_option("--colCell", action="store", default="cell_id", type="character", 
                help="cellIdColumn. Ignored if --heavyLight FALSE."),
    make_option("--colLocus", action="store", default="locus", type="character", 
                help="locusColumn. Ignored if --heavyLight FALSE.")
)
opt = parse_args(OptionParser(option_list=option_list))

subj_info = read.table(opt$pathCSV, header=T, sep=",", stringsAsFactors=F)

# check columns
# suffix for output filenames
if (opt$heavyLight) {
    stopifnot( all(c("path_db_heavy", "path_db_light") %in% colnames(subj_info)) )
    out_suffix = "_groupByHL"
} else {
    stopifnot( "path_db_heavy" %in% colnames(subj_info) )
    out_suffix = "_groupByHonly"
}

suppressPackageStartupMessages(library(shazam))

setwd(opt$pathWork)
sinkName = paste0("computingEnv_dtn_heavy_", Sys.Date(), "-", 
                  format(Sys.time(), "%H%M%S"), '.txt')
sink(sinkName)
cat("Partition using both heavy and light:", opt$heavyLight, "\n")
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
    
    sink_name = paste0("thresh_density", out_suffix, "_all.txt")
    
    for (i in 1:nrow(subj_info)) {
        
        subj = subj_info[["subj"]][i]
        
        if (opt$heavyLight) {
            # heavy and light
            
            # .tsv
            #db_h = read.table(subj_info[["path_db_heavy"]][i],
            #                  header=T, sep="\t", stringsAsFactors=F)
            #db_l = read.table(subj_info[["path_db_light"]][i],
            #                  header=T, sep="\t", stringsAsFactors=F)
            
            # .RData
            load(subj_info[["path_db_heavy"]][i])
            db_h = db; rm(db)
            load(subj_info[["path_db_light"]][i])
            db_l = db; rm(db)
            
            # columns should match
            stopifnot(all.equal(colnames(db_h), colnames(db_l)))
            
            stopifnot(all( c(opt$colCell, opt$colLocus) %in% colnames(db_h) ))
            
            # each cell should have 1 HC and 1 LC each
            cells_common = base::intersect(db_h[[opt$colCell]],
                                           db_l[[opt$colCell]])
            bool_common_h = db_h[[opt$colCell]] %in% cells_common
            bool_common_l = db_l[[opt$colCell]] %in% cells_common
            if (any(!bool_common_h)) {
                cat("\n", subj, "- excluded", sum(!bool_common_h), 
                    "heavy chain seqs for lacking light chain counterparts\n")
                db_h = db_h[bool_common_h, ]
            }
            if (any(!bool_common_l)) {
                cat("\n", subj, "- excluded", sum(!bool_common_l), 
                    "light chain seqs for lacking heavy chain counterparts\n")
                db_l = db_l[bool_common_l, ]
            }
            stopifnot(nrow(db_h)==nrow(db_l))
            
            db = rbind(db_h, db_l)
            
        } else {
            # heavy only
            
            # .tsv
            #db = read.table(subj_info[["path_db_heavy"]][i],
            #                header=T, sep="\t", stringsAsFactors=F)
            
            # .RData
            load(subj_info[["path_db_heavy"]][i])
            
        }
        
        nrow_bf = nrow(db)
        cat("\n", subj, "; nrow(db):", nrow_bf, "\n")
        
        fn = paste0("dtn", out_suffix, "_", subj, 
                    ifelse(is.null(opt$subsampleWithin), "", 
                           paste0("_subsample-", opt$subsampleWithin)),
                    ".RData")
        
        # adds $dist_nearest and $vjl_group columns
        
        if (opt$heavyLight) {
            # heavy and light
            db = distToNearest(db, 
                               sequenceColumn=opt$colSeq,
                               vCallColumn=opt$colV,
                               jCallColumn=opt$colJ,
                               cellIdColumn=opt$colCell,
                               locusColumn=opt$colLocus,
                               onlyHeavy=F,
                               subsample=opt$subsampleWithin,
                               model="ham", normalize="len", 
                               first=F, VJthenLen=F, keepVJLgroup=T,
                               nproc=opt$nproc, progress=T)
            
            # sanity check
            # heavy and light chains from the same cell should have the same partition
            stopifnot(all( sapply(unique(db[[opt$colCell]]), 
                                  function(s){
                                      # wrt db
                                      s_idx = which(db[[opt$colCell]]==s)
                                      s_bool = length(unique(db[["vjl_group"]][s_idx]))==1
                                      return(s_bool)
                                  }, USE.NAMES=F) ))
            # all counts of vjl_group should be even (heavy:light paired)
            stopifnot(all( table(db[["vjl_group"]]) %% 2 == 0 ))
            
        } else {
            # heavy only
            db = distToNearest(db, 
                               sequenceColumn=opt$colSeq,
                               vCallColumn=opt$colV,
                               jCallColumn=opt$colJ,
                               subsample=opt$subsampleWithin,
                               model="ham", normalize="len", 
                               first=F, VJthenLen=F, keepVJLgroup=T,
                               nproc=opt$nproc, progress=T)
        }
        
        nrow_af = nrow(db)
        # sanity check
        # number of rows should remain the same
        stopifnot( nrow_bf == nrow_af )
        
        save(db, file=fn)
        
        if (opt$calcThreshold) {
            
            # density estimate for threshold
            # if already subsampled in previous step, that will carry over 
            # (no need to re-subsample)
            thresh_obj = findThreshold(distances=db[["dist_nearest"]],
                                       method="density", progress=T)
            
            fn = paste0("thresh_density", out_suffix, "_", subj, ".RData")
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
            
            rm(thresh_obj)
        }
        
        rm(db)
        
    }
}


#### between-subject ####

if (opt$calcBetween) {
    
    # skip if only 1 donor
    if (nrow(subj_info)>1) {
        
        # concat all subjects
        if (opt$heavyLight) {
            cols_keep = c(opt$colSeqID, opt$colSeq, opt$colV, opt$colJ,
                          opt$colCell, opt$colLocus)
        } else {
            cols_keep = c(opt$colSeqID, opt$colSeq, opt$colV, opt$colJ)
        }
        
        for (i in 1:nrow(subj_info)) {
            
            subj = subj_info[["subj"]][i]
            
            if (opt$heavyLight) {
                # heavy and light
                
                # .tsv
                #db_tmp_h = read.table(subj_info[["path_db_heavy"]][i],
                #                      header=T, sep="\t", stringsAsFactors=F)
                #db_tmp_l = read.table(subj_info[["path_db_light"]][i],
                #                      header=T, sep="\t", stringsAsFactors=F)
                
                # .RData
                load(subj_info[["path_db_heavy"]][i])
                db_tmp_h = db; rm(db)
                load(subj_info[["path_db_light"]][i])
                db_tmp_l = db; rm(db)
                
                stopifnot(all.equal(colnames(db_tmp_h), colnames(db_tmp_l)))
                
                stopifnot(all( c(opt$colCell, opt$colLocus) %in% colnames(db_tmp_h) ))
                
                # each cell should have 1 HC and 1 LC each
                cells_common = base::intersect(db_tmp_h[[opt$colCell]],
                                               db_tmp_l[[opt$colCell]])
                bool_common_h = db_tmp_h[[opt$colCell]] %in% cells_common
                bool_common_l = db_tmp_l[[opt$colCell]] %in% cells_common
                if (any(!bool_common_h)) {
                    cat("\n", subj, "- excluded", sum(!bool_common_h), 
                        "heavy chain seqs for lacking light chain counterparts\n")
                    db_tmp_h = db_tmp_h[bool_common_h, ]
                }
                if (any(!bool_common_l)) {
                    cat("\n", subj, "- excluded", sum(!bool_common_l), 
                        "light chain seqs for lacking heavy chain counterparts\n")
                    db_tmp_l = db_tmp_l[bool_common_l, ]
                }
                stopifnot(nrow(db_tmp_h)==nrow(db_tmp_l))
                
                db_tmp = rbind(db_tmp_h, db_tmp_l)
            
            } else {
                # heavy only
                
                # .tsv
                #db_tmp = read.table(subj_info[["path_db_heavy"]][i],
                #                    header=T, sep="\t", stringsAsFactors=F)
                
                # .RData
                load(subj_info[["path_db_heavy"]][i])
                db_tmp = db; rm(db)
            }
            
            stopifnot(all(cols_keep %in% colnames(db_tmp)))
            db_tmp = db_tmp[, cols_keep]
            db_tmp[[opt$colSubj]] = subj_info[["subj"]][i]
            
            # can't name the combined data.frame `db` within this for loop
            # because of `rm(db)` when loading .RData files
            if (i==1) {
                # initiate db_all
                db_all = db_tmp
            } else {
                # append new rows to db_all
                db_all = rbind(db_all, db_tmp)
            }
            
            rm(db_tmp)
        }
        
        # rename
        db = db_all; rm(db_all)
        
        cat("\nBreakdown by subject:\n")
        print(table(db[[opt$colSubj]]))
        
        cat("\nCalculating between-subject distToNearest... \n")
        
        fn = paste0("dtn", out_suffix, "_btwSubj",
                    ifelse(is.null(opt$subsampleBetween), "", 
                           paste0("_subsample-", opt$subsampleBetween)),
                    ".RData")
        
        #* v1.0.2 stable release requires a bug fix from commit 47bffe0
        #* this bug fix is packed into julianqz/wu_cimm:main_0.1.1
        #* but if not using that docker image, v1.0.2 will fail next block
        if (opt$heavyLight) {
            # heavy and light
            db = distToNearest(db, 
                               sequenceColumn=opt$colSeq,
                               vCallColumn=opt$colV,
                               jCallColumn=opt$colJ,
                               cellIdColumn=opt$colCell,
                               locusColumn=opt$colLocus,
                               onlyHeavy=F,
                               cross=opt$colSubj,
                               subsample=opt$subsampleBetween,
                               model="ham", normalize="len", 
                               first=F, VJthenLen=F, keepVJLgroup=T,
                               nproc=opt$nproc, progress=T)
            
            # sanity check
            # heavy and light chains from the same cell should have the same partition
            stopifnot(all( sapply(unique(db[[opt$colCell]]), 
                                  function(s){
                                      # wrt db
                                      s_idx = which(db[[opt$colCell]]==s)
                                      s_bool = length(unique(db[["vjl_group"]][s_idx]))==1
                                      return(s_bool)
                                  }, USE.NAMES=F) ))
        } else {
            # heavy only
            db = distToNearest(db, 
                               sequenceColumn=opt$colSeq,
                               vCallColumn=opt$colV,
                               jCallColumn=opt$colJ,
                               cross=opt$colSubj,
                               subsample=opt$subsampleBetween,
                               model="ham", normalize="len", 
                               first=F, VJthenLen=F, keepVJLgroup=T,
                               nproc=opt$nproc, progress=T)
        }
        
        save(db, file=fn)
        
        rm(db, fn)
        
    } else {
        cat("\nOnly 1 subject found. Skipping between-subject.\n")
    }
    
}

