#!/opt/conda/bin/Rscript

# Julian Q. Zhou
# https://github.com/julianqz

# wrapper to partition sequences based on VJL combinations via `alakazam::groupGenes`
# - either based on heavy chains only
# - or based on both heavy and light chains

# assumes:
# - pathCSV points to a comma-separated file with the following headers
#   "subj", "path_db_heavy/light" where "path_db_*" points to a .tsv file
# - note that heavy and light chains are expected to be in separate files

# If the data is single-cell heavy:light paired, but partitioning based on
# heavy chains only is desired, simply set `--heavyLight` to `FALSE` and 
# light chains will be ignored.

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--pathCSV", action="store", default=NA, 
                type="character", help="Path to CSV containing subject list and paths to input files."),
    make_option("--pathWork", action="store", default=NA, type="character", 
                help="path_work."),
    make_option("--calcWithin", action="store", default=FALSE, type="logical", 
                help="Whether to calculate within-subject dtn."),
    make_option("--calcBetween", action="store", default=FALSE, type="logical", 
                help="Whether to calculate bewteen-subject dtn."),
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

suppressPackageStartupMessages(library(alakazam))

setwd(opt$pathWork)
sinkName = paste0("computingEnv_partition_", Sys.Date(), "-", 
                  format(Sys.time(), "%H%M%S"), '.txt')
sink(sinkName)
cat("Partition using both heavy and light:", opt$heavyLight, "\n")
cat("calcWithin:", opt$calcWithin, "\n")
cat("calcBetween:", opt$calcBetween, "\n")
sessionInfo()
sink()


#### within-subject ####

if (opt$calcWithin) {
    
    cat("\nPerforming within-subject partitioning... \n")
    
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
        
        # use `dtn` in filename because `bcr_infer_clone_wrapper.R` expects so
        fn = paste0("dtn", out_suffix, "_", subj, 
                    ".RData")
        
        # temporary column
        col_tmp_seq_len = "tmp_seq_len"
        db[[col_tmp_seq_len]] = nchar(db[[opt$colSeq]])
        
        # adds $vj_group columns
        # even if `junc_len` specified, added column is named `vj_group`
        
        if (opt$heavyLight) {
            # heavy and light
            db = groupGenes(data=db,
                            v_call=opt$colV,
                            j_call=opt$colJ,
                            junc_len=col_tmp_seq_len,
                            cell_id=opt$colCell,
                            locus=opt$colLocus,
                            only_heavy=F,
                            first=F)
        
            # sanity check
            # heavy and light chains from the same cell should have the same partition
            stopifnot(all( sapply(unique(db[[opt$colCell]]), 
                                  function(s){
                                      # wrt db
                                      s_idx = which(db[[opt$colCell]]==s)
                                      s_bool = length(unique(db[["vj_group"]][s_idx]))==1
                                      return(s_bool)
                                  }, USE.NAMES=F) ))
            # all counts of vj_group should be even (heavy:light paired)
            stopifnot(all( table(db[["vj_group"]]) %% 2 == 0 ))
            
        } else {
            # heavy only
            db = groupGenes(data=db,
                            v_call=opt$colV,
                            j_call=opt$colJ,
                            junc_len=col_tmp_seq_len,
                            cell_id=NULL,
                            first=F)
        }
        
        # rename `vj_group`; remove `vj_group` after renaming
        stopifnot("vj_group" %in% colnames(db))
        if ("vjl_group" %in% colnames(db)) {
            warning("A `vjl_group` column already exists; it will be overwritten.")
        }
        db[["vjl_group"]] = db[["vj_group"]]
        db[["vj_group"]] = NULL
        
        # remove temporary column
        db[[col_tmp_seq_len]] = NULL
        
        nrow_af = nrow(db)
        # sanity check
        # number of rows should remain the same
        # unless a row contained NA in any of v_call, j_call, or junc_len
        if (nrow_bf!=nrow_af) {
            warning("nrow_bf (", nrow_bf, ") != nrow_af (", nrow_af, ")\n")
        }
        
        save(db, file=fn)

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
            cols_keep = c(opt$colSeqID, opt$colSeq, opt$colV, opt$colJ,
                          opt$colLocus)
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
        
        nrow_bf = nrow(db)
        cat("\nnrow(db):", nrow_bf, "\n")
        
        cat("\nBreakdown by subject:\n")
        print(table(db[[opt$colSubj]]))
        
        cat("\nPerforming between-subject partitioning... \n")
        
        # use `dtn` in filename because `bcr_infer_clone_wrapper.R` expects so
        fn = paste0("dtn", out_suffix, "_btwSubj",
                    ".RData")
        
        # temporary column
        col_tmp_seq_len = "tmp_seq_len"
        db[[col_tmp_seq_len]] = nchar(db[[opt$colSeq]])
        
        # adds $vj_group columns
        # even if `junc_len` specified, added column is named `vj_group`
        
        if (opt$heavyLight) {
            # heavy and light
            db = groupGenes(data=db,
                            v_call=opt$colV,
                            j_call=opt$colJ,
                            junc_len=col_tmp_seq_len,
                            cell_id=opt$colCell,
                            locus=opt$colLocus,
                            only_heavy=F,
                            first=F)
            
            # sanity check
            # heavy and light chains from the same cell should have the same partition
            stopifnot(all( sapply(unique(db[[opt$colCell]]), 
                                  function(s){
                                      # wrt db
                                      s_idx = which(db[[opt$colCell]]==s)
                                      s_bool = length(unique(db[["vj_group"]][s_idx]))==1
                                      return(s_bool)
                                  }, USE.NAMES=F) ))
            # all counts of vj_group should be even (heavy:light paired)
            stopifnot(all( table(db[["vj_group"]]) %% 2 == 0 ))
            
        } else {
            # heavy only
            db = groupGenes(data=db,
                            v_call=opt$colV,
                            j_call=opt$colJ,
                            junc_len=col_tmp_seq_len,
                            cell_id=NULL,
                            first=F)
        }
        
        # rename `vj_group`; remove `vj_group` after renaming
        stopifnot("vj_group" %in% colnames(db))
        if ("vjl_group" %in% colnames(db)) {
            warning("A `vjl_group` column already exists; it will be overwritten.")
        }
        db[["vjl_group"]] = db[["vj_group"]]
        db[["vj_group"]] = NULL
        
        # remove temporary column
        db[[col_tmp_seq_len]] = NULL
        
        nrow_af = nrow(db)
        # sanity check
        # number of rows should remain the same
        # unless a row contained NA in any of v_call, j_call, or junc_len
        if (nrow_bf!=nrow_af) {
            warning("nrow_bf (", nrow_bf, ") != nrow_af (", nrow_af, ")\n")
        }
        
        save(db, file=fn)
        
        rm(db, fn)
        
    } else {
        cat("\nOnly 1 subject found. Skipping between-subject.\n")
    }
    
}

