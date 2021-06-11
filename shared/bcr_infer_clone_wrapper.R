#!/opt/conda/bin/Rscript

# wrapper to infer B cell clones for a list of individuals

# assumes:
# - alakazam::groupGene has been performed in the appropriate mode (bulk/single-cell)
#   and a VJL group column has been created in `db`
# - pathCSV points to a comma-separated file with the following headers
#   -- "subj"
#   -- "path_dtn" containing a dtn .RData file containing `db`
#   -- if requiring propagation to light chains, "path_db_light" pointing to light chain tsv

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--pathCSV", action="store", default=NA, type="character", 
                help="Path to CSV containing subject list and paths to input files."),
    make_option("--pathWork", action="store", default=NA, type="character", 
                help="Path to save outputs."),
    make_option("--pathHelper", action="store", default=NA, type="character", 
                help="Path to bcr_infer_clone.R."),
    make_option("--pathAux", action="store", default=NA, type="character", 
                help="Path to directory to output subj_info_createGermlines.csv."),
    make_option("--threshold", action="store", default=NA, 
                type="character", help="Clustering threshold."),
    make_option("--colSeq", action="store", default="cdr3", 
                type="character", help="sequenceColumn."),
    make_option("--colVJLgroup", action="store", default="vjl_group", 
                type="character", help="VJLgroupColumn."),
    make_option("--colClone", action="store", default="clone_id", 
                type="character", help="cloneColumn."),
    make_option("--maxmiss", action="store", default=0, 
                type="numeric", help="maxmiss."),
    make_option("--linkage", action="store", default="single", 
                type="character", help="linkage."),
    make_option("--heavyLight", action="store", default=FALSE, type="logical", 
                help="Whether partition was based on both heavy and light chains."),
    make_option("--propagateToLight", action="store", default=FALSE, type="logical", 
                help="Whether to propagate clone assignment based on heavy chain CDR3 to light chains (regardless of whether partition was based on heavy only or both heavy and light chains)."),
    make_option("--colLocus", action="store", default="locus", type="character", 
                help="locusColumn. Ignored if --heavyLight FALSE."),
    make_option("--colCell", action="store", default="cell_id", type="character", 
                help="cellCColumn. Ignored if --heavyLight FALSE.")
)
opt = parse_args(OptionParser(option_list=option_list))

source(opt$pathHelper)

subj_info = read.table(opt$pathCSV, header=T, sep=",", stringsAsFactors=F)

# data contains heavy only;    partition based on heavy only;    no need to propagate to light
# data contains heavy & light; partition based on heavy only;    choice of propagating to light
# data contains heavy & light; partition based on heavy & light; choice of propagating to light

# in all cases, must specify path to dtn .RData
stopifnot("path_dtn" %in% colnames(subj_info))

if (opt$heavyLight) {
    # partition based on both heavy and light
    # dtn .RData contains both heavy and light
    out_suffix = "_groupByHL"
    out_suffix_2 = "HL"
} else {
    # partition based on heavy only
    out_suffix = "_groupByHonly"
    out_suffix_2 = "Honly"
}

if (opt$propagateToLight) {
    # must specify path to light chain db
    stopifnot("path_db_light" %in% colnames(subj_info))
}

col_locus = opt$colLocus
col_cell = opt$colCell

# parse
# \s is space
# ? means preceding item is optional and will be matched at most once

# first parse into strings
thresh_str = strsplit(opt$threshold, "\\s?,\\s?")[[1]]
# check length
stopifnot( length(thresh_str) == nrow(subj_info) )
# convert into numeric
thresh_vec = as.numeric(thresh_str)
names(thresh_vec) = subj_info[["subj"]]
# sanity check
.check_thresh = function(v) { !is.na(v) & is.numeric(v) & v>=0 & v<=1 }
stopifnot( all(sapply(thresh_vec, .check_thresh)) )


setwd(opt$pathWork)
sinkName = paste0("computingEnv_cluster_", Sys.Date(), "-", 
                  format(Sys.time(), "%H%M%S"), '.txt')
sink(sinkName)
cat("Partition was based on both heavy and light:", opt$heavyLight, "\n")
cat("Propagate clone assignment based on heavy chain CDR3 to light chains:",
    opt$propagateToLight, "\n")
cat("threshold:", opt$threshold, "\n")
sessionInfo()
sink()

for (i in 1:nrow(subj_info)) {
    
    subj = subj_info[["subj"]][i]
    thresh = thresh_vec[subj]
    cat(subj, thresh, "\n")
    
    # load `db` from dtn .RData
    fn_in = paste0(subj_info[["path_dtn"]][i],
                   "dtn", out_suffix, "_", subj, ".RData")
    load(fn_in)
    
    # if partition was based on both heavy and light
    if (opt$heavyLight) {
        # check columns
        stopifnot(all( c(col_locus, col_cell) %in% colnames(db) ))
        # expect db to contain both heavy and light
        stopifnot( sum(c("IGH","IGL","IGK") %in% db[[col_locus]])>=2 )
        # subset to heavy 
        # (defineCloneDb expects heavy only)
        db = db[db[[col_locus]]=="IGH", ]
    }
    
    # infer clones
    # Within each partition (based on either heavy only, or both heavy & light),
    # cluster using heavy chain CDR3
    infer_lst = defineCloneDb(db, 
                              sequenceColumn=opt$colSeq,
                              VJLgroupColumn=opt$colVJLgroup,
                              cloneColumn=opt$colClone,
                              threshold=thresh,
                              maxmiss=opt$maxmiss,
                              linkage=opt$linkage, 
                              verbose=T)
    
    db_heavy_clust = infer_lst[["db_clust"]]
    db_heavy_fail = infer_lst[["db_fail"]]
    # don't remove db yet (needed for propagating to light)
    rm(infer_lst)
    
    # export heavy
    if (!is.null(db_heavy_clust)) {
        fn_out_pass_base = paste0("cluster-pass_partition-", out_suffix_2, 
                                  "_chain-heavy_", subj)
        fn_out_pass_tsv = paste0(fn_out_pass_base, ".tsv")
        fn_out_pass_r = paste0(fn_out_pass_base, ".RData")
        write.table(x=db_heavy_clust, file=fn_out_pass_tsv, quote=F, sep="\t",
                    row.names=F, col.names=T)
        save(db_heavy_clust, file=fn_out_pass_r)
    }
    
    if (!is.null(db_heavy_fail)) {
        fn_out_fail_base = paste0("cluster-fail_partition-", out_suffix_2, 
                                  "_chain-heavy_", subj)
        fn_out_fail_tsv = paste0(fn_out_pass_base, ".tsv")
        fn_out_fail_r = paste0(fn_out_pass_base, ".RData")
        write.table(x=db_heavy_fail, file=fn_out_fail_tsv, quote=F, sep="\t",
                    row.names=F, col.names=T)
        save(db_heavy_fail, file=fn_out_fail_r)
    }
    
    
    # if applicable, propagate clone assignment to light chains
    # using cell id as the common link
    if (opt$propagateToLight) {
        
        # load light
        db_light = read.table(subj_info[["path_db_light"]][i],
                              sep="\t", header=T, stringsAsFactors=F)
        
        # check columns exist
        stopifnot(all( c(col_locus, col_cell) %in% colnames(db_light) ))
        
        # add columns
        # dist_nearest, vjl_group (pass & fail)
        # opt$colClone (pass only)
        
        # dist_nearest only calculated for heavy (regardless of partitioning)
        # need a placeholder column in order to perform rbind etc.
        db_light[["dist_nearest"]] = NA
        
        if (opt$heavyLight) {
            # partitioning was based on heavy & light
            # makes sense to propagate to light
            
            # wrt db (not db_heavy_clust/fail)
            idx_cell_wrt_db = match(db_light[[col_cell]], db[[col_cell]])
            stopifnot(all.equal( db_light[[col_cell]], 
                                 db[[col_cell]][idx_cell_wrt_db] ))
            
            db_light[["vjl_group"]] = db[["vjl_group"]][idx_cell_wrt_db] 
        } else {
            # partitioning was based on heavy only
            # doesn't make sense to propagate to light
            # but need a placeholder column in order to perform rbind etc.
            db_light[["vjl_group"]] = NA
        }
        
        bool_pass = db_light[[col_cell]] %in% db_heavy_clust[[col_cell]]
        
        if (any(bool_pass)) {
            # if any passed
            db_light_clust = db_light[bool_pass, ]
            
            # wrt db_heavy_clust (not db)
            idx_cell_wrt_db_heavy_clust = match(db_light_clust[[col_cell]], 
                                                db_heavy_clust[[col_cell]])
            stopifnot(all.equal( db_light_clust[[col_cell]], 
                                 db_heavy_clust[[col_cell]][idx_cell_wrt_db_heavy_clust] ))
            
            db_light_clust[[opt$colClone]] = db_heavy_clust[[opt$colClone]][idx_cell_wrt_db_heavy_clust] 
            
            # check same columns exist in light as in heavy
            stopifnot(all.equal( colnames(db_light_clust), colnames(db_heavy_clust) ))
            
            # match order of cells in heavy
            # wrt db_light_clust
            idx_cell_wrt_db_light_clust = match(db_heavy_clust[[col_cell]], 
                                                db_light_clust[[col_cell]])
            stopifnot(all.equal( db_heavy_clust[[col_cell]],
                                 db_light_clust[[col_cell]][idx_cell_wrt_db_light_clust] ))
            
            db_light_clust = db_light_clust[idx_cell_wrt_db_light_clust, ]
            rownames(db_light_clust) = NULL
            
            cols_ck = c(col_cell, "vjl_group", opt$colClone)
            
            stopifnot(all.equal(db_heavy_clust[, cols_ck], 
                                db_light_clust[, cols_ck],
                                check.attributes=F))
            
            # export
            fn_out_pass_base = paste0("cluster-pass_partition-", out_suffix_2, 
                                      "_chain-light_", subj)
            fn_out_pass_tsv = paste0(fn_out_pass_base, ".tsv")
            fn_out_pass_r = paste0(fn_out_pass_base, ".RData")
            write.table(x=db_light_clust, file=fn_out_pass_tsv, quote=F, sep="\t",
                        row.names=F, col.names=T)
            save(db_light_clust, file=fn_out_pass_r)
        }
        
        if (any(!bool_pass)) {
            # if any failed
            db_light_fail = db_light[!bool_pass, ]
            
            # check same columns exist in light as in heavy
            stopifnot(all.equal( colnames(db_light_fail), colnames(db_heavy_fail) ))
            
            # match order of cells in heavy
            # wrt db_light_fail
            idx_cell_wrt_db_light_fail = match(db_heavy_fail[[col_cell]], 
                                               db_light_fail[[col_cell]])
            db_light_fail = db_light_fail[idx_cell_wrt_db_light_fail, ]
            rownames(db_light_fail) = NULL
            
            cols_ck = c(col_cell, "vjl_group")
            stopifnot(all.equal(db_heavy_fail[, cols_ck], 
                                db_light_fail[, cols_ck],
                                check.attributes=F))
            
            # export
            fn_out_fail_base = paste0("cluster-fail_partition-", out_suffix_2, 
                                      "_chain-light_", subj)
            fn_out_fail_tsv = paste0(fn_out_pass_base, ".tsv")
            fn_out_fail_r = paste0(fn_out_pass_base, ".RData")
            write.table(x=db_light_fail, file=fn_out_fail_tsv, quote=F, sep="\t",
                        row.names=F, col.names=T)
            save(db_light_fail, file=fn_out_fail_r)
        }
    }
    
}

# output subj_info_createGermlines.csv for CreateGermlines

cur_csv = cbind(subj=DONORS,
                path_db_heavy=paste0(opt$pathWork, 
                                     "cluster-pass_partition-", out_suffix_2, 
                                     "_chain-heavy_", DONORS, ".tsv"))
if (opt$propagateToLight) {
    cur_csv = cbind(cur_csv,
                    path_db_light=paste0(opt$pathWork, 
                                         "cluster-pass_partition-", out_suffix_2, 
                                         "_chain-light_", DONORS, ".tsv"))
}
cur_fn = "subj_info_createGermlines.csv"
setwd(opt$pathAux)
# note that headers are omitted for this file (col.names=F)
write.table(x=cur_csv, file=cur_fn, quote=F, sep=",", row.names=F, col.names=F)

