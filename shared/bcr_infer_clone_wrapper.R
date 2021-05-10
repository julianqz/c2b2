#!/opt/conda/bin/Rscript

# wrapper to infer B cell clones for a list of individuals

# assumes:
# - alakazam::groupGene has been performed in the appropriate mode (bulk/single-cell)
#   and a VJL group column has been created in `db`
# - pathCSV points to a comma-separated file with the following headers
#   "subj", "path_db" where "path_db" points to a .RData file containing `db`

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--pathCSV", action="store", default=NA, type="character", 
                help="Path to CSV containing subject list and paths to input files."),
    make_option("--pathWork", action="store", default=NA, type="character", 
                help="Path to save outputs."),
    make_option("--pathHelper", action="store", default=NA, type="character", 
                help="Path to bcr_infer_clone.R."),
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
                type="character", help="linkage.")
)
opt = parse_args(OptionParser(option_list=option_list))

source(opt$pathHelper)

subj_info = read.table(opt$pathCSV, header=T, sep=",", stringsAsFactors=F)

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
cat("threshold:", opt$threshold, "\n")
sessionInfo()
sink()

for (i in 1:nrow(subj_info)) {
    
    subj = subj_info[["subj"]][i]
    thresh = thresh_vec[subj]
    cat(subj, thresh, "\n")
    
    # load db
    fn_in = subj_info[["path_db"]][i]
    load(fn_in)
    
    # infer clones
    infer_lst = defineCloneDb(db, 
                              sequenceColumn=opt$colSeq,
                              VJLgroupColumn=opt$colVJLgroup,
                              cloneColumn=opt$colClone,
                              threshold=thresh,
                              maxmiss=opt$maxmiss,
                              linkage=opt$linkage, 
                              verbose=T)
    
    db_clust = infer_lst[["db_clust"]]
    db_fail = infer_lst[["db_fail"]]
    rm(infer_lst, db)
    
    # export
    if (!is.null(db_clust)) {
        fn_out_pass_tsv = paste0("cluster-pass_", subj, ".tsv")
        fn_out_pass_r = paste0("cluster-pass_", subj, ".RData")
        write.table(x=db_clust, file=fn_out_pass_tsv, quote=F, sep="\t",
                    row.names=F, col.names=T)
        save(db_clust, file=fn_out_pass_r)
    }
    
    if (!is.null(db_fail)) {
        fn_out_fail_tsv = paste0("cluster-fail_", subj, ".tsv")
        fn_out_fail_r = paste0("cluster-fail_", subj, ".RData")
        write.table(x=db_fail, file=fn_out_fail_tsv, quote=F, sep="\t",
                    row.names=F, col.names=T)
        save(db_fail, file=fn_out_fail_r)
    }
    
}
