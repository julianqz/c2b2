#!/opt/conda/bin/Rscript

# wrapper to collapse duplicates

# In each clone,
# remove duplicate V(D)J sequences (e.g. $sequence_alignment), 
# optionally retaining duplicates derived from different subsets, 
# assigned to different isotypes, etc.

# assumes:
# - pathCSV points to a comma-separated file with the following headers
#   "subj", "path_db_heavy/light" where "path_db_*" points to a .tsv file

# Known issue:
# If data is unpaird bulk light, previous steps in current pipeline does not 
# add clone assignment to the light chains. However, --colClone is required 
# by this function. As a hack-around, manually add a dummy --colClone to the
# light db with all light chain seqs assigned to a single clone

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--pathHelper", action="store", default=NA, 
                type="character", help="Path to bcr_collapse_duplicates."),
    make_option("--pathCSV", action="store", default=NA, 
                type="character", help="Path to CSV containing subject list and paths to input files."),
    make_option("--pathWork", action="store", default=NA, type="character", 
                help="path_work."),
    make_option("--nproc", action="store", default=1, type="numeric", 
                help="nproc."),
    make_option("--heavyLight", action="store", default=FALSE, 
                type="logical", 
                help="Whether to run separately for both heavy and light chains. Default is heavy only."),
    make_option("--colClone", action="store", default="clone_id", 
                type="character", help="col_clone."),
    make_option("--colSeq", action="store", default="sequence_alignment", 
                type="character", help="col_seq."),
    make_option("--colID", action="store", default="", 
                type="character", help="col_id."),
    make_option("--colTextFields", action="store", default=NULL, 
                type="character", 
                help="Comma-separated string parsed into col_text_fields. Do not specify via command line if NULL."),
    make_option("--colNumFields", action="store", default=NULL, 
                type="character", 
                help="Comma-separated string parsed into col_num_fields. Do not specify via command line if NULL."),
    make_option("--colSeqFields", action="store", default=NULL, 
                type="character", 
                help="Comma-separated string parsed into col_seq_fields. Do not specify via command line if NULL."),
    make_option("--colPreserve", action="store", default=NULL, 
                type="character", 
                help="col_preserve. Do not specify via command line if NULL."),
    make_option("--valPreserveVec", action="store", default=NULL, 
                type="character", 
                help="Comma-separated string parsed into val_preserve_vec. Do not specify via command line if NULL."),
    make_option("--colDistinctVec", action="store", default=NULL, 
                type="character", 
                help="Comma-separated string parsed into col_distinct_vec. Do not specify via command line if NULL."),
    make_option("--verify", action="store", default=FALSE, type="logical", 
                help="Whether to perform verification."),
    make_option("--verifyN", action="store", default=1000, type="numeric", 
                help="The top N biggest clones on which verification is to be performed.")
)
opt = parse_args(OptionParser(option_list=option_list))

# without next line sink won't be able to capture alakazam version
# (even if func in opt$pathHelper `require(alakazam)`)
suppressPackageStartupMessages(require(alakazam))

source(opt$pathHelper)

subj_info = read.table(opt$pathCSV, header=T, sep=",", stringsAsFactors=F)

if (opt$heavyLight) {
    run_mode = c("heavy", "light")
} else {
    run_mode = c("heavy")
}

# check presence of necessary columns
stopifnot( all( paste0("path_db_", run_mode) %in% colnames(subj_info) ) )


# parse
# \s is space
# ? means preceding item is optional and will be matched at most once

if (!is.null(opt$colTextFields)) {
    col_text_fields = strsplit(opt$colTextFields, "\\s?,\\s?")[[1]]
} else {
    col_text_fields = NULL
}

if (!is.null(opt$colNumFields)) {
    col_num_fields = strsplit(opt$colNumFields, "\\s?,\\s?")[[1]]
} else {
    col_num_fields = NULL
}

if (!is.null(opt$colSeqFields)) {
    col_seq_fields = strsplit(opt$colSeqFields, "\\s?,\\s?")[[1]]
} else {
    col_seq_fields = NULL
}

if (!is.null(opt$valPreserveVec)) {
    val_preserve_vec = strsplit(opt$valPreserveVec, "\\s?,\\s?")[[1]]
} else {
    val_preserve_vec = NULL
}

if (!is.null(opt$colDistinctVec)) {
    col_distinct_vec = strsplit(opt$colDistinctVec, "\\s?,\\s?")[[1]]
} else {
    col_distinct_vec = NULL
}


setwd(opt$pathWork)
sinkName = paste0("computingEnv_collapse_duplicates_", Sys.Date(), "-", 
                  format(Sys.time(), "%H%M%S"), '.txt')
sink(sinkName)
cat("run_mode:", run_mode, "\n")
# if NULL, will appear as "...: " (i.e. blank)
cat("col_seq:", opt$colSeq, "\n")
cat("col_id:", opt$colID, "\n")
cat("col_text_fields:", col_text_fields, "\n")
cat("col_num_fields:", col_num_fields, "\n")
cat("col_seq_fields:", col_seq_fields, "\n")
cat("col_preserve:", opt$colPreserve, "\n")
cat("val_preserve_vec:", val_preserve_vec, "\n")
cat("col_distinct_vec:", col_distinct_vec, "\n")
cat("nproc:", opt$nproc, "\n")
sessionInfo()
sink()


for (i in 1:nrow(subj_info)) {
    
    subj = subj_info[["subj"]][i]
    
    for (cur_run_mode in run_mode) {
        
        cat("\n", subj, ";", cur_run_mode, "\n")
        
        #### collapse ####
        
        # load data
        cur_col_db = paste0("path_db_", cur_run_mode)
        cur_fn_db = subj_info[[cur_col_db]][i]
        cat(" - loading", cur_fn_db, "\n")
        db = read.table(cur_fn_db, 
                        sep="\t", header=T, stringsAsFactors=F)
        
        # collapse
        db = run_collapse_duplicates(db, nproc=opt$nproc, 
                                     col_clone=opt$colClone,
                                     col_seq=opt$colSeq,
                                     col_id=opt$colID,
                                     col_text_fields=col_text_fields,
                                     col_num_fields=col_num_fields,
                                     col_seq_fields=col_seq_fields,
                                     col_preserve=opt$colPreserve,
                                     val_preserve_vec=val_preserve_vec,
                                     col_distinct_vec=col_distinct_vec)
        
        # export
        fn = paste0("collapse_dups_", cur_run_mode, "_", subj, ".RData")
        save(db, file=fn)
        
        #### verify ####
        
        if (opt$verify) {
            
            verify_collapse_duplicates(db=db, N=opt$verifyN,
                                       col_clone=opt$colClone,
                                       col_seq=opt$colSeq,
                                       col_preserve=opt$colPreserve,
                                       val_preserve_vec=val_preserve_vec,
                                       col_distinct_vec=col_distinct_vec)
        }
        
        rm(db)
    }
}

