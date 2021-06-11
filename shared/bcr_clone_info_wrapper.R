#!/opt/conda/bin/Rscript

# wrapper to summarize clonal compositions

# assumes:
# - pathCSV points to a comma-separated file with the following headers
#   "subj", "path_db_heavy" where "path_db_heavy" points to a .RData file

# intended to run based on heavy chains only

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--pathHelper", action="store", default=NA, 
                type="character", help="Path to bcr_collapse_duplicates."),
    make_option("--pathCSV", action="store", default=NA, 
                type="character", help="Path to CSV containing subject list and paths to input files."),
    make_option("--pathWork", action="store", default=NA, type="character", 
                help="path_work."),
    make_option("--colClone", action="store", default="clone_id", 
                type="character", help="col_clone."),
    make_option("--colVec", action="store", default=NA, 
                type="character", help="Parsed into col_vec."),
    make_option("--order", action="store", default="decreasing", 
                type="character", help="order_by_size.")
)
opt = parse_args(OptionParser(option_list=option_list))

# summarize_clone
source(opt$pathHelper)

subj_info = read.table(opt$pathCSV, header=T, sep=",", stringsAsFactors=F)

# parse
# \s is space
# ? means preceding item is optional and will be matched at most once
col_vec = strsplit(opt$colVec, "\\s?,\\s?")[[1]]

setwd(opt$pathWork)

for (i in 1:nrow(subj_info)) {
    
    subj = subj_info[["subj"]][i]
    cat("\n", subj, "\n")
    
    # load db
    load(subj_info[["path_db_heavy"]][i])
    
    # summarize
    clone_info = summarize_clone(db, 
                                 col_clone=opt$colClone,
                                 col_vec=col_vec,
                                 order_by_size=opt$order)
    
    # export
    fn = paste0("clone_info_", subj, ".RData")
    save(clone_info, file=fn)
    
    rm(db, clone_info)
}

