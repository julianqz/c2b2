#!/opt/conda/bin/Rscript

# wrapper to call run_tigger for a list of individuals from command line via Rscript

# assumes:
# - all input files are in .RData format
# - all input files are in the same pathDb location
# - all input files are named as [prefix]_[subj].RData

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--subjList", action="store", default=NA, 
                type="character", help="Comma-separated list of subjects."),
    make_option("--prefix", action="store", default=NA, 
                type="character", 
                help="Common filename prefix to input .RData files. Eg: 'db_combined_'."),
    make_option("--pathDb", action="store", default=NA, type="character", 
                help="Common path to input .RData files."),
    make_option("--pathIMGT", action="store", default=NA, type="character", 
                help="path_imgt."),
    make_option("--pathHelper", action="store", default=NA, type="character", 
                help="path_helper."),
    make_option("--pathWork", action="store", default=NA, type="character", 
                help="path_work."),
    make_option("--colSeq", action="store", default="sequence_alignment", 
                type="character", help="col_seq."),
    make_option("--colV", action="store", default="v_call", 
                type="character", help="col_v."),
    make_option("--colProd", action="store", default="productive", 
                type="character", help="col_prod."),
    make_option("--colJuncLen", action="store", default="junction_length", 
                type="character", help="col_junc_len."),
    make_option("--findUnmutated", action="store", default=TRUE, 
                type="logical", help="p_find_unmutated."),
    make_option("--textSize", action="store", default=12, 
                type="numeric", help="p_text_size."),
)
opt = parse_args(OptionParser(option_list=option_list))

# parse

# \s is space
# ? means preceding item is optional and will be matched at most once
vec_subj = strsplit(opt$subjList, "\\s?,\\s?")[[1]]

library(tigger)

setwd(opt$pathWork)
sinkName = paste0("computingEnv_tigger_", Sys.Date(), "-", 
                  format(Sys.time(), "%H%M%S"), '.txt')
sink(sinkName)
sessionInfo()
sink()


for (subj in vec_subj) {
    
    # load db
    setwd(opt$pathDb)
    fn = paste0(opt$prefix, subj, ".RData")
    load(fn)
    
    run_tigger(path_imgt=opt$pathIMGT, 
               path_helper=opt$pathHelper, 
               path_work=opt$pathWork, 
               subj=subj, db=db, 
               col_seq=opt$colSeq, 
               col_v=opt$colV, 
               col_prod=opt$colProd, 
               col_junc_len=opt$colJuncLen,
               p_find_unmutated=opt$findUnmutated,
               p_text_size=opt$textSize)
    
    rm(fn, db)
}

