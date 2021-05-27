#!/opt/conda/bin/Rscript

# wrapper to call run_tigger for a list of individuals from command line via Rscript

# assumes:
# - all input files are in .RData format
# - all input files are in the same pathDb location
# - all input files are named as [prefix]_[subj].RData

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--pathRunTigger", action="store", default=NA, type="character", 
                help="Path to bcr_tigger.R."),
    make_option("--pathCSV", action="store", default=NA, 
                type="character", help="Path to CSV containing subject list and paths to input files."),
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
    make_option("--colCDR3Len", action="store", default="cdr3_length", 
                type="character", help="col_cdr3_len."),
    make_option("--findUnmutated", action="store", default=TRUE, 
                type="logical", help="p_find_unmutated."),
    make_option("--textSize", action="store", default=12, 
                type="numeric", help="p_text_size.")
)
opt = parse_args(OptionParser(option_list=option_list))

# without next line sink won't be able to capture tigger version
# (even if func in opt$pathRunTigger `require(tigger)`)
suppressPackageStartupMessages(library(tigger))

# run_tigger()
source(opt$pathRunTigger)

subj_info = read.table(opt$pathCSV, header=T, sep=",", stringsAsFactors=F)


setwd(opt$pathWork)
sinkName = paste0("computingEnv_tigger_", Sys.Date(), "-", 
                  format(Sys.time(), "%H%M%S"), '.txt')
sink(sinkName)
sessionInfo()
sink()


for (i in 1:nrow(subj_info)) {
    
    subj = subj_info[["subj"]][i]
    cat("\n", subj, "\n")
    
    # db
    load(subj_info[["path_db"]][i])
    
    run_tigger(path_imgt=opt$pathIMGT, 
               path_helper=opt$pathHelper, 
               path_work=opt$pathWork, 
               subj=subj, db=db, 
               col_seq=opt$colSeq, 
               col_v=opt$colV, 
               col_prod=opt$colProd, 
               col_cdr3_len=opt$colCDR3Len,
               p_find_unmutated=opt$findUnmutated,
               p_text_size=opt$textSize)

}

