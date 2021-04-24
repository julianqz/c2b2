#!/opt/conda/bin/Rscript

# helpful: http://www.cureffi.org/2014/01/15/running-r-batch-mode-linux/

# wrapper to call perform_qc and split_db from command line via Rscript
# some parameters for the functions are hard-coded as AIRR-C column names

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--helper", action="store", default=NA, type="character", 
                help="Path to bcr_qc.R."),
    make_option("--qc", action="store", default=FALSE, type="logical", 
                help="Boolean. Whether to perform QC."),
    make_option("--qcSeq", action="store", default=FALSE, type="logical", 
                help="Boolean. Whether to perform sequence-level QC."),
    make_option("--qcCell", action="store", default=FALSE, type="logical", 
                help="Boolean. Whether to perform cell-level QC."),
    make_option("--qcDb", action="store", default=NA, type="character", 
                help="db_name for perform_qc."),
    make_option("--qcOutname", action="store", default=NA, type="character", 
                help="[outname]_qc.tsv."),
    make_option("--qcOutdir", action="store", default=NA, type="character", 
                help="Path to write [outname]_qc.tsv."),
    make_option("--qcMaxPercN", action="store", default=10, type="numeric", 
                help="max_perc_N."),
    make_option("--qcColPercN", action="store", default="sequence_alignment", type="character", 
                help="col_perc_N."),
    make_option("--qcMaxNumNonATGCN", action="store", default=4, type="numeric", 
                help="max_num_nonATGCN."),
    make_option("--qcColNoneEmpty", action="store", default="germline_alignment", type="character", 
                help="col_none_empty."),
    make_option("--qcColNA", action="store", 
                default="PRCONS, CREGION, germline_alignment, junction, productive", 
                type="character", 
                help="col_NA."),
    make_option("--sp", action="store", default=FALSE, type="logical", 
                help="Boolean. Whether to split db."),
    make_option("--spDb", action="store", default=NA, type="character", 
                help="db_name for split_db."),
    make_option("--spOutname", action="store", default=NA, type="character", 
                help="[outname]_[heavy|light]_[pr|npr].tsv."),
    make_option("--spOutdir", action="store", default=NA, type="character", 
                help="Path to write [outname]_[heavy|light]_[pr|npr].tsv.")
)
opt = parse_args(OptionParser(option_list=option_list))

# source helper functions
source(opt$helper)

# parse 

# \s is space
# ? means preceding item is optional and will be matched at most once
col_perc_N = strsplit(opt$qcColPercN, "\\s?,\\s?")[[1]]
col_none_empty = strsplit(opt$qcColNoneEmpty, "\\s?,\\s?")[[1]]
col_NA = strsplit(opt$qcColNA, "\\s?,\\s?")[[1]]

# for debugging
#cat("col_perc_N:", col_perc_N, "; len:", length(col_perc_N), "\n")
#cat("col_none_empty:", col_none_empty, "; len:", length(col_none_empty), "\n")
#cat("col_NA:", col_NA, "; len:", length(col_NA), "\n")

# perform QC

if (opt$qc) {
    
    perform_qc(db_name=opt$qcDb, seq_level=opt$qcSeq, cell_level=opt$qcCell, 
               outname=opt$qcOutname, outdir=opt$qcOutdir,
               chain_type="IG",
               col_v_call="v_call", col_d_call="d_call", 
               col_j_call="j_call", col_c_call="c_gene",
               check_valid_vj=T, 
               check_chain_consistency=T, 
               check_perc_N=T, max_perc_N=opt$qcMaxPercN, 
               col_perc_N=col_perc_N, last_pos_N=312,
               check_num_nonATGCN=T, col_obsv="sequence_alignment", 
               col_germ="germline_alignment",
               max_num_nonATGCN=opt$qcMaxNumNonATGCN, last_pos_nonATGCN=312,
               check_none_empty=T, 
               col_none_empty=col_none_empty,
               check_NA=T, 
               col_NA=col_NA,
               check_len_mod3=T, col_len_mod3="junction")
    
}

# split db

if (opt$sp) {
    
    split_db(db_name=opt$spDb, col_v_call="v_call", 
             col_prod="productive", val_prod=TRUE, 
             outname=opt$spOutname, outdir=opt$spOutdir)
    
}
