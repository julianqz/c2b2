#!/opt/conda/bin/Rscript

# Julian Q. Zhou
# https://github.com/julianqz

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
    make_option("--qcSequential", action="store", default=FALSE, type="logical", 
                help="Boolean. Whether to perform sequence- and cell-level QCs sequentially."),
    make_option("--qcDb", action="store", default=NA, type="character", 
                help="db_name for perform_qc."),
    make_option("--qcOutname", action="store", default=NA, type="character", 
                help="[outname]_qc.tsv."),
    make_option("--qcOutdir", action="store", default=NA, type="character", 
                help="Path to write [outname]_qc.tsv."),
    make_option("--qcColV", action="store", default="v_call", type="character", 
                help="col_v_call"),
    make_option("--qcColD", action="store", default="d_call", type="character", 
                help="col_d_call."),
    make_option("--qcColJ", action="store", default="j_call", type="character", 
                help="col_j_call."),
    make_option("--qcColC", action="store", default=NA, type="character", 
                help="col_c_call."),
    make_option("--qcColObsv", action="store", default="sequence_alignment", 
                type="character", 
                help="col_obsv."),
    make_option("--qcColGerm", action="store", default="germline_alignment", 
                type="character", 
                help="col_germ."),
    make_option("--qcMaxN", action="store", default=10, type="numeric", 
                help="max_N."),
    make_option("--qcColN", action="store", default="cdr3", 
                type="character", 
                help="col_N."),
    # Set type fr --qcLastPosN as character so that the command line input can be
    # parsed (this enables multiple values separated by "," to be inputted)
    # After parsing, the values are typecast to integer
    make_option("--qcLastPosN", action="store", default=NA, type="character", 
                help="last_pos_N."),
    make_option("--qcAsPercN", action="store", default=FALSE, type="logical", 
                help="as_perc_N."),
    make_option("--qcMaxNonATGC", action="store", default=4, type="numeric", 
                help="max_nonATGC."),
    make_option("--qcAsPercNonATGC", action="store", default=FALSE, type="logical", 
                help="as_perc_nonATGC."),
    make_option("--qcColNoneEmpty", action="store", default="germline_alignment", 
                type="character", 
                help="col_none_empty."),
    make_option("--qcColNA", action="store", 
                default="PRCONS, CREGION, germline_alignment, cdr3, productive", 
                type="character", 
                help="col_NA."),
    make_option("--qcColLenMod3", action="store", default="cdr3", 
                type="character", 
                help="col_len_mod3."),
    make_option("--qcColLocus", action="store", default=NA, 
                type="character", 
                help="col_locus."),
    make_option("--qcColCell", action="store", default=NA, 
                type="character", 
                help="col_cell."),
    make_option("--qcColUMI", action="store", default=NA, 
                type="character", 
                help="col_umi."),
    make_option("--qcLogicNumHL", action="store", default=NA, 
                type="character", 
                help="logic_num_HL."),
    make_option("--sp", action="store", default=FALSE, type="logical", 
                help="Boolean. Whether to split db."),
    make_option("--spDb", action="store", default=NA, type="character", 
                help="db_name for split_db."),
    make_option("--spOutname", action="store", default=NA, type="character", 
                help="[outname]_[heavy|light]_[pr|npr].tsv."),
    make_option("--spOutdir", action="store", default=NA, type="character", 
                help="Path to write [outname]_[heavy|light]_[pr|npr].tsv."),
    make_option("--spColV", action="store", default="v_call", type="character", 
                help="col_v_call."),
    make_option("--spColProd", action="store", default="productive", type="character", 
                help="col_prod."),
    make_option("--spValProd", action="store", default=TRUE, type="logical", 
                help="val_prod.")
)
opt = parse_args(OptionParser(option_list=option_list))

# source helper functions
source(opt$helper)

# parse 

# if --qcColC is set via `-D` in command line as `"NA"`
# is.na(opt$qcColC) will be TRUE

# \s is space
# ? means preceding item is optional and will be matched at most once
col_N = strsplit(opt$qcColN, "\\s?,\\s?")[[1]]
col_none_empty = strsplit(opt$qcColNoneEmpty, "\\s?,\\s?")[[1]]
col_NA = strsplit(opt$qcColNA, "\\s?,\\s?")[[1]]
col_len_mod3 = strsplit(opt$qcColLenMod3, "\\s?,\\s?")[[1]]

last_pos_N = as.integer(strsplit(opt$qcLastPosN, "\\s?,\\s?")[[1]])

# for debugging
#cat("col_N:", col_N, "; len:", length(col_N), "\n")
#cat("col_none_empty:", col_none_empty, "; len:", length(col_none_empty), "\n")
#cat("col_NA:", col_NA, "; len:", length(col_NA), "\n")

# perform QC

if (opt$qc) {
    
    perform_qc(db_name=opt$qcDb, seq_level=opt$qcSeq, cell_level=opt$qcCell, 
               sequential=opt$qcSequential,
               outname=opt$qcOutname, outdir=opt$qcOutdir,
               chain_type="IG",
               col_v_call=opt$qcColV, col_d_call=opt$qcColD, 
               col_j_call=opt$qcColJ, col_c_call=opt$qcColC,
               check_valid_vj=T, 
               check_chain_consistency=T, 
               check_N=T, max_N=opt$qcMaxN, 
               col_N=col_N, last_pos_N=last_pos_N, as_perc_N=opt$qcAsPercN,
               check_nonATGC=T, col_obsv=opt$qcColObsv, 
               col_germ=opt$qcColGerm,
               max_nonATGC=opt$qcMaxNonATGC, last_pos_nonATGC=312,
               as_perc_nonATGC=opt$qcAsPercNonATGC,
               check_none_empty=T, 
               col_none_empty=col_none_empty,
               check_NA=T, 
               col_NA=col_NA,
               check_len_mod3=T, col_len_mod3=col_len_mod3,
               col_locus=opt$qcColLocus, 
               col_cell=opt$qcColCell,
               col_umi=opt$qcColUMI,
               check_locus=T,
               check_num_HL=T, logic_num_HL=opt$qcLogicNumHL
               )
    
}

# split db

if (opt$sp) {
    
    split_db(db_name=opt$spDb, col_v_call=opt$spColV, 
             col_prod=opt$spColProd, val_prod=opt$spValProd, 
             outname=opt$spOutname, outdir=opt$spOutdir)
    
}
