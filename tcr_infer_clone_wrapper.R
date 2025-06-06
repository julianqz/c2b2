#!/opt/conda/bin/Rscript

# Julian Q. Zhou
# https://github.com/julianqz

# wrapper to infer T cell clones for a list of individuals

# assumes:
# - If single-cell, data has been QC'ed so that each cell has 1 VDJ and 1 VJ  
# - pathCSV points to a comma-separated file with the following headers
#   -- "subj"
#   -- "path_vdj" containing a .tsv or .RData (db) file containing VDJ seqs
#   -- "path_vj" containing a .tsv or .RData (db) file containing VJ seqs [single-cell mode]

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--pathCSV", action="store", default=NA, type="character", 
                help="Path to CSV containing subject list and paths to input files."),
    make_option("--pathWork", action="store", default=NA, type="character", 
                help="Path to save outputs."),
    make_option("--pathHelper", action="store", default=NA, type="character", 
                help="Path to tcr_infer_clone.R."),
    make_option("--pathAux", action="store", default=NA, type="character", 
                help="Path to directory to output subj_info_createGermlines.csv."),
    make_option("--colV", action="store", default="v_call", 
                type="character", help="v_call."),
    make_option("--colJ", action="store", default="j_call", 
                type="character", help="j_call."),
    make_option("--colJuncLen", action="store", default="cdr3_length", 
                type="character", help="junc_len."),
    make_option("--colJunc", action="store", default="cdr3", 
                type="character", help="junc."),
    make_option("--colCell", action="store", default="cell_id", type="character", 
                help="cell_id. Ignored if --singleCellMode FALSE."),
    make_option("--colLocus", action="store", default="locus", type="character", 
                help="locus. Ignored if --singleCellMode FALSE."),
    make_option("--colClone", action="store", default="clone_id", 
                type="character", help="clone_id."),
    make_option("--singleCellMode", action="store", default=FALSE, type="logical", 
                help="Whether performing inference in single-cell mode."),
    make_option("--useOnlyVDJ", action="store", default=FALSE, type="logical", 
                help="Whether partition was based on both VDJ and VJ TCR chains."),
    make_option("--createGermlinesCSV", action="store", default=TRUE, type="logical",
                help="Whether to output auxiliary CSV for running CreateGermlines wrapper.")
)
opt = parse_args(OptionParser(option_list=option_list))

col_v_call = opt$colV
col_j_call = opt$colJ
col_junc_len = opt$colJunclen
col_junc = opt$colJunc
col_locus = opt$colLocus
col_cell_id = opt$colCell
col_clone_id = opt$colClone
single_cell_mode = opt$singleCellMode
use_only_vdj = opt$useOnlyVDJ

source(opt$pathHelper)

subj_info = read.table(opt$pathCSV, header=T, sep=",", stringsAsFactors=F)

# in all cases, must specify path to VDJ .RData
stopifnot("path_vdj" %in% colnames(subj_info))

if (single_cell_mode) {
    
    stopifnot(all(c("path_vdj", "path_vj") %in% colnames(subj_info)))
    
    if (use_only_vdj) {
        # partition based on VDJ only
        out_suffix = "onlyVDJ"
    } else {
        # partition based on both VDJ and VJ chains
        out_suffix = "VDJandVJ"
    }
    
} else {
    # partition based on VDJ only
    out_suffix = "onlyVDJ"
}


setwd(opt$pathWork)
sinkName = paste0("computingEnv_cluster_", Sys.Date(), "-", 
                  format(Sys.time(), "%H%M%S"), '.txt')
sink(sinkName)
cat("Single cell mode:", single_cell_mode, "\n")
if (single_cell_mode) {
    cat("Partition was based on both VDJ and VJ:", use_only_vdj, "\n")
}
sessionInfo()
sink()

vec_loci_vdj = c("TRB", "TRD")
vec_loci_vj = c("TRA", "TRG")

for (i in 1:nrow(subj_info)) {
    
    subj = subj_info[["subj"]][i]
    cat("\n************************\n")
    cat("\n", subj, "\n")

    # load db_vdj
    fn_in = subj_info[["path_vdj"]][i]
    
    if (grepl(pattern="\\.tsv$", x=fn_in)) {
        db_vdj = read.table(fn_in, header=T, sep="\t", stringsAsFactors=F)
        
    } else if (grepl(pattern="\\.RData$", x=fn_in)) {
        load(fn_in)
        db_vdj = db; rm(db)
    }
    
    if (single_cell_mode) {
        
        # load db_vj
        fn_in = subj_info[["path_vj"]][i]
        
        if (grepl(pattern="\\.tsv$", x=fn_in)) {
            db_vj = read.table(fn_in, header=T, sep="\t", stringsAsFactors=F)
            
        } else if (grepl(pattern="\\.RData$", x=fn_in)) {
            load(fn_in)
            db_vj = db; rm(db)
        }
        
        # combine VDJ and VJ
        stopifnot(all.equal(colnames(db_vdj), colnames(db_vj)))
        db_input = rbind(db_vdj, db_vj)
        
        rm(db_vdj, db_vj)
        
    } else {
        db_input = db_vdj
        rm(db_vdj)
    }
    
    lst_results = define_tcr_clone(db=db_input, 
                                   v_call=col_v_call, j_call=col_j_call, 
                                   junc_len=col_junc_len, junc=col_junc, 
                                   cell_id=col_cell_id, locus=col_locus, 
                                   single_cell_mode=single_cell_mode, 
                                   use_only_vdj=use_only_vdj, 
                                   clone_id=col_clone_id)
    
    db_passed = lst_results[["db_passed"]]
    db_failed = lst_results[["db_failed"]]
    
    bool_db_passed_vdj = db_passed[[col_locus]] %in% vec_loci_vdj
    bool_db_passed_vj = db_passed[[col_locus]] %in% vec_loci_vj
    
    if (any(bool_db_passed_vdj)) {
        db_clust_vdj = db_passed[bool_db_passed_vdj, ]
        fn_out_pass_base = paste0("cluster-pass_partition-", out_suffix, 
                                  "_chain-vdj_", subj)
        fn_out_pass_tsv = paste0(fn_out_pass_base, ".tsv")
        fn_out_pass_r = paste0(fn_out_pass_base, ".RData")
        write.table(x=db_clust_vdj, 
                    file=fn_out_pass_tsv, quote=F, sep="\t",
                    row.names=F, col.names=T)
        save(db_clust_vdj, file=fn_out_pass_r)
        rm(db_clust_vdj)
    }
    
    if (any(bool_db_passed_vj)) {
        db_clust_vj = db_passed[bool_db_passed_vj, ]
        fn_out_pass_base = paste0("cluster-pass_partition-", out_suffix, 
                                  "_chain-vj_", subj)
        fn_out_pass_tsv = paste0(fn_out_pass_base, ".tsv")
        fn_out_pass_r = paste0(fn_out_pass_base, ".RData")
        write.table(x=db_clust_vj, 
                    file=fn_out_pass_tsv, quote=F, sep="\t",
                    row.names=F, col.names=T)
        save(db_clust_vj, file=fn_out_pass_r)
        rm(db_clust_vj)
    }
    
    if (!is.null(db_failed)) {
        
        bool_db_failed_vdj = db_failed[[col_locus]] %in% vec_loci_vdj
        bool_db_failed_vj = db_failed[[col_locus]] %in% vec_loci_vj
        
        if (any(bool_db_failed_vdj)) {
            db_failed_vdj = db_failed[bool_db_failed_vdj, ]
            fn_out_fail_base = paste0("cluster-fail_partition-", out_suffix, 
                                      "_chain-vdj_", subj)
            fn_out_fail_tsv = paste0(fn_out_fail_base, ".tsv")
            fn_out_fail_r = paste0(fn_out_fail_base, ".RData")
            write.table(x=db_failed_vdj, 
                        file=fn_out_fail_tsv, quote=F, sep="\t",
                        row.names=F, col.names=T)
            save(db_failed_vdj, file=fn_out_fail_r)
            rm(db_failed_vdj)
        }
        
        if (any(bool_db_failed_vj)) {
            db_failed_vj = db_failed[bool_db_failed_vj, ]
            fn_out_fail_base = paste0("cluster-fail_partition-", out_suffix, 
                                      "_chain-vj_", subj)
            fn_out_fail_tsv = paste0(fn_out_fail_base, ".tsv")
            fn_out_fail_r = paste0(fn_out_fail_base, ".RData")
            write.table(x=db_failed_vj, 
                        file=fn_out_fail_tsv, quote=F, sep="\t",
                        row.names=F, col.names=T)
            save(db_failed_vj, file=fn_out_fail_r)
            rm(db_failed_vj)
        }
        
    }

}

# output subj_info_createGermlines.csv for CreateGermlines

if (opt$createGermlinesCSV) {
    
    all_subj = subj_info[["subj"]]
    
    cur_csv = cbind(subj=all_subj,
                    path_db_vdj=paste0(opt$pathWork, 
                                         "cluster-pass_partition-", out_suffix, 
                                         "_chain-vdj_", all_subj, ".tsv"))
    if (single_cell_mode) {
        cur_csv = cbind(cur_csv,
                        path_db_vj=paste0(opt$pathWork, 
                                             "cluster-pass_partition-", out_suffix, 
                                             "_chain-vj_", all_subj, ".tsv"))
    }
    cur_fn = "subj_info_createGermlines.csv"
    setwd(opt$pathAux)
    # note that headers are omitted for this file (col.names=F)
    write.table(x=cur_csv, file=cur_fn, quote=F, sep=",", row.names=F, col.names=F)
    
}

