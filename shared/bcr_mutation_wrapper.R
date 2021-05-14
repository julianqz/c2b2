#!/opt/conda/bin/Rscript

# wrapper to calculate mutation frequency

# assumes:
# - pathCSV points to a comma-separated file with the following headers
#   "subj", "path_db" where "path_db" points to a .RData file

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--pathHelper", action="store", default=NA, 
                type="character", help="Path to bcr_mutation_functions."),
    make_option("--pathCSV", action="store", default=NA, 
                type="character", help="Path to CSV containing subject list and paths to input files."),
    make_option("--pathWork", action="store", default=NA, type="character", 
                help="Path to working directory where outputs will be saved."),
    make_option("--nproc", action="store", default=1, type="numeric", 
                help="nproc for foreach."),
    make_option("--colObsv", action="store", default="sequence_alignment", 
                type="character", help="Column to observed sequences."),
    make_option("--colGerm", action="store", default="", 
                type="character", help="Column to germline sequences."),
    make_option("--useFull", action="store", default=TRUE, type="logical", 
                help="Whether to calculate using entire/full sequences. Boolean. If FALSE, segmentLimits must be specified."),
    make_option("--segmentLimits", action="store", default=NULL, type="character", 
                help="Comma-separated string parsed into lst_limits. Do not specify via command line if NULL. If NULL, useFull must be TRUE."),
    make_option("--colRemove", action="store", default=NULL, type="character", 
                help="Comma-separated string parsed into col_rmv. Do not specify via command line if NULL.")
)
opt = parse_args(OptionParser(option_list=option_list))

# without next line sink won't be able to capture shazam version
# (even if func in opt$pathHelper `require(shazam)`)
suppressPackageStartupMessages(require(alakazam))
suppressPackageStartupMessages(require(shazam))
suppressPackageStartupMessages(require(doParallel))
suppressPackageStartupMessages(require(foreach))

source(opt$pathHelper)

subj_info = read.table(opt$pathCSV, header=T, sep=",", stringsAsFactors=F)

nproc = opt$nproc
col_obsv = opt$colObsv
col_germ = opt$colGerm

# allowed combinations:
# - useFull TRUE,  segmentLimits NULL
# - useFull TRUE,  segmentLimits not NULL
# - useFull FALSE, segmentLimits not NULL
if (!opt$useFull & is.null(opt$segmentLimits)) {
    stop("At least one of these two must be true: {useFull is TRUE, segmentLimits is not NULL}.")
}

# parse
# \s is space
# ? means preceding item is optional and will be matched at most once

if (!is.null(opt$segmentLimits)) {
    # examples of expected input: 
    # "1,312"
    # "1,312;19,312"
    
    # first, separate different segments
    # the expected separator between segments is @
    # this works with just one segment without ; too
    vec_segments = strsplit(opt$segmentLimits, "\\s?@\\s?")[[1]]
    
    # next, within each segment, get the two limits
    # the expected separator between limits is comma
    # upper bound limit must be greater than lower limit
    lst_limits = sapply(vec_segments, function(s){
                            limits = as.integer(strsplit(s, "\\s?,\\s?")[[1]])
                            stopifnot(limits[2] > limits[1])
                            return(limits)
                        }, USE.NAMES=F, simplify=F)
    
    # use c(NA,NA) to signal "use entire/full sequences"
    if (opt$useFull) {
        lst_limits = c(lst_limits, list(c(NA,NA)))
    }
    
} else {
    if (opt$useFull) {
        lst_limits = list(c(NA,NA))
    }
}

stopifnot(length(lst_limits)>=1)


if (!is.null(opt$colRemove)) {
    col_rmv = strsplit(opt$colRemove, "\\s?;\\s?")[[1]]
} else {
    col_rmv = NULL
}


setwd(opt$pathWork)
sinkName = paste0("computingEnv_mutation_", Sys.Date(), "-", 
                  format(Sys.time(), "%H%M%S"), '.txt')
sink(sinkName)
# if NULL, will appear as "...: " (i.e. blank)
cat("Column to observed sequence:", col_obsv, "\n")
cat("Column to germline sequence:", col_germ, "\n")
cat("Compute using full sequence:", opt$useFull, "\n")
cat("Compute using specified segments:", opt$segmentLimits, "\n")
cat("Column to remove:", col_rmv, "\n")
cat("nproc:", nproc, "\n")
sessionInfo()
sink()


for (i in 1:nrow(subj_info)) {
    
    subj = subj_info[["subj"]][i]
    cat("\n", subj, "\n")
    
    # db
    load(subj_info[["path_db"]][i])
    
    # holder
    lst_add = vector(mode="list", length=length(lst_limits))
    
    for (j in 1:length(lst_limits)) {
        
        ### setup
        
        cur_limits = lst_limits[[j]]
        
        if (all(is.na(cur_limits))) {
            # use full sequence
            vec_input = db[[col_obsv]]
            vec_germ = db[[col_germ]]
            # suffix to append to colnames
            suffix = "_full"
        } else {
            stopifnot(!any(is.na(cur_limits)))
            # use substrings
            vec_input = substr(db[[col_obsv]], cur_limits[1], cur_limits[2])
            vec_germ = substr(db[[col_germ]], cur_limits[1], cur_limits[2])
            # suffix to append to colnames
            suffix = paste0("_", as.character(cur_limits[1]), 
                            "_", as.character(cur_limits[2]))
        }
        
        cat(" - suffix:", suffix, "\n")
        
        
        ### calculate per sequence
        
        if (nproc==1) {
            # If needed to run on a single core/cpu then, registerDoSEQ
            # Without doing this, foreach will give warning (though will still run)
            registerDoSEQ()
        } else if (nproc>1) {
            cluster = parallel::makeCluster(nproc, type="PSOCK")
            registerDoParallel(cluster)
            
            # export to cluster
            export_functions <- list("vec_input",
                                     "vec_germ",
                                     "calcObservedMutations",
                                     "parse_shazam_mutation",
                                     "calc_aa_mutation"
                                     )
            parallel::clusterExport(cluster, export_functions, envir=environment())
        }
    
        lst_calc = foreach(i_seq = 1:length(vec_input)) %dopar% {
            
            # nt, calculate
            # regionDefinition set to NULL as this has been handled via segmentLimits
            obj_nuc = calcObservedMutations(inputSeq=vec_input[i_seq],
                                            germlineSeq=vec_germ[i_seq],
                                            regionDefinition=NULL,
                                            returnRaw=T, frequency=F)
            
            # nt, parse into a named numeric vector
            obj_nuc_parsed = parse_shazam_mutation(obj_nuc)
            
            # aa, calculate
            # upper_bound set to NULL as this has been handled via segmentLimits
            # no need further parsing
            obj_aa = calc_aa_mutation(germ=vec_germ[i_seq], 
                                      obsv=vec_input[i_seq], 
                                      upper_bound=NULL)
            
            # return as a list
            lst_nuc_aa = list(obj_nuc_parsed=obj_nuc_parsed,
                              obj_aa=obj_aa)
            
            return(lst_nuc_aa)
        }
        
        # stop the cluster
        if (nproc>1) { parallel::stopCluster(cluster) }
        
        
        ### aggregate across sequences
        
        ## nt
        
        # extract
        lst_nuc_parsed = lapply(lst_calc, function(x){ x[["obj_nuc_parsed"]] })
        
        # e.g.
        #   nuc_denom nuc_R nuc_S nuc_RS nuc_RS_freq
        # 1       297     8     3     11 0.037037037
        db_nuc_parsed = data.frame(do.call(rbind, lst_nuc_parsed))
        
        # append calculation config to colnames
        colnames(db_nuc_parsed) = paste0(colnames(db_nuc_parsed), suffix)
        

        ## aa
        
        # extract
        lst_aa_parsed = lapply(lst_calc, function(x){ x[["obj_aa"]] })
        
        # aggregate across sequences
        db_aa_parsed = data.frame(do.call(rbind, lst_aa_parsed))
        # rename and rearrange to be consistent with db_nuc_parsed
        colnames(db_aa_parsed) = c("aa_R", "aa_R_freq", "aa_denom")
        db_aa_parsed = db_aa_parsed[, c("aa_denom", "aa_R", "aa_R_freq")]
        
        # append calculation config to colnames
        colnames(db_aa_parsed) = paste0(colnames(db_aa_parsed), suffix)
        
        
        ### record
        
        lst_add[[j]] = cbind(db_nuc_parsed, db_aa_parsed)
        
        rm(cur_limits, suffix, vec_input, vec_germ, lst_calc, cluster,
           lst_nuc_parsed, db_nuc_parsed, 
           lst_aa_parsed, db_aa_parsed)
    }
    
    
    # put everything together
    db = cbind(db, do.call(cbind, lst_add))
    
    # remove extra columns?
    if (!is.null(col_rmv)) {
        
        idx_rmv = match(col_rmv, colnames(db))
        
        if (any(!is.na(idx_rmv)) & !all(is.na(idx_rmv))) {
            idx_rmv = idx_rmv[!is.na(idx_rmv)]
            db = db[, -idx_rmv]
        }
    }
    
    # save
    fn = paste0("mu_freq_", subj, ".RData")
    save(db, file=fn)
    
    rm(db, fn, lst_add)
}

