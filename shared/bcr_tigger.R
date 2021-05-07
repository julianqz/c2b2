# copied from tigger/R/functions.R 
# last commit 2becdde on 2021-04-07; pulled 2021-05-07
# contains bug fix for `genotypeFasta` pushed in commit 9b023e2 on 2021-03-16
# bug fix not in stable release v1.0.0
genotypeFasta_fix = function(genotype, germline_db, novel=NA){
    if(!is.null(nrow(novel))){
        # Extract novel alleles if any and add them to germline_db
        novel <- filter(novel, !is.na(!!rlang::sym("polymorphism_call"))) %>%
            select(!!!rlang::syms(c("germline_call", "polymorphism_call", "novel_imgt")))
        if(nrow(novel) > 0){
            novel_gl <- novel$novel_imgt
            names(novel_gl) <- novel$polymorphism_call
            germline_db <- c(germline_db, novel_gl)
        }
    }
    
    genotype$gene <- getGene(genotype$gene, first = T, strip_d = T)
    g_names <- names(germline_db)
    names(g_names) <- getAllele(names(germline_db), first = T, strip_d = T)
    
    table_calls <- mapply(paste, genotype$gene, strsplit(genotype$alleles, ","),
                          sep="*")
    table_calls_names <- unlist(table_calls)
    seq_names <- g_names[names(g_names) %in% table_calls_names]
    seqs <- germline_db[seq_names]
    not_found <- !table_calls_names %in% names(g_names)
    
    if ( any(not_found) ) {
        stop("The following genotype alleles were not found in germline_db: ",
             paste(table_calls_names[not_found], collapse = ", "))
    }
    
    return(seqs)
}


#' Run tigger to infer genotype for an individual (without inferring novel alleles)
#'
#' @param   path_imgt              Path to fasta file containing IMGT germline reference.
#'                                 E.g. "./IGHV_no_dup.fasta".
#' @param   path_helper            Path to helper R script containing `read_multiline_fasta`.
#' @param   path_work              Path to save tigger outputs.
#' @param   subj                   Name of individual. Used in output filenames.
#' @param   db                     Input `data.frame`.
#' @param   col_seq                Column name of sequence column.
#' @param   col_v                  Column name of V call.
#' @param   col_prod               Column name of productive vs. non-productive.
#' @param   col_junc_len           Column name of junction length.
#' @param   p_fraction_to_explain  Parameter passed to `inferGenotype`.
#' @param   p_gene_cutoff          Parameter passed to `inferGenotype`.
#' @param   p_find_unmutated       Parameter passed to `inferGenotype`.
#' @param   p_text_size            Parameter passed to `plotGenotype`. 
#' 
#' @return  Various tigger objects in as .RData. Genotype plot.
#'          Updated db with `v_call_genotyped` column in .tsv and .RData formats.
#'          
#' @details Pay close attention to the genotype inferred and adjust 
#'          `p_find_unmutated` as necessary.
#'          
run_tigger = function(path_imgt, path_helper, path_work, 
                      subj, db, 
                      col_seq="sequence_alignment", 
                      col_v="v_call", 
                      col_prod="productive", 
                      col_junc_len="junction_length",
                      p_fraction_to_explain=0.875,
                      p_gene_cutoff=1e-04,
                      p_find_unmutated=T,
                      p_text_size=12) {
    
    # compatible with v1.0.0
    suppressPackageStartupMessages(require(tigger))
    suppressPackageStartupMessages(require(alakazam))
    
    ### read_multiline_fasta
    source(path_helper) 
    
    ### IMGT germline reference 
    germIMGT = toupper(read_multiline_fasta(path_imgt))
    # parse allele names
    names(germIMGT) = sapply(names(germIMGT), 
                             function(name){
                                 unlist(strsplit(x=name, split="\\|"))[2]
                             })
    
    cat("\n", subj, "\n")
    
    ### Sanity check
    # All sequences should be either heavy or light chain
    # relies on V call to perform this check
    # Assumes that within-sequence consistency btw V/D/J has been checked
    bool_all_heavy = !any(grepl(pattern="IG[KL]", x=toupper(db[[col_v]])))
    bool_all_light = !any(grepl(pattern="IGH", x=toupper(db[[col_v]])))
    stopifnot( bool_all_heavy+bool_all_light == 1 )
    
    chain_type = ifelse(bool_all_heavy, "heavy", "light")
    cat("\nchain_type:", chain_type, "\n")
    
    cat("\nnrow, initial:", nrow(db), "\n")
    
    cat("\nproductive vs. non-productive:\n")
    print(table(db[[col_prod]]))
    
    
    ### remove junction_length==0/NA (tigger would fail)
    cat("\nsummary on distribution of junction lengths:\n")
    print(summary(db[[col_junc_len]]))
    
    # non-productive could contribute junction_length of NA
    bool_jl_0 = db[[col_junc_len]]==0
    num_jl_0 = sum(bool_jl_0, na.rm=T)
    bool_jl_NA = is.na(db[[col_junc_len]])
    num_jl_NA = sum(bool_jl_NA)
    cat("\n# junction lengths being 0:", num_jl_0, "\n")
    cat("\n# junction lengths being NA:", num_jl_NA, "\n")
    
    # subset
    db = db[(!bool_jl_0 & !bool_jl_NA), ]
    cat("\nnrow, after removing junction_length of 0 and NA:", nrow(db), "\n")
    
    if (nrow(db)>0) {
        ### genotyping
        # Infer the individual's genotype, 
        # - using only unmutated sequences (default), and
        # - checking for the use of the novel alleles inferred (skipped here)
        
        novelDf = NA
        
        cat("\ninferGenotype()\n")
        # geno: a data.frame
        geno = inferGenotype(data=db, germline_db=germIMGT, novel=novelDf,
                             v_call=col_v, seq=col_seq,
                             fraction_to_explain=p_fraction_to_explain,
                             gene_cutoff=p_gene_cutoff,
                             find_unmutated=p_find_unmutated) 
        
        # genotype sequences to a vector
        cat("\ngenotypeFasta() \n")
        # genoSeqs: a named vector
        #* bug fixed in commit 9b023e2 but not in stable release yet
        genoSeqs = genotypeFasta_fix(genotype=geno, 
                                     germline_db=germIMGT, novel=novelDf)
        
        setwd(path_work)
        fn = paste0("geno_", subj, "_", chain_type, ".RData")
        save(geno, genoSeqs, file=fn)
        
        # visualize - bars indicate presence, not proportion.
        cat("\nplotGenotype()\n")
        fn = paste0("genotype_", subj, "_", chain_type, ".pdf")
        pdf(fn, width=6, height=6)
        plotGenotype(genotype=geno, gene_sort="name", text_size=p_text_size)
        dev.off()
        
        ### correct allele calls
        # Use the personlized genotype to determine corrected allele assignments
        # Currently only hamming and keep_gene=TRUE are implemented
        # Appends in $v_call_genotyped the best allele call from genotype_db
        # Bug-free if version >=v0.3.1
        
        cat("\nreassignAlleles()\n")
        db = reassignAlleles(data=db, genotype_db=genoSeqs, 
                             v_call=col_v, seq=col_seq, 
                             method="hamming", keep_gene="gene")
        
        fn = paste0("db_reassign_", subj, "_", chain_type, ".pdf")
        save(db, file=fn)
        
        ### split & export
        
        # productive
        tmp = db; rm(db)
        db = tmp[tmp[[col_prod]], ]
        cat("\n", chain_type, "productive, genotyped - nrow:",nrow(db), "\n")
        if (nrow(db)>0) {
            writeChangeoDb(data=db, file=paste0(subj, "_", chain_type, "_pr_genotyped.tsv"))
            save(db, file=paste0(subj, "_", chain_type, "_pr_genotyped.RData"))
            rm(db)
        }
        rm(db)
        
        # non-productive
        tmp = db; rm(db)
        db = tmp[!tmp[[col_prod]], ]
        cat("\n", chain_type, "non-productive, genotyped - nrow:",nrow(db), "\n")
        if (nrow(db)>0) {
            writeChangeoDb(data=db, file=paste0(subj, "_", chain_type, "_npr_genotyped.tsv"))
            save(db, file=paste0(subj, "_", chain_type, "_npr_genotyped.RData"))
            rm(db)
        }
        rm(db)
        
        rm(tmp)
        
        # prints blank if none
        warnings()
        
        cat("\nFinished for", subj, chain_type, "\n")
        
    } else {
        cat("\nNo sequence in db; skipped\n")
    }
    
}
