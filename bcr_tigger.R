# Julian Q. Zhou
# https://github.com/julianqz

#' Run tigger to infer genotype for an individual (without inferring novel alleles)
#'
#' @param   path_imgt              Path to fasta file containing IMGT germline reference.
#'                                 E.g. "./IGHV_no_dup.fasta".
#' @param   path_helper            Path to helper R script containing `read_multiline_fasta`.
#' @param   path_work              Path to save tigger outputs.
#' @param   subj                   Name of individual. Used in output filenames.
#' @param   db                     Input `data.frame`.
#' @param   no_split_version       If `TRUE`, export a no-split version combining
#'                                 productive and non-productive, in addition to
#'                                 the split version. Defaults to `FALSE`.
#' @param   col_seq                Column name of sequence column.
#' @param   col_v                  Column name of V call.
#' @param   col_prod               Column name of productive vs. non-productive.
#' @param   col_cdr3_len           Column name of cdr3 length.
#' @param   p_fraction_to_explain  Parameter passed to `inferGenotype`.
#' @param   p_gene_cutoff          Parameter passed to `inferGenotype`.
#' @param   p_find_unmutated       Parameter passed to `inferGenotype`.
#' @param   p_text_size            Parameter passed to `plotGenotype`.
#' @param   p_keep_gene            Parameter passed to `reassignAlleles`. 
#' 
#' @return  Various tigger objects in as .RData. Genotype plot.
#'          Updated db with `v_call_genotyped` column in .tsv and .RData formats.
#'          
#' @details Pay close attention to the genotype inferred and adjust 
#'          `p_find_unmutated` as necessary.
#'          
#'          Depending on the data, and settings of `p_fraction_to_explain` 
#'          and/or `p_keep_gene`, `v_call_genotyped` outputted by 
#'          `reassignAlleles` may or may not end up containing empty strings.
#'          When `p_keep_gene="gene"` and this happens, those `v_call_genotyped`
#'          values are padded with the corresponding `col_v` values and a message
#'          is printed.
#'          
#'          Chain type (heavy or light) is automatically identified and noted in
#'          the output filenames. Assumes that input `db` contains either all 
#'          heavy chains or all light chains, but not a mix of heavy and light.
#'          
run_tigger = function(path_imgt, path_helper, path_work, 
                      subj, db, no_split_version=F,
                      col_seq="sequence_alignment", 
                      col_v="v_call", 
                      col_prod="productive", 
                      col_cdr3_len="cdr3_length",
                      p_fraction_to_explain=0.875,
                      p_gene_cutoff=1e-04,
                      p_find_unmutated=T,
                      p_text_size=12,
                      p_keep_gene="gene") {
    
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
    
    
    ### remove cdr3_length==0/NA (tigger would fail)
    cat("\nsummary on distribution of", col_cdr3_len, ":\n")
    print(summary(db[[col_cdr3_len]]))
    
    # non-productive could contribute cdr3_length of NA
    bool_jl_0 = db[[col_cdr3_len]]==0
    num_jl_0 = sum(bool_jl_0, na.rm=T)
    bool_jl_NA = is.na(db[[col_cdr3_len]])
    num_jl_NA = sum(bool_jl_NA)
    cat("\n# cdr3 lengths being 0:", num_jl_0, "\n")
    cat("\n# cdr3 lengths being NA:", num_jl_NA, "\n")
    
    # subset
    db = db[(!bool_jl_0 & !bool_jl_NA), ]
    cat("\nnrow, after removing", col_cdr3_len, "of 0 and NA:", nrow(db), "\n")
    
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
        #* v1.0.0 stable release requires a bug fix from commit 9b023e2
        #* this bug fix is packed into julianqz/wu_cimm:main_0.1.1
        #* but if not using that docker image, v1.0.0 may fail next line
        genoSeqs = genotypeFasta(genotype=geno, 
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
                             method="hamming", keep_gene=p_keep_gene)
        
        fn = paste0("db_reassign_", subj, "_", chain_type, ".RData")
        save(db, file=fn)
        
        
        ## sanity check
        # there should be no empty string or NA in $v_call_genotyped
        col_new_v = "v_call_genotyped"
        
        bool_ck = db[[col_new_v]]!="" & !is.na(db[[col_new_v]]) 
        
        if (p_keep_gene=="gene" & !all(bool_ck)) {
            # pad
            idx_pad = which(!bool_ck)
            
            for (i_idx_pad in idx_pad) {
                cat("\npadding", col_new_v, "of row", i_idx_pad,
                    "with", col_v, ":", db[[col_v]][i_idx_pad], "\n")
            }
            
            db[[col_new_v]][idx_pad] = db[[col_v]][idx_pad]
            
            # recompute
            bool_ck = db[[col_new_v]]!="" & !is.na(db[[col_new_v]]) 
        }
        
        stopifnot(all(bool_ck))
        
        
        ### split & export
        
        # no-split version too?
        if (no_split_version) {
            cat("\n", chain_type, "productive + non-productive, genotyped - nrow:", nrow(db), "\n")
            if (nrow(db)>0) {
                #writeChangeoDb(data=db, file=paste0(subj, "_", chain_type, "_genotyped.tsv"))
                save(db, file=paste0(subj, "_", chain_type, "_genotyped.RData"))
            }
        }
        
        tmp = db; rm(db)
        
        # productive
        db = tmp[tmp[[col_prod]], ]
        cat("\n", chain_type, "productive, genotyped - nrow:",nrow(db), "\n")
        if (nrow(db)>0) {
            #writeChangeoDb(data=db, file=paste0(subj, "_", chain_type, "_pr_genotyped.tsv"))
            save(db, file=paste0(subj, "_", chain_type, "_pr_genotyped.RData"))
        }
        rm(db)
        
        # non-productive
        db = tmp[!tmp[[col_prod]], ]
        cat("\n", chain_type, "non-productive, genotyped - nrow:",nrow(db), "\n")
        if (nrow(db)>0) {
            #writeChangeoDb(data=db, file=paste0(subj, "_", chain_type, "_npr_genotyped.tsv"))
            save(db, file=paste0(subj, "_", chain_type, "_npr_genotyped.RData"))
        }
        rm(db)
        
        rm(tmp)
        
        # prints NULL if none
        cat("\nwarnings():\n")
        print(warnings())
        
        cat("\nFinished for", subj, chain_type, "\n")
        
    } else {
        cat("\nNo sequence in db; skipped\n")
    }
    
}
