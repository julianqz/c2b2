#' Read in a multi-line fasta file
#' 
#' Given a fasta file in which each entry is represented by a fasta header that
#' with a `>` and one or more lines of the actual sequence (as opposed to always
#' by a fasta header followed by a single line of sequence), create a vector
#' in which each item corresponds to a fasta entry.
#' 
#' @param   filename         filename of multi-line fasta file.
#' 
#' @return  Returns a character vector. Each vector item corresponds to a fasta
#'          entry. The names of the items correspond to the fasta headers.
#'          
#' @details This function works with fasta files downloaded from IMGT Ig germline
#'          references, in which the headers are very long with multiple `|`s. 
#'          In contrast, `seqinr::read.fasta` does not read in the entire header
#'          in those cases. 
#'          
#'          This function is renamed from `readIMGTfasta` originally in 
#'          `spatial/getGermIMGT.R`.         

read_multiline_fasta = function(filename) {
    suppressPackageStartupMessages(require(stringi))
    rawLines = readLines(filename)
    headerIdx = which(grepl(pattern=">", x=rawLines))
    nSeqs = length(headerIdx)
    fasta = vector("character", length=nSeqs)
    for (i in 1:nSeqs) {
        
        idxFirst = headerIdx[i]+1
        
        # first through second last entry
        if (i<nSeqs) {
            idxLast = headerIdx[i+1]-1
        } else {
            # last entry
            idxLast = length(rawLines)
        }
        
        # if last entry is followed by a "" on the very last line in file, 
        # that "" goes away automatically when running next line
        fasta[i] = paste0(rawLines[idxFirst:idxLast], collapse="")
    }
    
    names(fasta) = rawLines[headerIdx]
    return(fasta)
}


#' Export sequences as a fasta file
#' 
#' Given a vector of sequences and their headers, write the sequences as a fasta
#' file.
#' 
#' @param   sequences          a vector of sequences
#' @param   headers            a vector of headers for the sequences
#' @param   add_header_symbol  a Boolean value indicating whether the `>` symbol
#'                             is to be added as prefix to the headers
#' @param   filename           filename (include ".fasta") to be written
#' 
#' @return  Writes a fasta file to the working directory. 
#' 
#' @details A double-line fasta is written, as opposed to a multi-line fasta

export_fasta = function(sequences, headers, add_header_symbol, filename) {
    
    stopifnot(length(sequences)==length(headers))
    
    sink(filename)
    for (i in 1:length(sequences)) {
        # header
        cat(ifelse(add_header_symbol, ">", ""), headers[i], sep="", "\n")
        # sequence
        cat(sequences[i], sep="", "\n")
    }
    sink()
}

#' Remove duplicate fasta entries
#' 
#' Given a fasta file, remove duplicate entries based on sequence identity.
#' 
#' @param   filename       filename of fasta file. Multi-line fasta file supported.
#' 
#' @return  A named character vector containing unique fasta entries. Fasta 
#'          headers are in the names of the vector.
#'          
#' @details For each unique sequence corresponding to multiple duplicate fasta
#'          entries, only one fasta entry is kept. The entry to be kept is the 
#'          one that appears the first amongst the duplicate entries in the 
#'          input fasta file.          

remove_duplicate_fasta = function(filename) {
    cat("\n", filename, "\n")
    
    # read in fasta entries
    vec = read_multiline_fasta(filename)
    
    # are there duplicate entries in terms of sequence identity?
    bool_dup = length(unique(vec))<length(vec)
    
    # if there are duplicate entries
    if (bool_dup) {
        # tabulate sequences
        tab = table(vec)
        
        # sequences of duplicate entries
        seqs_dup = names(tab[tab>1])
        
        # Each list entry corresponds to a single sequence
        # For each sequence, make a table showing
        # index wrt vec, and header, of duplicate entries
        lst_dup = sapply(1:length(seqs_dup),
                         function(i) {
                             s_idx = which(vec==seqs_dup[i])
                             s_headers = names(vec)[s_idx]
                             
                             df = data.frame(matrix(NA, nrow=length(s_idx), ncol=3))
                             colnames(df) = c("idx_dup", "idx_vec", "header")
                             df[["idx_dup"]] = rep(i, length(s_idx))
                             df[["idx_vec"]] = s_idx
                             df[["header"]] = s_headers
                             
                             return(df)
                         }, simplify=F, USE.NAMES=F)
        
        # compile tables into a single data.frame
        df_dup = do.call(rbind, lst_dup)
        
        # For each duplicate seq, keep only one fasta entry
        # Here, the first entry (as it appears in `vec`) is kept
        
        # idx wrt vec of fasta entries to be kept
        idx_keep = rep(NA, length(seqs_dup))
        
        for (i in 1:length(seqs_dup)) {
            idx_df_dup = which(df_dup[["idx_dup"]]==i)
            
            # first entry is kept
            idx_keep[i] = df_dup[["idx_vec"]][idx_df_dup[1]]
            
            # verbose
            cat("---------------", i, "---------------\nIn:", 
                df_dup[["header"]][idx_df_dup[1]],
                ";\nOut:", df_dup[["header"]][idx_df_dup[-1]], "\n")
        }
        # idx_keep should all have been filled with non-NA
        stopifnot(!any(is.na(idx_keep)))
        
        # idx wrt vec of fasta entries to be excluded
        idx_throw = df_dup[["idx_vec"]][ ! df_dup[["idx_vec"]] %in% idx_keep ]
        
        # idx_keep and idx_throw should have no overlap
        stopifnot( length(intersect(idx_keep, idx_throw)) == 0 )
        # union of idx_keep and idx_throw should match df_dup[["idx_vec]]
        stopifnot( length(setdiff(c(idx_keep, idx_throw), df_dup[["idx_vec"]])) == 0 )
        
        # remove duplicate fasta entries
        vec = vec[-idx_throw]
        
        # all fasta entries should now be unique
        stopifnot( length(unique(vec)) == length(vec) )
        # number of unique sequences should be the same before & after
        stopifnot( length(vec) == length(tab) )
        
    } else {
        cat("No duplicate fasta entries found.\n")
    }
    
    return(vec)
}


#' Remove IMGT gaps in germline sequence and optionally also 
#' in observed sequence based on gaps present in germline sequence.
#' 
#' Adapted from helpers.R for JI 2020
#' 
#' @param    germ  IMGT-gapped germline sequence.
#' @param    obsv  IMGT-gapped observed sequence(s); optional.
#' 
#' @return   A vector containing `germ_no_gaps` and `obsv_no_gaps`.
#' 
#' @details  IMGT gaps are removed from germline. 
#'           Corresponding gap positions in observed are also removed.
#'        
#' @examples
#' germ = "ABC...EDF..GTGH...JHK.....LOP....W"
#' obsv = "ABCXXXEDFXXGTGHXXXJHKXXXXXLOPXXXXW"
#' remove_imgt_gaps(germ, NULL)
#' remove_imgt_gaps(obsv, NULL)
#' remove_imgt_gaps(germ, obsv)
#' 
remove_imgt_gaps = function(germ, obsv=NULL) {
    
    require(stringi)
    
    # only IMGT gaps (triple dots "..." are removed)
    if (stri_detect(str=germ, regex="\\.\\.\\.")) {
        # remove IMGT gaps from germline
        germ_no_gaps = stri_replace(str=germ, replacement="", 
                                    regex="\\.\\.\\.", mode="all")
        
        if (!is.null(obsv)) {
            # locate IMGT gaps in germline
            # matrix: rows: instances; cols: start & end
            gap_pos = stri_locate(str=germ, regex="\\.\\.\\.", mode="all")[[1]]
            # remove IMGT gaps from observed sequence
            for (i in 1:nrow(gap_pos)) {
                stri_sub(str=obsv, from=gap_pos[i, "start"], to=gap_pos[i, "end"]) <- ""
                gap_pos = gap_pos - 3
            }
        }
        obsv_no_gaps = obsv
    } else {
        germ_no_gaps = germ
        obsv_no_gaps = obsv
    }
    
    return(c(germ_no_gaps=germ_no_gaps, obsv_no_gaps=obsv_no_gaps))
}

#' Clean an observed sequence in preparation for expression 
#' 
#' @param   obsv   IMGT-gapped observed sequence.
#' @param   germ   IMGT-gapped germline sequence.
#' @param   trim   Boolean. Whether to trim output nt length to a multiple of 3.
#' 
#' @return  A vector containing `germ_clean` and `obsv_clean`.
#' 
#' @details `obsv` undergoes the following cleaning steps:
#'          1) IMGT gaps are removed.
#'          2) If there's any non-ATGC position, patch those positions with
#'             the corresponding germline positions, provided that the 
#'             germline positions are non-ATGC.
#'          3) If `trim=TRUE`, trim nt length to a multiple of 3.
#' 
clean_obsv = function(obsv, germ, trim) {
    require(seqinr)
    
    stopifnot(nchar(obsv)==nchar(germ))
    
    # remove IMGT gaps
    vec_no_gaps = remove_imgt_gaps(obsv=obsv, germ=germ)
    obsv_no_gaps = vec_no_gaps["obsv_no_gaps"]
    germ_no_gaps = vec_no_gaps["germ_no_gaps"]
    
    # if there's any non-ATGC in obsv, 
    # patch with corresponding positions from germ
    vec_atgc = c("A","T","G","C","a","t","g","c")
    if (grepl(pattern="[^ATGCatgc]", x=obsv_no_gaps)) {
        # convert to vector of single characters
        obsv_no_gaps_c = s2c(obsv_no_gaps)
        germ_no_gaps_c = s2c(germ_no_gaps)
        
        # which positions in obsv is non-ATGC
        idx_obsv_nonATGC = which(!obsv_no_gaps_c %in% vec_atgc)
        stopifnot(length(idx_obsv_nonATGC)>0)
        
        # are the corresponding positions ATGC in germ?
        bool_germ = germ_no_gaps_c[idx_obsv_nonATGC] %in% vec_atgc
        
        if (any(bool_germ)) {
            # at least one position can be patched by germline
            
            if (any(!bool_germ)) {
                # one or more of the corresponding positions in germ is non-ATGC
                cat("germline positions (", idx_obsv_nonATGC[!bool_germ], 
                    ") non-ATGC (", germ_no_gaps_c[idx_obsv_nonATGC[!bool_germ]],
                    "); patching skipped for these positions\n")
            }
            
            # exclude any non-ATGC germline positions for patching
            idx_obsv_nonATGC_patch = idx_obsv_nonATGC[bool_germ]
            
            # patch
            cat("patching performed for (", obsv_no_gaps_c[idx_obsv_nonATGC_patch],
                ") at positions (", idx_obsv_nonATGC_patch, 
                ") in obsv with germline (",
                germ_no_gaps_c[idx_obsv_nonATGC_patch], ")\n")
            obsv_no_gaps_c[idx_obsv_nonATGC_patch] = germ_no_gaps_c[idx_obsv_nonATGC_patch]
            
        } else {
            # none of the corresponding positions in germ is ATGC
            cat("all germline positions of concern (", idx_obsv_nonATGC, 
                ") non-ATGC (", germ_no_gaps_c[idx_obsv_nonATGC],
                "); patching skipped entirely\n")
        }
        
        # sanity check
        # the number of non-ATGC position left in obsv after patching, if any, 
        # should match the number of non-ATGC in germline
        stopifnot( sum(!obsv_no_gaps_c %in% vec_atgc) == sum(!bool_germ) )
        
        # convert back to string
        obsv_no_gaps = c2s(obsv_no_gaps_c)
    }
    
    # if trim=TRUE, trim length to multiple of 3
    if (trim) {
        stopifnot(nchar(obsv_no_gaps)==nchar(germ_no_gaps))
        len_obsv_no_gaps = nchar(obsv_no_gaps)
        if (len_obsv_no_gaps %% 3 != 0) {
            new_len = len_obsv_no_gaps - len_obsv_no_gaps%%3
            stopifnot(new_len %% 3 == 0)
            obsv_no_gaps = substr(obsv_no_gaps, 1, new_len)
            germ_no_gaps = substr(germ_no_gaps, 1, new_len)
        }
    }
    
    vec_return = c(obsv_no_gaps, germ_no_gaps)
    names(vec_return) = c("obsv_clean", "germ_clean")
    return(vec_return)
}

