# Julian Q. Zhou
# https://github.com/julianqz

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
#' @param   vec_fasta      a named character vector; names correspond to fasta IDs
#' 
#' @return  A named character vector containing unique fasta entries. Fasta 
#'          headers are in the names of the vector.
#'          
#' @details For each unique sequence corresponding to multiple duplicate fasta
#'          entries, only one fasta entry is kept. The entry to be kept is the 
#'          one that appears the first amongst the duplicate entries in the 
#'          input fasta file or vector.
#'          
#'          Exactly one of `filename` or `vec_fasta` should be specified, with
#'          the other one being left as `NULL`.
#'                    
remove_duplicate_fasta = function(filename=NULL, vec_fasta=NULL) {
    
    # check input
    # exactly one of file and vec_fasta must be NULL
    stopifnot(sum(c(is.null(filename), is.null(vec_fasta)))==1)
    
    if (!is.null(filename)) {
        cat("\n", filename, "\n")
        
        # read in fasta entries
        vec = read_multiline_fasta(filename)
    } else {
        vec = vec_fasta
    }

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
            cat("---------------", i, "---------------\n In:", 
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
#' @param   vdj_obsv        IMGT-gapped observed sequence.
#' @param   vdj_germ        IMGT-gapped germline sequence.
#' @param   full_seq        Input sequence. The `sequence` column in AIRR format.
#' @param   vdj_obsv_start  Position in `full_seq` where `vdj_obsv` starts.
#'                          The `v_sequence_start` column in AIRR format. 
#'                          See details for finer points.
#' @param   vdj_obsv_end    Position in `full_seq` where `vdj_obsv` ends.
#'                          The `j_sequence_end` column in AIRR format.
#'                          See details for finer points.                          
#' 
#' @return  A vector containing 
#'          - `germ_clean`: germline sequence with IMGT gaps removed. 
#'          - `obsv_clean`: "cleaned" observed sequence with IMGT gaps removed,
#'                          non-ATGC positions corrected by germline, and nt length
#'                          trimmed to a multiple of 3.
#'          - `spacer_5`:   a 9-bp extracted upstream of IMGT-numbered nt pos 1.
#'          - `obsv_end`:   the position wrt `full_seq` of the last nt pos in 
#'                          `obsv_clean`.
#'          - `correction`: notes where germline-guided correction was performed.
#'          - `review`:     notes where manual review may be required.
#'          
#' @details `obsv` undergoes the following cleaning steps:
#' 
#'          1) IMGT gaps are removed. 
#'             IMGT gaps are identified as positions that contain in-frame 
#'             triplet `...` in `vdj_germ`.
#'             
#'          2) If there's any non-ATGC position, patch those positions with
#'             the corresponding germline positions, provided that the 
#'             germline positions are non-ATGC.
#'             
#'             If any of the positions that need patching is non-ATGC in the 
#'             germline, patching is skipped for that position. A note is 
#'             added to `review` with the format 
#'             `germline_non-ATGC=[position]=[char]`.
#'             
#'             Any position that receives patching is noted in `correction` 
#'             with the format 
#'             `[char before patching]|position|[char after patching]`.
#'             
#'             Presence of `-` before patching indicates potential deletion. 
#'             A note is added to `review`.
#'             
#'          3) Trim nt length to a multiple of 3.
#'             If trimming is performed, a note is added to `correction` to 
#'             record the position(s) trimmed. 
#'          
#'          In addition, a 9-bp 5' spacer is extracted from `full_seq`.
#'          A note is added to `review` if there's insufficient positions in 
#'          `full_seq` for extracting the spacer, or if there's any non-ATGC
#'          position in the spacer. A message is also printed.
#'          
#'          Prep work for extracting a 3' spacer that extends from the 3' end of
#'          the junction onwards is also carried out. Actual 
#'          extraction of the 3' spacer is to be performed by a separate function,
#'          after the results returned by this function has been reviewed. 
#'          This is because the review could result in changes in `obsv_clean`,
#'          which would affect the part of the 3' spacer that lies within the 
#'          VDJ.
#'            
#'          If the IMGT-aligned observed sequence contains `.` or `N` at the 
#'          very beginning or the very end, `v_sequence_start` and 
#'          `j_sequence_end` do not include those positions. Instead, they 
#'          start after and end before those positions. These positions are
#'          noted in `reivew` as `lead=` and `trail=` respectively.
#'          This is taken into account when extracting the spacers. 
#'          
#'          The 5' spacer is extracted 9-bp upstream of IMGT-numbered nt position 1. 
#'          
#'          The part of the 3' spacer that lies outside the VDJ (if any) is 
#'          to be extracted downstream of `obsv_end`, which is the position in
#'          `full_seq` of the last nt, post-trimming, of `obsv_clean`.
#'          
clean_obsv = function(vdj_obsv, vdj_germ,
                      full_seq, vdj_obsv_start, vdj_obsv_end) {
    require(seqinr)
    require(stringi)
    
    stopifnot(nchar(vdj_obsv)==nchar(vdj_germ))
    
    
    review = ""
    
    
    # 1) remove IMGT gaps
    vec_no_gaps = remove_imgt_gaps(obsv=vdj_obsv, germ=vdj_germ)
    vdj_obsv_no_gaps = vec_no_gaps["obsv_no_gaps"]
    vdj_germ_no_gaps = vec_no_gaps["germ_no_gaps"]
    len_vdj_obsv_no_gaps = nchar(vdj_obsv_no_gaps)
    
    # check if there's leading N or non-IMGT "."
    bool_lead = grepl(pattern="^[Nn\\.]", x=vdj_obsv_no_gaps)
    if (bool_lead) {
        # count number of leading N or non-IMGT "."
        lead_dot = stri_extract_first(str=vdj_obsv_no_gaps, regex="^[Nn\\.]+")
        stopifnot(all( s2c(lead_dot) %in% c(".", "N", "n") ))
        n_lead = nchar(lead_dot)
        review = paste(review, paste0("lead=", lead_dot),
                       sep=ifelse(review=="", "", ";"))
    } else {
        n_lead = 0
    }
    
    # check if there's trailing non-IMGT "."
    bool_trail = grepl(pattern="[Nn\\.]$", x=vdj_obsv_no_gaps)
    if (bool_trail) {
        # count number of trailing non-IMGT "."
        trail_dot = stri_extract_last(str=vdj_obsv_no_gaps, regex="[Nn\\.]+$")
        stopifnot(all( s2c(trail_dot) %in% c(".", "N", "n") ))
        n_trail = nchar(trail_dot)
        review = paste(review, paste0("trail=", trail_dot),
                       sep=ifelse(review=="", "", ";"))
    } else {
        n_trail = 0
    }
    
    # the length without IMGT gaps should usually match vdj_obsv_end-start+1
    # push further upstream 5' if there're leading dots
    # push further downstream 3' if there're trailing dots
    # exceptions:
    # - insertions (excluded in IMGT-aligned seq)
    #   (eg. 368-22, h, i=23, 159; l, i=22)
    len_ck = (vdj_obsv_end+n_trail)-(vdj_obsv_start-n_lead)+1
    bool_len_ck = len_vdj_obsv_no_gaps == len_ck
    if (!bool_len_ck) {
        cat("len_vdj_obsv_no_gaps =", len_vdj_obsv_no_gaps, "vs len_ck =", len_ck, "\n")
        review = paste(review, 
                       paste0("vdj_length_potential_indel=", 
                              as.character(len_ck-len_vdj_obsv_no_gaps)),
                       sep=ifelse(review=="", "", ";"))
    }
    
    
    # 2) if there's any non-ATGC in obsv, 
    # patch with corresponding positions from germ
    vec_atgc = c("A","T","G","C","a","t","g","c")
    # [^] means any character NOT in the list
    if (grepl(pattern="[^ATGCatgc]", x=vdj_obsv_no_gaps)) {
        # convert to vector of single characters
        vdj_obsv_no_gaps_c = s2c(vdj_obsv_no_gaps)
        vdj_germ_no_gaps_c = s2c(vdj_germ_no_gaps)
        
        # which positions in obsv is non-ATGC
        idx_vdj_obsv_nonATGC = which(!vdj_obsv_no_gaps_c %in% vec_atgc)
        stopifnot(length(idx_vdj_obsv_nonATGC)>0)
        
        # are the corresponding positions ATGC in germ?
        bool_germ = vdj_germ_no_gaps_c[idx_vdj_obsv_nonATGC] %in% vec_atgc
        
        if (any(bool_germ)) {
            # at least one position can be patched by germline
            
            if (any(!bool_germ)) {
                # one or more of the corresponding positions in germ is non-ATGC
                cat("germline positions (", idx_vdj_obsv_nonATGC[!bool_germ], 
                    ") non-ATGC (", vdj_germ_no_gaps_c[idx_vdj_obsv_nonATGC[!bool_germ]],
                    "); patching skipped for these positions\n")
                review = paste(review, 
                               paste0("germline_non-ATGC=", 
                                      paste(as.character(idx_vdj_obsv_nonATGC[!bool_germ]), 
                                            collapse=","),
                                      "=", 
                                      paste(vdj_germ_no_gaps_c[idx_vdj_obsv_nonATGC[!bool_germ]], 
                                            collapse=",")
                                      ),
                               sep=ifelse(review=="", "", ";"))
            }
            
            # exclude any non-ATGC germline positions for patching
            idx_vdj_obsv_nonATGC_patch = idx_vdj_obsv_nonATGC[bool_germ]
            
            # patch
            bf_patch = vdj_obsv_no_gaps_c[idx_vdj_obsv_nonATGC_patch]
            af_patch = vdj_germ_no_gaps_c[idx_vdj_obsv_nonATGC_patch]
            
            cat("patching performed for (", bf_patch,
                ") at positions (", idx_vdj_obsv_nonATGC_patch, 
                ") in obsv with germline (", af_patch, ")\n")
            
            correction = paste0(paste(bf_patch, collapse=""), "|",
                                paste(as.character(idx_vdj_obsv_nonATGC_patch),
                                      collapse=","), "|",
                                paste(af_patch, collapse=""))
                
            if (any(bf_patch=="-")) {
                review = paste(review, "observed_potentital_del", 
                               sep=ifelse(review=="", "", ";"))
            }
            
            vdj_obsv_no_gaps_c[idx_vdj_obsv_nonATGC_patch] = af_patch
            
        } else {
            # none of the corresponding positions in germ is ATGC
            cat("all germline positions of concern (", idx_vdj_obsv_nonATGC, 
                ") non-ATGC (", vdj_germ_no_gaps_c[idx_vdj_obsv_nonATGC],
                "); patching skipped entirely\n")
            
            correction = ""
            
            review = paste(review, 
                           paste0("germline_non-ATGC=", 
                                  paste(as.character(idx_vdj_obsv_nonATGC), 
                                        collapse=","),
                                  "=", 
                                  paste(vdj_germ_no_gaps_c[idx_vdj_obsv_nonATGC], 
                                        collapse=",")
                           ),
                           sep=ifelse(review=="", "", ";"))
        }
        
        # sanity check
        # the number of non-ATGC position left in obsv after patching, if any, 
        # should match the number of non-ATGC in germline
        stopifnot( sum(!vdj_obsv_no_gaps_c %in% vec_atgc) == sum(!bool_germ) )
        
        # convert back to string
        vdj_obsv_no_gaps = c2s(vdj_obsv_no_gaps_c)
        
        # sanity check
        # vdj_obsv len should stay unchanged
        stopifnot(nchar(vdj_obsv_no_gaps)==len_vdj_obsv_no_gaps)
    } else {
        correction = ""
    }
    
    # 3) # trim if length is not multiple of 3
    mod3 = len_vdj_obsv_no_gaps %% 3
    if (mod3!=0) {
        new_len = len_vdj_obsv_no_gaps - mod3
        stopifnot(new_len %% 3 == 0)
        
        vdj_obsv_no_gaps = substr(vdj_obsv_no_gaps, 1, new_len)
        vdj_germ_no_gaps = substr(vdj_germ_no_gaps, 1, new_len)
        
        correction = paste(correction, 
                           paste0("trim=",
                                  paste(as.character( (new_len+1):len_vdj_obsv_no_gaps ), 
                                        collapse=",")),
                           sep=ifelse(correction=="", "", ";"))
        
        # update len_vdj_obsv_no_gaps
        len_vdj_obsv_no_gaps = nchar(vdj_obsv_no_gaps)
        stopifnot(len_vdj_obsv_no_gaps==new_len)
    }
    
    
    # 4) extract 5' spacer and prep for extracting 3' spacer
    
    # modify vdj start and end positions for extracting spacers
    # wrt full_seq
    
    # e.g.
    # 122222222 pos in sequence
    # 901234567 19,20,...,27
    # ....tcacc sequence_alignment; vdj_start 23
    # xxxxtcacc padded; vdj_start_mod 23-4=19
    vdj_obsv_start_mod = vdj_obsv_start-n_lead
    
    # e.g.
    # 122222222 pos in sequence
    # 901234567 19,20,...,27
    # atccc...  sequence_alignment; vdj_end 23
    # atcccxxx  padded
    # atcccx    padded & trimmed; vdj_end_mod 23+3-2=24
    vdj_obsv_end_mod = vdj_obsv_end+n_trail-mod3
    
    #* numbers specified by WK
    spacer_5_len = 9
    
    # 5'
    # upstream of IMGT-numbered pos 1
    spacer_5 = substr(full_seq, 
                      start=(vdj_obsv_start_mod-spacer_5_len), 
                      stop=(vdj_obsv_start_mod-1))
    if (nchar(spacer_5)!=spacer_5_len) {
        cat("nchar(spacer_5) =", nchar(spacer_5), "; expected ", spacer_5_len, "\n")
        review = paste(review, 
                       paste0("spacer_5_length=", nchar(spacer_5)), 
                       sep=ifelse(review=="", "", ";"))
    }
    
    
    # sanity check
    # no non-ATGC in fina obsv
    stopifnot(!grepl(pattern="[^ATGCatgc]", x=vdj_obsv_no_gaps))
    # print message if there's non-ATGC in spacers
    if (grepl(pattern="[^ATGCatgc]", x=spacer_5)) {
        cat("non-ATGC in 5' spacer:", spacer_5, "\n")
        review = paste(review, "spacer_5_non-ATGC", 
                       sep=ifelse(review=="", "", ";"))
    }
    
    # return
    vec_return = c(vdj_obsv_no_gaps, vdj_germ_no_gaps, 
                   spacer_5, vdj_obsv_end_mod,
                   correction, review)
    names(vec_return) = c("obsv_clean", "germ_clean",
                          "spacer_5", "obsv_end",
                          "correction", "review")
    return(vec_return)
}

# for testing clean_obsv
RUN=F
if (RUN) {
    i=23
    vdj_obsv = db_selected[["h_sequence_alignment"]][i]
    vdj_germ = db_selected[["h_germline_alignment"]][i]
    locus=db_selected[["h_locus"]][i]
    full_seq=db_selected[["h_sequence"]][i]
    vdj_obsv_start=db_selected[["h_v_sequence_start"]][i]
    vdj_obsv_end=db_selected[["h_j_sequence_end"]][i]
    
    i=22
    vdj_obsv = db_selected[["l_sequence_alignment"]][i]
    vdj_germ = db_selected[["l_germline_alignment"]][i]
    locus=db_selected[["l_locus"]][i]
    full_seq=db_selected[["l_sequence"]][i]
    vdj_obsv_start=db_selected[["l_v_sequence_start"]][i]
    vdj_obsv_end=db_selected[["l_j_sequence_end"]][i]
}

# for testing extract_3_primer_spacer
RUN=F
if (RUN) {
    i=2
    vdj_obsv = db_selected[["l_obsv"]][i]
    vdj_obsv_end = db_selected[["l_obsv_end"]][i]
    full_seq = db_selected[["l_sequence"]][i]
    locus = db_selected[["l_locus"]][i]
    cdr3 = db_selected[["l_cdr3"]][i]
}

#' Extract a 3' spacer for expression of Ab
#' 
#' @param  vdj_obsv      "Cleaned" observed VDJ sequence. The expectation is that 
#'                       this was returned as `obsv_clean` by `clean_obsv`.
#' @param  vdj_obsv_end  The position wrt `full_seq` of the last nt pos in `vdj_obsv`.
#'                       The expectation is that this was returned as `obsv_end` 
#'                       by `clean_obsv`.
#' @param  full_seq      Input sequence. The `sequence` column in AIRR format.
#' @param  locus         Locus. One of `IGH`, `IGK`, or `IGL`.
#' @param  cdr3          Nt sequence of CDR3 (note: NOT IMGT-defined junction).
#' 
#' @return  A vector containing
#'          - spacer_start: starting position wrt `vdj_obsv` of the 3' spacer.
#'                          Namely, the 1st nt pos downstream of the 3' end of
#'                          the IMGT-defined junction (the 4th nt pos downstream 
#'                          of the 3' end of the CDR3). 
#'          - pt_out: if any, the nt sequence of the part of the 3' spacer that 
#'                    lies outside the VDJ.        
#'          - review: notes where manual review may be required.
#' 
#' @details  The total length of the 3' spacer is specified as 30/27/66-bp for
#'           IGH/IGK/IGL respectively, downstream of the 3' end of the IMGT-defined
#'           junction.
#'           
#'           If part of the 3' spacer lies outside the VDJ, that part is extracted
#'           and returned. If the length of `full_seq` outside the VDJ is not 
#'           sufficient, a note is added to `review` with the format 
#'           `spacer_3_outside_vdj_length=[actual len]/[expected len]`. 
#'           A message is also printed.
#'           
#'           If the 3' spacer lies entirely inside the VDJ but actually ends 
#'           before the end of the VDJ, a note is added to `review` with the 
#'           format `spacer_3_inside_vdj_length=[len]`.
#'           
#'           If there's any non-ATGC character in the part of the 3' spacer 
#'           that lies outside the VDJ, a note is added to `review` with the 
#'           format `spacer_3_outside_vdj_non-ATGC=[char]`. A message 
#'           is also printed.
#'           
extract_3_prime_spacer = function(vdj_obsv, vdj_obsv_end, full_seq, 
                                   locus, cdr3) {
    require(stringi)
    
    stopifnot(locus %in% c("IGH","IGK","IGL"))
    
    total_len_vec = c("IGH"=30, "IGK"=27, "IGL"=66)
    
    total_len = total_len_vec[locus]
    
    review = ""
    
    # locate cdr3
    # returns a matrix
    # no match:
    #      start end
    # [1,]    NA  NA
    # match:
    #      start end
    # [1,]   289 330
    loc_cdr3 = stri_locate(str=vdj_obsv, fixed=cdr3)
    # expect match
    stopifnot(!any(is.na(loc_cdr3)))
    # expect exactly 1 match
    stopifnot(nrow(loc_cdr3)==1)
    
    # part from within obsv
    # +4 to skip anchor AA (downstream of junction)
    spacer_start = loc_cdr3[1, 2]+4
    pt_obsv = substring(vdj_obsv, first=spacer_start)
    len_pt_obsv = nchar(pt_obsv)
    
    if (len_pt_obsv==total_len) {
        # no need for spacer downstream vdj_obsv_end
        # just enough
        cat("no need for spacer downstream vdj_obsv_end\n")
        pt_out = ""
        
    } else if (len_pt_obsv > total_len) {
        # no need for spacer downstream vdj_obsv_end
        # in fact, vdj_obsv seems too long -- problematic; make a note
        cat("no need for spacer downstream vdj_obsv_end; part within vdj too long\n")
        pt_out = ""
        
        review = paste(review, 
                       paste0("spacer_3_inside_vdj_length=", len_pt_obsv),
                       sep=ifelse(review=="", "", ";"))
        
    } else if (len_pt_obsv < total_len) {
        
        # expected length of spacer downstream of vdj_obsv_end
        pt_out_len_expected = total_len - len_pt_obsv
        # extract
        pt_out = substr(full_seq, vdj_obsv_end+1, vdj_obsv_end+pt_out_len_expected)
        # actual length
        pt_out_len_actual = nchar(pt_out)
        
        if (pt_out_len_expected != pt_out_len_actual) {
            
            # eg. not enough positions in full_seq for spacer
            cat("nchar(pt_out) =", pt_out_len_actual, 
                "; expected", pt_out_len_expected, "\n")
            
            review = paste(review, 
                           paste0("spacer_3_outside_vdj_length=", pt_out_len_actual,
                                  "/", pt_out_len_expected),
                           sep=ifelse(review=="", "", ";"))
        } 
        
    }
    
    if (grepl(pattern="[^ATGCatgc]", x=pt_out)) {
        cat("non-ATGC in 3' spacer downstream of vdj_obsv_end:", pt_out, "\n")
        review = paste(review, "spacer_3_outside_vdj_non-ATGC", 
                       sep=ifelse(review=="", "", ";"))
    }
    
    return_vec = c(spacer_start, pt_out, review)
    names(return_vec) = c("spacer_start", "pt_out", "review")
    return(return_vec)
}



#' Prepare 3' spacer for mm lambda light chain for Ab expression
#' 
#' @param  vdj_obsv      "Cleaned" observed VDJ sequence. The expectation is that 
#'                       this was returned as `obsv_clean` by `clean_obsv`.
#' @param  full_seq      Input sequence. The `sequence` column in AIRR format.
#' @param  cdr3          Nt sequence of CDR3 (note: NOT IMGT-defined junction).
#' 
#' @return  A vector containing
#'          - cdr3_end_vdj_obsv: ending nt position of `cdr3` wrt `vdj_obsv`.
#'          - cdr3_end_full_seq: ending nt position of `cdr3` wrt `full_seq`.
#'          - junc_anchor_3_prime_nt: 3' junction anchor in nt.
#'          - junc_anchor_3_prime_aa: 3' junction anchor in aa.
#'          - bits_after_spacer_pt_1_actual_nt: bits of sequence in nt immediately
#'                                              after `spacer_pt_1` in `full_seq`.  
#'          - bits_after_spacer_pt_1_actual_aa: bits of sequence in aa immediately
#'                                              after `spacer_pt_1` in `full_seq`.
#'          - spacer_pt_1_vdj_obsv: the part of 3' spacer extracted from `vdj_obsv` 
#'                                  (immeidately after `cdr3` ).
#'          - spacer_pt_1_full_seq: 3' spacer extracted from `full_seq`.
#'          - vdj_obsv_up_to_cdr3_end: `vdj_obsv` up to and including the end of `cdr3`.
#'          - review: notes where manual review may be required.
#' 
#' @details  Only for mm lambda light chains. 
#'           Designed to prepare the following construct:
#'           
#'           CDR3 + 3' junction anchor + spacer of a fixed length + human IGLC2 CDS + XhoI
#' 
#'           in which,
#'           - CDR3 is the same as `cdr3`;
#'           - 3' junction anchor is captured by `junc_anchor_3_prime_nt`; and
#'           - spacer of a fixed length has an expected length specified by `spacer_pt_1_len_expected`
#'             and is captured by `spacer_pt_1_vdj_obsv` and `spacer_pt_1_full_seq`.
#'             
#'           (human IGLC2 CDS and XhoI are handled outside the scope of this function)
#'           
#'           Note: When deriving `bits_after_spacer_pt_1_actual_nt` and `bits_after_spacer_pt_1_actual_aa`,
#'           the full-length of `spacer_pt_1` as specified by `spacer_pt_1_len_expected` is used, 
#'           disregarding the actual lengths of `spacer_pt_1_vdj_obsv` and `spacer_pt_1_full_seq`.
#'           
#'           If the 3' junction anchor is not Phe (TTT/TTC) or Trp (TGG), 
#'           a review note is added.
#'           
#'           If there're not enough positions in `vdj_obsv` for `spacer_pt_1_vdj_obsv`, 
#'           a review note is added. 
#'           If `spacer_pt_1_vdj_obsv` has nt length that is not a multiple of 3, 
#'           a review note is added.
#'           
#'           If there're not enough positions in `full_seq` for `spacer_pt_1_full_seq`, 
#'           a review note is added. 
#'           If `spacer_pt_1_full_seq` has nt length that is not a multiple of 3, 
#'           a review note is added.
#'           
#'           If the bits immediately after `spacer_pt_1` in `full_seq` do not match 
#'           `bits_after_spacer_pt_1_expected_aa`, a review note is added. 
#'           If there're not enough positions in `full_seq` for such bits, 
#'           a review note is added.
#'   
prep_3_prime_spacer_mm_lambda = function(vdj_obsv, full_seq, cdr3) {
    
    require(stringi)
    require(alakazam)
    
    #* expected 3' anchor of junction
    #  TTT/TTC = Phe = F; TGG = Trp = W
    junc_anchor_3_prime_nt_expected = c("TTT", "TTC", "TGG")

    #* nt length of part of spacer to be derived from obsv
    #  9 aa
    spacer_pt_1_len_expected = 9*3
    
    #* expected bits immediately after spacer_pt_1
    #  3 aa
    bits_len_nt = 3*3
    bits_after_spacer_pt_1_expected_aa = "GQP"
    
    review = ""
    
    ### locate cdr3
    # returns a matrix
    # no match:
    #      start end
    # [1,]    NA  NA
    # match:
    #      start end
    # [1,]   289 330
    
    ## wrt vdj_obsv
    
    loc_cdr3_vdj_obsv = stri_locate(str=vdj_obsv, fixed=cdr3)
    # expect match
    stopifnot(!any(is.na(loc_cdr3_vdj_obsv)))
    # expect exactly 1 match
    stopifnot(nrow(loc_cdr3_vdj_obsv)==1)
    
    cdr3_end_vdj_obsv = loc_cdr3_vdj_obsv[1, 2]
    vdj_obsv_up_to_cdr3_end = substr(vdj_obsv, 1, cdr3_end_vdj_obsv)
    
    ## wrt full_seq
    
    loc_cdr3_full_seq = stri_locate(str=full_seq, fixed=cdr3)
    # expect match
    stopifnot(!any(is.na(loc_cdr3_full_seq)))
    # expect exactly 1 match
    stopifnot(nrow(loc_cdr3_full_seq)==1)
    
    cdr3_end_full_seq = loc_cdr3_full_seq[1, 2]
    
    
    ### 3' anchor
    
    junc_anchor_3_prime_nt = substr(vdj_obsv, 
                                    cdr3_end_vdj_obsv+1, 
                                    cdr3_end_vdj_obsv+3)
    stopifnot(nchar(junc_anchor_3_prime_nt)==3)
    junc_anchor_3_prime_aa = translateDNA(junc_anchor_3_prime_nt, trim=F)
    
    if (!junc_anchor_3_prime_nt %in% junc_anchor_3_prime_nt_expected) {
        
        cat("junc_anchor_3_prime_nt =", junc_anchor_3_prime_nt, 
            " (unexpected)", "\n")
        
        review = paste(review, 
                       paste0("junc_anchor_3_prime_nt=", junc_anchor_3_prime_nt,
                              "!=", 
                              paste(junc_anchor_3_prime_nt_expected, collapse="/")),
                       sep=ifelse(review=="", "", ";"))
        
    }
    
    ### part of spacer from within vdj
    # derive from vdj_obsv
    
    # between 3' anchor AA of junction and start of human IGLC2 CDS
    # length determined by spacer_pt_1_len_expected
    # +4 to skip anchor AA (downstream of junction)
    spacer_start_vdj_obsv = loc_cdr3_vdj_obsv[1, 2]+4
    
    spacer_pt_1_vdj_obsv = substr(vdj_obsv, spacer_start_vdj_obsv, 
                                  spacer_start_vdj_obsv+spacer_pt_1_len_expected-1)
    spacer_pt_1_vdj_obsv_len_actual = nchar(spacer_pt_1_vdj_obsv)
    stopifnot(spacer_pt_1_vdj_obsv_len_actual<=spacer_pt_1_len_expected)
    
    if (spacer_pt_1_vdj_obsv_len_actual < spacer_pt_1_len_expected) {
        
        cat("nchar(spacer_pt_1_vdj_obsv) =", spacer_pt_1_vdj_obsv_len_actual, 
            "; expected", spacer_pt_1_len_expected, "\n")
        
        review = paste(review, 
                       paste0("spacer_pt_1_vdj_obsv_length=", spacer_pt_1_vdj_obsv_len_actual,
                              "/", spacer_pt_1_len_expected),
                       sep=ifelse(review=="", "", ";"))
    }
    if (spacer_pt_1_vdj_obsv_len_actual%%3!=0) {
        
        cat("nchar(spacer_pt_1_vdj_obsv) =", spacer_pt_1_vdj_obsv_len_actual, 
            "%%3!=0", "\n")
        
        review = paste(review, 
                       paste0("spacer_pt_1_vdj_obsv_length=", spacer_pt_1_vdj_obsv_len_actual,
                              "%%3!=0"),
                       sep=ifelse(review=="", "", ";"))
    }
    
    ### part of spacer from within vdj
    # derive from full_seq
    # This could come in handy in case there's trimming in vdj_obsv (as part of 
    # clean_obsv() due to the germline seq), which could lead to there being not
    # enough positions in vdj_obsv for extracting the full-length spacer; whereas 
    # at the same time if one were to extract from full_seq for the full-length spacer,
    # the bits immediately after spacer are still of expected aa identifies
    # E.g. 
    # spacer derived from vdj_obsv has 21 nt; that from full_seq has 27 nt;
    # bits in full_seq immediately after the 27-nt spacer are GQP
    
    spacer_start_full_seq = loc_cdr3_full_seq[1, 2]+4
    
    spacer_pt_1_full_seq = substr(full_seq, spacer_start_full_seq, 
                                  spacer_start_full_seq+spacer_pt_1_len_expected-1)
    spacer_pt_1_full_seq_len_actual = nchar(spacer_pt_1_full_seq)
    stopifnot(spacer_pt_1_full_seq_len_actual<=spacer_pt_1_len_expected)
    
    if (spacer_pt_1_full_seq_len_actual < spacer_pt_1_len_expected) {
        
        cat("nchar(spacer_pt_1_full_seq) =", spacer_pt_1_full_seq_len_actual, 
            "; expected", spacer_pt_1_len_expected, "\n")
        
        review = paste(review, 
                       paste0("spacer_pt_1_full_seq_length=", spacer_pt_1_full_seq_len_actual,
                              "/", spacer_pt_1_len_expected),
                       sep=ifelse(review=="", "", ";"))
    }
    if (spacer_pt_1_full_seq_len_actual%%3!=0) {
        
        cat("nchar(spacer_pt_1_full_seq) =", spacer_pt_1_full_seq_len_actual, 
            "%%3!=0", "\n")
        
        review = paste(review, 
                       paste0("spacer_pt_1_full_seq_length=", spacer_pt_1_full_seq_len_actual,
                              "%%3!=0"),
                       sep=ifelse(review=="", "", ";"))
    }
    
    
    ### first 3 aa after spacer_pt_1
    # derive from full_seq in case vdj_obsv has run out
    
    # important to include +3 to account for 3' junction anchor
    bits_start = cdr3_end_full_seq+3+spacer_pt_1_len_expected+1
    bits_end = bits_start+bits_len_nt-1
    bits_after_spacer_pt_1_actual_nt = substr(full_seq, bits_start, bits_end)
    # translateDNA turns "" into NA
    bits_after_spacer_pt_1_actual_aa = translateDNA(bits_after_spacer_pt_1_actual_nt, trim=F)
    
    if (nchar(bits_after_spacer_pt_1_actual_nt)!=bits_len_nt) {
        
        cat("nchar(bits_after_spacer_pt_1_actual_nt) =", 
            nchar(bits_after_spacer_pt_1_actual_nt), 
            "; expected", bits_len_nt, "\n")
        
        review = paste(review, 
                       paste0("bits_after_spacer_pt_1_actual_nt_length=", 
                              nchar(bits_after_spacer_pt_1_actual_nt),
                              "/", bits_len_nt),
                       sep=ifelse(review=="", "", ";"))
        
    }
    if (is.na(bits_after_spacer_pt_1_actual_aa) || 
        bits_after_spacer_pt_1_actual_aa!=bits_after_spacer_pt_1_expected_aa) {
        
        cat("bits_after_spacer_pt_1_actual_aa =", 
            bits_after_spacer_pt_1_actual_aa,
            "; expected", bits_after_spacer_pt_1_expected_aa, "\n")
        
        review = paste(review, 
                       paste0("bits_after_spacer_pt_1_actual_aa=", 
                              bits_after_spacer_pt_1_actual_aa,
                              "!=", bits_after_spacer_pt_1_expected_aa),
                       sep=ifelse(review=="", "", ";"))
        
    }
    
    return_vec = c(cdr3_end_vdj_obsv, cdr3_end_full_seq, 
                   junc_anchor_3_prime_nt, junc_anchor_3_prime_aa,
                   bits_after_spacer_pt_1_actual_nt, bits_after_spacer_pt_1_actual_aa,
                   spacer_pt_1_vdj_obsv, spacer_pt_1_full_seq,
                   vdj_obsv_up_to_cdr3_end,
                   review)
    names(return_vec) = c("cdr3_end_vdj_obsv", "cdr3_end_full_seq", 
                          "junc_anchor_3_prime_nt", "junc_anchor_3_prime_aa",
                          "bits_after_spacer_pt_1_actual_nt", "bits_after_spacer_pt_1_actual_aa",
                          "spacer_pt_1_vdj_obsv", "spacer_pt_1_full_seq",
                          "vdj_obsv_up_to_cdr3_end",
                          "review")
    return(return_vec)
}


# `ambiguous=TRUE` translates degeneracy-based ambiguous bases (instead of "X")
#   e.g. GGN -> G (Gly, encoded by GG[ATGC])
#        GAN -> X (GA[TC] -> D; GA[AG] -> E)

# in-frame triplet "..." is turned into "." (instead of "X")
# in-frame triplet "---" is turned into "-" (instead of "X")
# any other triplet containing one or more of "." and/or "-" is turned into "X"
# triplet containing no "." or "-" is translated per seqinr::translate

# sequence input with length less than 3 is returned as NA
# NA is returned as NA

# helper function to translate a single triplet
translate_triplet = function(triplet, ambiguous) {
    require(seqinr)
    
    if (triplet=="...") {
        return(".")
    } else if (triplet=="---") {
        return("-")
    } else if (grepl(pattern="[-.]", x=triplet)) {
        return("X")
    } else {
        return(seqinr::translate(seq=s2c(triplet), 
                                 numcode=1, NAstring="X",
                                 ambiguous=ambiguous))
    }
}

translate_one_seq = function(sequence, ambiguous=TRUE) {
    require(stringi)
    
    if (is.na(sequence)) {
        return(NA)
    } else {
        stopifnot(is.character(sequence))
        
        nchar_seq = stri_length(sequence)
        
        if (nchar_seq >= 3) {
            # trim to length being a multiple of 3
            nchar_init_mod3 = nchar_seq %% 3
            
            if (nchar_init_mod3 != 0) {
                sequence = substr(sequence, 1, nchar_seq-nchar_init_mod3)
                # recalculate length
                nchar_seq = stri_length(sequence)
            }
            
            # split into triplets
            vec_triplet_start_index = seq(from=1, to=nchar_seq, by=3)
            vec_triplet = stri_sub(str=sequence, 
                                   from=vec_triplet_start_index, 
                                   to=(vec_triplet_start_index+2))
            
            # translate multiple triplets
            vec_triplet_trans = sapply(vec_triplet, translate_triplet, 
                                       ambiguous=ambiguous,
                                       USE.NAMES=F)
            # concat
            str_trans = stri_join(vec_triplet_trans, collapse="")
            
            return(str_trans)
            
        } else {
            # input length <3
            return(NA)
        }
    }
}

# translate_one_seq(NA)       # NA
# translate_one_seq("A")      # NA
# translate_one_seq("AA")     # NA
# translate_one_seq("AAA")    # K
# translate_one_seq("AAAA")   # K (trimmed)
# translate_one_seq("AAAAA")  # K (trimmed)
# translate_one_seq("AAAAAA") # KK
# translate_one_seq("AAA...") # K.
# translate_one_seq("AAA..A") # KX
# translate_one_seq("AAA..-") # KX
# translate_one_seq("AAA.A-") # KX
# translate_one_seq("AAA---") # K-
# translate_one_seq("AAAA--") # KX
# translate_one_seq("GGNGGN") # GG
# translate_one_seq("GGNGAN") # GX

# run translate_one_seq for a vector of nucleotide sequences
translate_dna = function(vec_sequence, ambiguous=TRUE) {
    require(seqinr)
    require(stringi)
    
    vec_trans = sapply(vec_sequence, translate_one_seq, ambiguous=ambiguous,
                       USE.NAMES=F)
    return(vec_trans)
}


# translate_dna(vec_sequence=c(NA,       # NA
#                              "A",      # NA
#                              "AA",     # NA
#                              "AAA",    # K
#                              "AAAA",   # K (trimmed)
#                              "AAAAA",  # K (trimmed)
#                              "AAAAAA", # KK
#                              "AAA...", # K.
#                              "AAA..A", # KX
#                              "AAA..-", # KX
#                              "AAA.A-", # KX
#                              "AAA---", # K-
#                              "AAAA--", # KX
#                              "GGNGGN", # GG
#                              "GGNGAN"),# GX
#               ambiguous=T)


# calculate all three possible reading frames
# Input
# - nt_seq: nucleotide sequence; a string
# - frames: an integer vector; values can be one or more of {1,2,3}
# Output
# A data.frame, with rows corresponding a frame

get_frames = function(nt_seq, frames=c(1:3)) {
    stopifnot(length(frames)>0 & 
                  length(frames)<=3 & 
                  all(frames %in% 1:3))
    
    nt_seq_len = nchar(nt_seq)
    
    df_cols = c("frame", "seq")
    df = data.frame(matrix(NA, nrow=length(frames), ncol=length(df_cols)))
    colnames(df) = df_cols
    
    df[["frame"]] = frames
    df[["seq"]] = sapply(frames, 
                         function(fr) {
                             if (fr<=nt_seq_len) {
                                 return(substring(nt_seq, first=fr))
                             } else {
                                 return(NA)
                             }
                         },
                         simplify=T, USE.NAMES=F)
    return(df)
}

# get_frames("A", frames=c(1:3)) # expect NA for frames 2:3
# get_frames("AT", frames=c(1:3)) # expect NA for frame 3
# get_frames("ATG", frames=c(1:3)) 
# get_frames("ATGC", frames=c(1:3)) 
# get_frames("ATGC", frames=c(1)) 
# get_frames("ATGC", frames=c(1:4)) # expects error
