# script that helpts to prepare auxiliary files used for launching 
# cellranger runs on RIS

# see main/z_cr_aux_prep.R for SOP using this function

# Given a specific type of csv file from MGI,
# generate a data.frame containing columns that can be used to 
# easily populate aux files necessary for launching cellranger runs on RIS
# Input
#  - fn_csv: path to csv file from MGI
#  - library_type_csv: library type as specific in the corresponding column in the csv
#  - library_type_remark: a simplified label for the same library type to be used
#                         in the remark note that is to be generated
# Output
#  - A data.frame in which
#      Each row is a unique combination of {Library.Name, ESP.ID, Index.Sequence} (derived from CSV)
#      Each row also has theses added columns: aux_column; remark
# 
# Expected input csv file from MGI
# "Samplemap2.csv"
# Each row represents a library
# Useful columns:
#  - FASTQ Path - Read 1
#  - FASTQ Path - Read 2
#  - Flowcell ID
#  - Index Sequence
#  - Flowcell Lane
#  - ESP ID
#  - Library Type
#  - Library Name

generate_aux = function(fn_csv, library_type_csv, library_type_remark) {

    # columns expected to be present
    # " ", "-" gets converted to "."
    vec_expected = c("FASTQ.Path...Read.1", "FASTQ.Path...Read.2",
                     "Flowcell.ID", "Index.Sequence", "Flowcell.Lane",
                     "ESP.ID", "Library.Type", "Library.Name")
    
    if (file.exists(fn_csv)) {
        df_csv = read.table(fn_csv, sep=",", header=T, stringsAsFactors=F)
    } else {
        stop(fn_csv, " does not exist")
    }
    
    bool_vec_expected = vec_expected %in% colnames(df_csv)
    if (!all(bool_vec_expected)) {
        msg = paste0("Column missing: ", colnames(df_csv)[!bool_vec_expected])
        for (m in msg) {
            message(m, "\n")
            stop("Not all expected columns are present.")
        }
    }
    
    # subset columns
    df_csv = df_csv[, vec_expected]
    
    # subset to specified library_type
    if (!library_type_csv %in% df_csv[["Library.Type"]]) {
        stop(library_type_csv, " not found in csv")
    }
    
    bool_library_type = df_csv[["Library.Type"]]==library_type_csv
    
    df_csv = df_csv[bool_library_type, ]
    
    # These 3 columns should be redundant
    df_csv_uniq = unique(df_csv[, c("Library.Name", "ESP.ID", "Index.Sequence")])
    stopifnot(!any(duplicated(df_csv_uniq[["Library.Name"]])))
    stopifnot(!any(duplicated(df_csv_uniq[["ESP.ID"]])))
    stopifnot(!any(duplicated(df_csv_uniq[["Index.Sequence"]])))
    
    # check FASTQ file names
    # expected file name format: [ESP.ID]_[Flowcell.ID]_S*_L*_R[12]_001.fastq.gz
    .parse_fastq_fn = function(s) {
        s_split = strsplit(s, "_")[[1]]
        stopifnot(length(s_split)==6)
        return(s_split)
    }
    
    # each list entry corresponds to a row in df_csv
    lst_fastq_fn_parsed = sapply(df_csv[["FASTQ.Path...Read.1"]],
                                 .parse_fastq_fn, 
                                 simplify=F, USE.NAMES=F)
    df_fastq_fn_parsed = do.call(rbind, lst_fastq_fn_parsed)
    
    stopifnot(all.equal( df_fastq_fn_parsed[, 1], df_csv[["ESP.ID"]] ))
    stopifnot(all.equal( df_fastq_fn_parsed[, 2], df_csv[["Flowcell.ID"]] ))
    
    # unique library names
    vec_library_name_uniq = sort(df_csv_uniq[["Library.Name"]])
    # placeholder
    df_csv_uniq[["aux_column"]] = NA
    df_csv_uniq[["remark"]] = NA
    
    for (i_lib in 1:length(vec_library_name_uniq)) {
        
        cur_library_name = vec_library_name_uniq[i_lib]
        cat(cur_library_name, "\n")
        cur_idx = which(df_csv[["Library.Name"]]==cur_library_name)
        
        # sanity check 
        # actually redundant to the checks based on the 3 columns earlier 
        stopifnot(length(unique(df_csv[["ESP.ID"]][cur_idx]))==1)
        stopifnot(length(unique(df_csv[["Index.Sequence"]][cur_idx]))==1)
        
        # there could be multiple flowcells
        cur_vec_esp_id_flowcell_id = paste0(df_fastq_fn_parsed[cur_idx, 1], "_",
                                            df_fastq_fn_parsed[cur_idx, 2])
        cur_vec_esp_id_flowcell_id_uniq = sort(unique(cur_vec_esp_id_flowcell_id))
        
        cur_n_flowcell = length(unique(df_csv[["Flowcell.ID"]][cur_idx]))
        stopifnot(length(cur_vec_esp_id_flowcell_id_uniq) == cur_n_flowcell)
        
        # string together
        cur_str_esp_id_flowcell_id_uniq = paste0(cur_vec_esp_id_flowcell_id_uniq, 
                                                 collapse=",")
        
        df_csv_uniq[["aux_column"]][i_lib] = cur_str_esp_id_flowcell_id_uniq
        
        cur_n_lane = nrow(unique(df_csv[cur_idx, c("Flowcell.ID", "Flowcell.Lane")]))
        
        cur_remark = paste0(library_type_remark, "_split_over_",
                            cur_n_flowcell, "_flow_",
                            ifelse(cur_n_flowcell>1, "cells", "cell"),
                            "_", cur_n_lane, "_",
                            ifelse(cur_n_lane>1, "lanes", "lane"))
        
        df_csv_uniq[["remark"]][i_lib] = cur_remark
    }
    
    stopifnot(!any(is.na(df_csv_uniq[["aux_column"]])))
    
    return(df_csv_uniq)
}


