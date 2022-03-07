#' Collect metrics_summary.csv's from 10x `outs/` directories
#' 
#' @param  path_csv     A header-less, two-column CSV containing paths to
#'                      metrics_summary.csv's. The first and second columns 
#'                      indicate library names and paths respectively.
#' @param  output_tsv   Filename for output TSV (note: not CSV). Can be an 
#'                      absolute path, a relative path, or just a filename.
#' @param  overwrite    Boolean. Whether to overwrite existing `output_tsv`.
#'                      Defaults to `FALSE`.
#'                                            
#' @return A TSV file containing information from all the input metrics_summary.csv's.
#'         Each row corresponds to a 10x metrics_summary.csv.
#'         
#'         If `output_tsv` is just a filename (as opposed to an absolute or a relative
#'         path), the TSV file is exported to the directory from which the function
#'         is called.
#'         
#'         If `output_tsv` already exists and `overwrite` is set to `FALSE`, a
#'         warning will be printed and no TSV exported.
#' 
#' @details 
#' An example of what the content in `path_csv` might look like:
#' s104,/storage1/cr_bcr/s104/outs/metrics_summary.csv
#' s105,/storage1/cr_bcr/s105/outs/metrics_summary.csv
#' 
#' @example 
#' collect_metrics_summary(path_csv="cr_path_metrics_bcr.csv", 
#'                         output_tsv="all_metrics_bcr.tsv")
#' 
collect_metrics_summary = function(path_csv, output_tsv, overwrite=FALSE) {
    
    # check that path_csv exists
    stopifnot(file.exists(path_csv))
    
    tab_input = read.table(path_csv, header=F, sep=",", stringsAsFactors=F)
    
    # two-column CSV
    stopifnot(ncol(tab_input)==2)
    
    # library names
    libs = tab_input[, 1]
    
    # check that all netrucs_summary.csv's exist
    bool_exist = sapply(tab_input[, 2], file.exists)
    if (!all(bool_exist)) {
        cat("CSV(s) not found for:\n")
        for (s in libs[!bool_exist]) {
            cat(s, "\n")
        }
        stop("Missing CSV(s)")
    }
    
    # holder 
    lst = vector(mode="list", length=length(libs))
    
    # read in individual metrics_summary.csv's
    for (i_row in 1:nrow(tab_input)) {
        cur_csv = read.table(tab_input[i_row, 2], sep=",", header=T, stringsAsFactors=F)
        # store in holder
        lst[[i_row]] = cur_csv
    }
    
    # rbind
    all_csv = do.call(rbind, lst)
    
    # append library names as first column 
    df = cbind(Library=libs, all_csv, stringsAsFactors=F)
    
    # sanity check: expect no NA anywhere
    stopifnot(!any(is.na(df)))
    
    # export
    if (file.exists(output_tsv)) {
        if (overwrite) {
            # overwrite if `output_tsv` already exists
            cat("`overwrite` set to TRUE; overwriting", output_csv, "\n")
            
            write.table(x=df, file=output_tsv, quote=F, sep="\t",
                        row.names=F, col.names=T)
        } else {
            warning("`overwrite` set to FALSE; already exists: ", output_tsv, " ; halted\n")
        }
    } else {
        write.table(x=df, file=output_tsv, quote=F, sep="\t",
                    row.names=F, col.names=T)
    }
    
}

