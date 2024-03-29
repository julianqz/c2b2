# Julian Q. Zhou
# https://github.com/julianqz

#' Circos plot for visualizating B cell clonal overlap
#'
#' @params vec_sectors        A vector specifying what each arc of the circos plot represents.
#' @params vec_sectors_type   A vector specifying the "type" of the arc. This "type" is
#'                           used to determine whether there's clonal overlap.
#' @params vec_sectors_color  The color of each arc.
#' @params vec_gaps           The amount of space between arcs.
#' @params compute_overlap    Boolean specifying whether to compute clonal overlap.
#' @params vec_overlap_types  At least 2 unique types from `vec_sectors_type` that define 
#'                           clonal overlap.
#' @params col_sector         Column name in `df_seq_data` that specifies `vec_sectors`.
#' @params col_clone_id       Column name in `df_seq_data` and `df_clone_info` that specifies
#'                           clone ID.
#' @params df_seq_data        Data.frame containing sequence data.
#' @params df_clone_info      Data.frame containing summary info on clones.
#' @params color_overlap      Color for chords where there's clonal overlap between "types".
#'                            Ignored if `compute_overlap` is `FALSE`.
#' @params color_no_overlap   Color for chords where there's no clonal overlap between "types".
#'
#' @details  The lengths of `vec_sectors`, `vec_sectors_type`, `vec_sectors_color`, and
#'           `vec_gaps` must match.
#'
#'           `df_seq_data` and `df_clone_info` must contain the same set of clones. Only
#'            clones that should be included in the overlap analysis and for visualization 
#'            should be included.
#'
#'            For each entry in `vec_sectors` and `vec_sectors_type`, there must be a column
#'            of the same name in `df_clone_info`.
#'
#'            It is assumed that there are at least 2 "types" of arcs/sectors. 
#'            Clonal overlap is deemed to exist if there is connection between arcs of 
#'            exactly 2 different "types". 
#'            In practice, "type" often corresponds to compartment. When sectors/arcs 
#'            correspond to compartment-timepoint combinations, knowing the "types" of each
#'            sector/arc makes it easy (internally for the function) to determine if there's 
#'            clonal overlap (btw compartments) in the presence of an additional variable, timepoint.
#'            
#'            If `compute_overlap=FALSE`, no clonal overlap will be quantified or visualized.
#'            Accordingly, all chords will be drawn in `color_no_overlap`. The number and
#'            percentage of clonal overlap will be reported as 0 in the print out.
#'            `color_overlap` and `vec_overlap_types` must still be specified as an argument 
#'            but both will be ignored.
#'
circos_clonal_overlap = function(vec_sectors, vec_sectors_type, 
                                 vec_sectors_color, vec_gaps, 
                                 compute_overlap=TRUE,
                                 vec_overlap_types,
                                 col_sector, col_clone_id, 
                                 df_seq_data, df_clone_info, 
                                 color_overlap, color_no_overlap) {
    require(circlize) # v0.4.13
    
    # pre-checks
    
    stopifnot(length(vec_sectors) == length(vec_sectors_type))
    stopifnot(length(vec_sectors) == length(vec_sectors_color))
    stopifnot(length(vec_gaps) == length(vec_sectors))
    
    stopifnot(length(vec_overlap_types) <= length(vec_sectors))
    stopifnot(!any(duplicated(vec_overlap_types)))
    
    stopifnot( col_sector %in% colnames(df_seq_data) )
    stopifnot( col_clone_id %in% colnames(df_seq_data) )
    
    stopifnot( col_clone_id %in% colnames(df_clone_info) )
    stopifnot(all( vec_sectors_type %in% colnames(df_clone_info) ))
    stopifnot(all( vec_sectors %in% colnames(df_clone_info) ))
    stopifnot(all( vec_overlap_types %in% vec_sectors_type ))
    
    # every clone in df_seq_data should have an entry in df_clone_info 
    stopifnot( nrow(df_clone_info) == length(unique(df_seq_data[[col_clone_id]])) )
    stopifnot( all( unique(df_seq_data[[col_clone_id]]) %in% df_clone_info[[col_clone_id]] ) )
    
    # expect at least 2 types of sectors
    stopifnot(length(unique(vec_sectors_type))>=2)
    # expect exactly 2 types of sectors for determining whether there's overlap
    stopifnot(length(unique(vec_overlap_types))==2)
    
    # overlap
    
    if (compute_overlap) {
        # each list item corresponds to a sector type in vec_overlap_types
        # each list item contains a vector; each vector item corresponds to a clone
        bool_types_lst = sapply(vec_overlap_types, function(cur_type){
            cur_sectors = vec_sectors[vec_sectors_type==cur_type]
            if (length(cur_sectors)==1) {
                return(df_clone_info[[cur_sectors]]>0)
            } else {
                return(rowSums(df_clone_info[, cur_sectors])>0)
            }
        }, simplify=F)
        # row: clone; col: sector type
        bool_types_mtx = do.call(cbind, bool_types_lst)
        # whether there's overlap (whether a clone contains both sector types)
        bool_overlap = rowSums(bool_types_mtx)==length(vec_overlap_types) 
        vec_overlap_clones = df_clone_info[[col_clone_id]][bool_overlap]
        cat("\ncompute_overlap set to TRUE\n\n")
    } else {
        vec_overlap_clones = c()
        cat("\ncompute_overlap set to FALSE\n\n")
    }
    
    
    # sequences for each sector
    lst = vector(mode="list", length=length(vec_sectors))
    names(lst) = vec_sectors
    for (s in vec_sectors) {
        
        # use which() in case there's NA in $col_sector
        cur_db = df_seq_data[which(df_seq_data[[col_sector]]==s), ]
        
        if (nrow(cur_db)>0) {
            
            col_cur_sector = s
            col_other_sectors_of_interest = vec_sectors[-which(vec_sectors==s)]
            
            cols_sectors = c(list(col_cur_sector), col_other_sectors_of_interest)
            names(cols_sectors)[1] = s
            
            # for each clone, whether it has seq in each sector of interest
            bools_sectors = lapply(cols_sectors, function(cols){
                if (length(cols)==1) {
                    return(df_clone_info[[cols]]>0)
                } else {
                    return(rowSums(df_clone_info[, cols])>0)
                }
            })
            # row: clone
            # col: sector
            bools_sectors = do.call(cbind, bools_sectors)
            # for each clone, whether it has seq in current sector, 
            # AND seq in any of the other sectors of interest
            # first col always corresponds to current sector
            if (nrow(bools_sectors)==1) {
                if (ncol(bools_sectors)==2) {
                    bools_sectors = bools_sectors[1] & bools_sectors[2]
                } else {
                    bools_sectors = bools_sectors[1] & any(bools_sectors[-1])
                }
            } else {
                if (ncol(bools_sectors)==2) {
                    bools_sectors = bools_sectors[, 1] & bools_sectors[, 2]
                } else {
                    bools_sectors = bools_sectors[, 1] & apply(bools_sectors[, -1], 1, any)
                }
            }
            
            
            if (any(bools_sectors)) {
                # initialize temporary $overlap column
                # do this before splitting data.frame in case cur_db_one_sector_only has 0 row
                # use "A" for TRUE, "B" for FALSE
                # order(c(TRUE, FALSE)) would put FALSE ahead of TRUE; want TRUE ahead of FALSE
                cur_db[["overlap"]] = "B"
                
                clones_multi_sectors = df_clone_info[[col_clone_id]][bools_sectors]
                
                cur_db_one_sector_only = cur_db[!cur_db[[col_clone_id]] %in% clones_multi_sectors, ]
                
                cur_db_multi_sectors = cur_db[cur_db[[col_clone_id]] %in% clones_multi_sectors, ]
                
                cur_db_multi_sectors[["overlap"]] = ifelse(cur_db_multi_sectors[[col_clone_id]] %in% vec_overlap_clones,
                                                           "A", "B")
                
                # if multi-sector, order first by whether there's clonal overlap, then by clone id
                cur_db_multi_sectors_order = cur_db_multi_sectors[order(cur_db_multi_sectors[["overlap"]], 
                                                                        cur_db_multi_sectors[[col_clone_id]]), ]
                
                cur_db_order = rbind(cur_db_one_sector_only, cur_db_multi_sectors_order)
            } else {
                cur_db_order = cur_db[order(cur_db[[col_clone_id]]), ]
            }
            
            lst[[s]] = cur_db_order
            
            # print out stat about overlap
            cur_db_uniq_clones = unique(cur_db_order[[col_clone_id]])
            cur_db_uniq_clones_bool_overlap = cur_db_uniq_clones %in% vec_overlap_clones
            cat(s, ": total # clones =", 
                length(cur_db_uniq_clones), 
                "; # clones w/ overlap =", 
                sum(cur_db_uniq_clones_bool_overlap),
                "(",  
                round(mean(cur_db_uniq_clones_bool_overlap)*100, 3),
                "%)", "\n")
        }
    }
    
    lst_bool = !unlist(lapply(lst, is.null)) 
    lst = lst[lst_bool]
    
    mtx_sectors_xlim = do.call(rbind, lapply(lst, function(l){ return(c(1, nrow(l))) } ))
    
    #* hack
    # c(1,1) will cause failure
    if (any(mtx_sectors_xlim[, 2]==1)) {
        mtx_sectors_xlim[, 2][ mtx_sectors_xlim[, 2]==1 ] = 1.5
    }
    
    print(mtx_sectors_xlim)
    
    vec_sectors = vec_sectors[lst_bool]
    vec_sectors_type = vec_sectors_type[lst_bool]
    vec_gaps = vec_gaps[lst_bool]
    vec_sectors_color = vec_sectors_color[lst_bool]
    
    # factorize
    vec_sectors_factor = factor(x=vec_sectors, levels=vec_sectors)
    
    # initialize plot
    circos.par(start.degree=90, clock.wise=T, gap.after=vec_gaps, 
               cell.padding = c(0.02, 0, 0.02, 0))
    
    circos.initialize(factors=vec_sectors_factor, xlim=mtx_sectors_xlim)
    
    circos.trackPlotRegion(ylim=c(0, 1), bg.col=vec_sectors_color, 
                           track.height=uh(3.5, "mm"), bg.border=NA)
    
    # draw links
    for (i in 1:length(vec_sectors)) {
        for (j in 1:length(vec_sectors)) {
            if (j>i) {
                #cat(i, j, "\n")
                i_sect = vec_sectors[i]
                j_sect = vec_sectors[j]
                
                # wrt df_clone_info columns
                i_sect_col_idx = which(colnames(df_clone_info)==i_sect)
                j_sect_col_idx = which(colnames(df_clone_info)==j_sect)
                
                # expect exactly 1 match
                stopifnot(length(i_sect_col_idx)==1)
                stopifnot(length(j_sect_col_idx)==1)
                
                i_bool = df_clone_info[, i_sect_col_idx]>0
                j_bool = df_clone_info[, j_sect_col_idx]>0
                
                clone_bool = i_bool & j_bool
                
                if (any(clone_bool)) {
                    cur_clones = df_clone_info[[col_clone_id]][clone_bool]
                    
                    for (cl in cur_clones) {
                        
                        if (cl %in% vec_overlap_clones) {
                            color_link = color_overlap
                        } else {
                            color_link = color_no_overlap
                        }
                        
                        idx_i = which( lst[[i_sect]][[col_clone_id]]==cl )
                        idx_j = which( lst[[j_sect]][[col_clone_id]]==cl )
                        
                        if (length(idx_i)==1) {
                            rg_i = c(idx_i, idx_i)
                        } else {
                            rg_i = c(idx_i[1], idx_i[length(idx_i)])
                        }
                        if (length(idx_j)==1) {
                            rg_j = c(idx_j, idx_j)
                        } else {
                            rg_j = c(idx_j[1], idx_j[length(idx_j)])
                        }
                        
                        circos.link(i_sect, rg_i, j_sect, rg_j, 
                                    col=scales::alpha(color_link, 0.5),
                                    #border=NA
                                    border=scales::alpha(color_link, 0.7), 
                                    lwd=0.3)
                        
                    }
                }
                
            }
        }
    }
    
    # clean up
    circos.clear()
}

