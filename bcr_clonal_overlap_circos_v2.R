# Julian Q. Zhou
# https://github.com/julianqz

#' Circos plot for visualizating B cell clonal overlap
#'
#' @params vec_sectors        A vector specifying what each arc of the circos plot represents.
#' @params vec_sectors_type   A vector specifying the "type" of the arc. This "type" is
#'                            used to determine whether there's clonal overlap.
#' @params vec_sectors_color  The color of each arc.
#' @params vec_gaps           The amount of space between arcs.
#' @params col_sector         Column name in `df_seq_data` that specifies `vec_sectors`.
#' @params col_clone_id       Column name in `df_seq_data` and `df_clone_info` that specifies
#'                           clone ID.
#' @params df_seq_data        Data.frame containing sequence data.
#' @params df_clone_info      Data.frame containing summary info on clones.
#' @params color_overlap      Color for chords where there's clonal overlap between "types".
#'                            Named vector. The names specificy overlaps. 
#' @params color_rest         Color for chords where there's no clonal overlap between "types".
#' @params color_mono         If `TRUE`, use `color_rest` for coloring chords throughout 
#'                            regardless of overlap. Default is `FALSE`.
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
#'            different "types". 
#'            More than 2 "types" may be involved in a overlap, and such overlaps will be
#'            reported accordingly in the stats print out. 
#'            However, currently, the visualization only supports coloring chords based on
#'            overlap between exactly 2 "types". In other words, it doesn't support, for 
#'            instance, coloring a 3-way overlap with a specific color that is different 
#'            from the colors specified for pairwise 2-way overlaps between the three "types" 
#'            involved.
#'            
#'            When reporting overlaps, a 2-way overlap that is part of a 3-way overlap is 
#'            counted under the 2-way overlap as well. That is, if both "GC@MBC" and "GC@MBC@PB"
#'            are specified, a clone spanning GC, MBC, and PB will be counted under "GC@MBC"
#'            as well as under "GC@MBC@PB". In other words, "GC@MBC" does not mean GC and MBC
#'            only and nothing else.
#'            
#'            In practice, "type" often corresponds to compartment. When sectors/arcs 
#'            correspond to compartment-timepoint combinations, knowing the "types" of each
#'            sector/arc makes it easy (internally for the function) to determine if there's 
#'            clonal overlap (btw compartments) in the presence of an additional variable, timepoint.
#'            
#'            If `color_mono=TRUE`, the color specified by `color_rest` will be used to
#'            color all chords, regardless of overlap. The number of percentage of clonal
#'            overlap, as specified by `color_overlap`, will still be computed and reported. 
#'
circos_clonal_overlap_v2 = function(vec_sectors, vec_sectors_type, 
                                    vec_sectors_color, vec_gaps, 
                                    col_sector, col_clone_id, 
                                    df_seq_data, df_clone_info, 
                                    color_overlap,
                                    color_rest,
                                    color_mono=FALSE) {
    require(circlize) # v0.4.13
    
    # pre-checks
    
    stopifnot(length(vec_sectors) == length(vec_sectors_type))
    stopifnot(length(vec_sectors) == length(vec_sectors_color))
    stopifnot(length(vec_gaps) == length(vec_sectors))
    
    # color_overlap should be a named vector
    stopifnot(!is.null(names(color_overlap)))
    # all sector types in names(color_overlap) should be in vec_sectors_type
    stopifnot(all( unlist(sapply(names(color_overlap), 
                                 function(s){strsplit(s, "@")[[1]]}, 
                                 USE.NAMES=F, simplify=F)) %in% vec_sectors_type ))
    
    stopifnot( col_sector %in% colnames(df_seq_data) )
    stopifnot( col_clone_id %in% colnames(df_seq_data) )
    
    stopifnot( col_clone_id %in% colnames(df_clone_info) )
    stopifnot(all( vec_sectors_type %in% colnames(df_clone_info) ))
    stopifnot(all( vec_sectors %in% colnames(df_clone_info) ))
    
    # every clone in df_seq_data should have an entry in df_clone_info 
    stopifnot( nrow(df_clone_info) == length(unique(df_seq_data[[col_clone_id]])) )
    stopifnot( all( unique(df_seq_data[[col_clone_id]]) %in% df_clone_info[[col_clone_id]] ) )
    
    # expect at least 2 types of sectors
    stopifnot(length(unique(vec_sectors_type))>=2)
    
    
    ### sector types involved in coloring
    vec_sectors_type_to_color = sort(unique(unlist(sapply(names(color_overlap), 
                                                          function(x){strsplit(x, "@")[[1]]},
                                                          simplify=F))))
    
    ### binary matrix
    # for each clone, whether the clone spans a sector type
    # row: clone
    # col: vec_sectors_type_to_color
    
    # each list item corresponds to a sector type in vec_sectors_type_to_color
    # each list item contains a boolean vector; each vector item corresponds to a clone
    # bool vectors indicate whether each clone spans a sector type
    bool_types_lst = sapply(vec_sectors_type_to_color, function(cur_type){
        cur_sectors = vec_sectors[vec_sectors_type==cur_type]
        if (length(cur_sectors)==1) {
            return(df_clone_info[[cur_sectors]]>0)
        } else {
            return(rowSums(df_clone_info[, cur_sectors])>0)
        }
    }, simplify=F)
    # convert list to mtx
    # row: clone; col: sector type
    bool_types_mtx = do.call(cbind, bool_types_lst)
    
    
    RUN=F
    if (RUN) {
        
        # not used: this could result in "GC@MBC_PBMC@PB", which ends up being
        # treated separately from "GC@MBC_PBMC"
        # As a result, GC MBC_PBMC overlap is undercounted as it's not counted 
        # in cases of GC@MBC_PBMC@PB 
        
        # each item corresponds to a clone
        # observed overlap
        # if clone contains none of vec_sectors_type_to_color, returned value is ""
        vec_cl_overlap_def = rep("", nrow(bool_types_mtx))
        for (i_cl in 1:nrow(bool_types_mtx)) {
            cur_cl_bool = bool_types_mtx[i_cl, ]
            cur_cl_types_to_color = vec_sectors_type_to_color[cur_cl_bool]
            vec_cl_overlap_def[i_cl] = paste(sort(cur_cl_types_to_color), collapse="@")
        }
        table(vec_cl_overlap_def, useNA="ifany")
    
    
    
        # map from clone_info to db
        # wrt df_clone_info
        idx_cl_id = match(df_seq_data[[col_clone_id]],
                          df_clone_info[[col_clone_id]])
        stopifnot(!any(is.na(idx_cl_id)))
        stopifnot(all.equal( df_seq_data[[col_clone_id]], df_clone_info[[col_clone_id]][idx_cl_id] ))
        
        col_overlap_def = "OVERLAP_DEF"
        stopifnot(!col_overlap_def %in% colnames(df_seq_data))
        df_seq_data[[col_overlap_def]] = vec_cl_overlap_def[idx_cl_id]
        
    }
    
    
    ### binary matrix
    # for each clone, whether the clone contains a given overlap 
    # row: clone
    # col: clonal overlap definitions
    
    # clonal overlap definitions derived from color_overlap
     
    bool_types_mtx_pt_2 = matrix(NA, nrow=nrow(bool_types_mtx), 
                                 ncol=length(color_overlap))
    colnames(bool_types_mtx_pt_2) = names(color_overlap)
    
    # for each overlap definition
    for (cur_ov in names(color_overlap)) {
        # check just in case weird character exists in color_overlap that
        # prevents unchanged colnames assignment
        stopifnot(cur_ov %in% colnames(bool_types_mtx_pt_2))
        
        # derive definition
        cur_ov_sector_types = strsplit(cur_ov, "@")[[1]]
        stopifnot(all(cur_ov_sector_types %in% colnames(bool_types_mtx)))
        
        # whether ALL sector types in cur_ov_sector_types present
        # works for nrow=1 too
        bool_types_mtx_pt_2[, cur_ov] = rowSums(bool_types_mtx[, cur_ov_sector_types])==length(cur_ov_sector_types)

    }
    
    # check that there's no clash in colnames
    stopifnot(length(intersect(colnames(bool_types_mtx), 
                               colnames(bool_types_mtx_pt_2)))==0)
    
    ### observed overlap(s) for each clone
    
    # each item corresponds to a clone
    # If clone contains no overlap, returned value is ""
    # Overlap definitions come from names(color_overlap)
    # Possible for a clone to fulfill multiple overlap definitions
    # Multiple overlaps are separated by "%"
    
    vec_cl_overlap_def = rep("", nrow(bool_types_mtx_pt_2))
    
    for (i_cl in 1:nrow(bool_types_mtx_pt_2)) {
        cur_cl_bool = bool_types_mtx_pt_2[i_cl, ]
        cur_cl_overlaps = colnames(bool_types_mtx_pt_2)[cur_cl_bool]
        vec_cl_overlap_def[i_cl] = paste(sort(cur_cl_overlaps), collapse="%")
    }
    
    #print( table(vec_cl_overlap_def, useNA="ifany") )
    
    ### map from clone_info to db
    # wrt df_clone_info
    # recall that bool_types_mtx and vec_cl_overlap_def have 
    # the same clone order as df_clone_info
    idx_cl_id = match(df_seq_data[[col_clone_id]],
                      df_clone_info[[col_clone_id]])
    stopifnot(!any(is.na(idx_cl_id)))
    stopifnot(all.equal( df_seq_data[[col_clone_id]], df_clone_info[[col_clone_id]][idx_cl_id] ))
    
    col_overlap_def = "OVERLAP_DEF"
    stopifnot(!col_overlap_def %in% colnames(df_seq_data))
    df_seq_data[[col_overlap_def]] = vec_cl_overlap_def[idx_cl_id]
    
    
    ### compute overlap and report stats

    # each list entry corresponds to a sector
    # each list item is a db containing (reordered) sequences for each sector
    
    lst = vector(mode="list", length=length(vec_sectors))
    names(lst) = vec_sectors
    
    for (s in vec_sectors) {
        
        # use which() in case there's NA in $col_sector
        cur_db = df_seq_data[which(df_seq_data[[col_sector]]==s), ]
        
        # if there's any seq in current sector
        if (nrow(cur_db)>0) {
            # whether any seq in current sector has any overlap
            bool_overlap_to_color = grepl(pattern="@", x=cur_db[[col_overlap_def]], fixed=T)
            
            if (any(bool_overlap_to_color)) {
                
                # reorder seqs according to overlap definitions
                
                # one-sector
                # multi-sector
                #   color, then clone_id
                
                cur_db_multi_sectors_to_color = cur_db[bool_overlap_to_color, ]
                cur_db_no_sector_to_color = cur_db[!bool_overlap_to_color, ]
                
                # order first by overlap color, then by clone id
                cur_db_multi_sectors_to_color_order = cur_db_multi_sectors_to_color[order(cur_db_multi_sectors_to_color[[col_overlap_def]],
                                                                                          cur_db_multi_sectors_to_color[[col_clone_id]]), ]
            
                cur_db_order = rbind(cur_db_no_sector_to_color, cur_db_multi_sectors_to_color_order)
                
                #print( unique(cur_db_order[, c(col_clone_id, col_overlap_def)]) )
                
            } else {
                cur_db_order = cur_db[order(cur_db[[col_clone_id]]), ]
            }
            
            lst[[s]] = cur_db_order
            
            # print out overall stats about overlap
            cur_db_uniq_clones = unique(cur_db_order[[col_clone_id]])
            cur_db_ncl = length(cur_db_uniq_clones)
            cat(s, ": total # seq = ", nrow(cur_db_order),
                "; # clones = ", cur_db_ncl, 
                sep="", "\n")
            
            # de-collapse col_overlap_def
            # list order matches cur_db_order
            lst_overlap_def_decollapsed_cur_db_order = sapply(cur_db_order[[col_overlap_def]],
                                                              function(s){
                                                                  strsplit(s, "%")[[1]]
                                                              }, USE.NAMES=F, simplify=F)
            
            # print out stats for each overlap definition
            for (cur_def in names(color_overlap)) {
                
                # for each seq in sector
                # whether its clone fulfills the current overlap definition
                
                cur_def_bool = unlist(lapply(lst_overlap_def_decollapsed_cur_db_order,
                                      function(s){ return(cur_def %in% s) }))
                
                # number of seqs whose clones fulfill the current overlap definition
                cur_def_nseq = sum(cur_def_bool)
                # number of unique clones fulfilling the current overlap definition
                cur_def_ncl = length(unique(cur_db_order[[col_clone_id]][cur_def_bool]))
                
                cat(" - ", cur_def, 
                    ": # seq = ", cur_def_nseq, " (",
                    # % is of total number of seqs in current sector
                    round(cur_def_nseq/nrow(cur_db_order)*100, 3), "%)", 
                    "; # clone = ", cur_def_ncl, " (",
                    # % is of total number of unique clones spanning current sector
                    round(cur_def_ncl/cur_db_ncl*100, 3), "%)", 
                    sep="", "\n")
                
            }
            
            cat("\n")
            
        }
    }
    
    # remove sector with no seq
    lst_bool = !unlist(lapply(lst, is.null)) 
    lst = lst[lst_bool]
    
    
    ### sequence count by sector
    
    mtx_sectors_xlim = do.call(rbind, lapply(lst, function(l){ return(c(1, nrow(l))) } ))
    
    #* hack
    # c(1,1) [only 1 seq in sector] will cause failure in circilize plotting function
    if (any(mtx_sectors_xlim[, 2]==1)) {
        mtx_sectors_xlim[, 2][ mtx_sectors_xlim[, 2]==1 ] = 1.05
    }
    
    print(mtx_sectors_xlim)
    
    # subset associated variables also to only sectors with seq 
    vec_sectors = vec_sectors[lst_bool]
    vec_sectors_type = vec_sectors_type[lst_bool]
    vec_gaps = vec_gaps[lst_bool]
    vec_sectors_color = vec_sectors_color[lst_bool]
    
    # factorize
    vec_sectors_factor = factor(x=vec_sectors, levels=vec_sectors)
    
    ### visualization
    
    ## initialize plot
    circos.par(start.degree=90, clock.wise=T, gap.after=vec_gaps, 
               cell.padding = c(0.02, 0, 0.02, 0))
    
    ## draw arcs
    # order of arcs determined by vec_sectors_factor
    # length of arcs determined by mtx_sectors_xlim
    # color of arcs determined by vec_sectors_color
    circos.initialize(factors=vec_sectors_factor, xlim=mtx_sectors_xlim)
    
    circos.trackPlotRegion(ylim=c(0, 1), bg.col=vec_sectors_color, 
                           track.height=uh(3.5, "mm"), bg.border=NA)
    
    ## draw links
    # pairwise between sectors
    for (i in 1:length(vec_sectors)) {
        for (j in 1:length(vec_sectors)) {
            if (j>i) {
                #cat(i, j, "\n")
                i_sect = vec_sectors[i]
                j_sect = vec_sectors[j]
                
                i_sect_type = vec_sectors_type[i]
                j_sect_type = vec_sectors_type[j]
                
                # determine color of link/chord
                # pairwise according to connected sector types
                
                # NOTE: this approach does not support coloring 3-way overlap
                # like GC@MBC@PB with a specific color
                
                if (!color_mono) {
                    # not using a single color throughout regardless of overlap
                    
                    # overlap definition between the two current sector types
                    cur_overlap_def = paste(sort(c(i_sect_type, j_sect_type)), collapse="@")
                    
                    # if this overlap definition is specified in color_overlap
                    # use its specified color
                    if (cur_overlap_def %in% names(color_overlap)) {
                        color_link = color_overlap[cur_overlap_def]
                    } else {
                        color_link = color_rest
                    }
                } else {
                    # using a single color throughout regardless of overlap
                    color_link = color_rest
                }
                
                ### identify clones spanning both sectors
                
                # wrt df_clone_info columns
                i_sect_col_idx = which(colnames(df_clone_info)==i_sect)
                j_sect_col_idx = which(colnames(df_clone_info)==j_sect)
                
                # expect exactly 1 match
                stopifnot(length(i_sect_col_idx)==1)
                stopifnot(length(j_sect_col_idx)==1)
                
                i_bool = df_clone_info[, i_sect_col_idx]>0
                j_bool = df_clone_info[, j_sect_col_idx]>0
                
                clone_bool = i_bool & j_bool
                
                # if there's any such clone spanning both sectors
                if (any(clone_bool)) {
                    
                    cur_clones = df_clone_info[[col_clone_id]][clone_bool]
                    
                    for (cl in cur_clones) {
                        
                        # find the seqs corresponding to current clone in lst
                        # this locates the seqs along the arcs
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
                        
                        # draw the link/chord
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

