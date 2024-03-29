# courtesy of Susanna Marquez

#' Plot heatmap of one or two db features overlap
#' 
#' Plot a matrix to visualize the the number of objects (e.g. clones and sequences) shared 
#' between groups.
#'
#' @details Can be used to visualize the number of clones and/or sequences shared between 
#' compartments, which are potentially antigen-specific
#' @param db              Changeo db data.frame  
#' @param group           Vector with column names for grouping. Overlap will
#'                        be calculated across the groups. e.g \code{SAMPLE}
#' @param features        Column name of the feature column(s) shared across group. 
#'                        e.g. \code{CLONE}, \code{SEQUENCE_INPUT}
#' @param heatmap_colors  Vector of colors representing low and high values on 
#'                        the heatmap and the diagonal. Default is 
#'                        c("white","orange", "grey80")
#' @param print_zero      Show labels on zero overlap cells or not.             
#' @param long_x_angle    Angle to rotate x axis labels when any of the labels is longer
#'                        than 6 characters       
#' @param title           A string that will be used for the heatmap title. If \code{NULL},
#'                        \code{group} and \code{featureS} will be used.
#' @param xlab            Text to be used as x axis title
#' @param ylab            Text to be used as y axis title                      
#' @param plot_order      A vector to reorder the grouped columns. Can contain either 
#'                        the specific names of the grouped columns (e.g. 
#'                        \code{c("Donor1","Donor2")} ); the names 
#'                        of columns in \code{db} (e.g. \code{"DONOR"}, which
#'                        has factor levels "Donor1" and "Donor2" or \code{NULL} 
#'                        if ordering is not relevant.
#' @param silent          If T, the plot will not be printed, it will be found in the
#'                        returned list.
#' @param similarity      vector of the same length as \code{features}. For each
#'                        \code{feature}, method used to quantify the overlap. 
#'                        "min" will use the number of shared features over the 
#'                        number of features in the smaller set. "jaccard" will 
#'                        use the jaccard index, expressed as a percent, 
#'                        and defined as the intersection of features over the 
#'                        union. Can also be one of the possible values in the last grouping
#'                        column. For example to make the percent always relative to
#'                        the memory samples, use \code{group="SAMPLE", "SORT")} and 
#'                        \code{similarity=c("memory")}
#' @param na.rm           logical. If TRUE, NA values will be removed and not 
#'                        considered      
#' @param exact           Logical vector of the same length as \code{features}.
#'                        For each \code{features}, compare the exact value of 
#'                        the features (TRUE) or allow ambiguous characters 
#'                        (FALSE) using the function \code{seqEqual}. See \code{details}.
#' @return A list with the plot object and a data.frame with the values.
#' 
#' @details 
#'          
#' \strong{Note}: When using \code{exact=FALSE} and \code{distance="min"} the overlap 
#' will be calculated using the number of shared features from the smallest set. In case 
#' of ties, the smallest number of shared features.
#' 
#' @examples
#' 
#' \dontrun{
#' data("ExampleDb", package="alakazam")
#' db <- ExampleDb
#' 
#' ## Plot the number of sequences that overlap across samples
#' overlap <- plotDbOverlap(db, group="SAMPLE", features="CLONE")
#' overlap <- plotDbOverlap(db, group="SAMPLE", features=c("CLONE","JUNCTION"))
#' 
#' ## The returned plot can be modified
#' ## To edit the axis labels. the title and change the color scale and change
#' ## the theme
#' overlap$p + xlab("Sequence") + ylab("Clone") + ggtitle("New title") +
#'  scale_fill_gradient(low="white", high="orange", na.value="black") + theme_bw()
#' }
#' 
#' @export
plotDbOverlap <- function(db, group="SAMPLE_ID", features=c("CLONE","SEQUENCE_IMGT"), 
                          heatmap_colors=c("white","orange", "grey80"), 
                          print_zero=FALSE, long_x_angle=90,
                          title=NULL,xlab=NULL, ylab=NULL,
                          plot_order=NULL, silent=F, similarity=c("min", "jaccard"),
                          na.rm=FALSE, exact=c(TRUE, FALSE), 
                          geom_text_size=3, axis_text_size=10, axis_title_size=10, title_text_size=10){
    
    valid_similarities <- eval(formals(plotDbOverlap)$similarity)
    
    if (!is.data.frame(db)) { stop('Must submit a data frame') }
    
    if (!length(features) %in% c(1,2)) { stop('`features` must be of length 1 or 2')}
    if (length(exact) != length(features)) { stop('`exact` must be of the same length a `features`')}
    if (length(similarity) != length(features)) { stop('`similarity` must be of the same length a `features`')}
    if (length(heatmap_colors)<3) { stop('`heatmap_colors` must be of length 3')}
    
    # Check for valid columns
    columns <- c(group,features, plot_order)
    columns <- columns[!is.null(columns)]
    
    check <- alakazam:::checkColumns(db, columns)
    if (check != TRUE) { stop(check) }
    
    db$GROUPS <- db %>%
        dplyr::group_by_(.dots=group) %>% group_indices()
    
    group_dots <- c("GROUPS",group)
    
    group_tables <- lapply(features, function(feature) {
        
        group_table <- db %>%
            dplyr::group_by_(.dots=group_dots) %>%
            dplyr::distinct_(.dots=feature)
        
        if (na.rm) {
            group_table$NA_RM_COUNT <- as.numeric(!is.na(group_table[[feature]]))
        } else {
            group_table$NA_RM_COUNT <- 1
        }
        
        group_table <- group_table %>%
            dplyr::summarise(NUM_FEATURES=sum(NA_RM_COUNT)) %>% data.frame()
        
        group_table$GROUP_NAME <- sapply(1:nrow(group_table), function(g) {
            use <- colnames(group_table) %in% c("GROUPS","NUM_FEATURES") == F
            paste(as.matrix(group_table[g,use,drop=F]),collapse="_")
        } )
        
        NUM_GROUPS = nrow(group_table) 
        
        if (!is.null(plot_order)) {
            
            if (all(plot_order %in% group_table$GROUP_NAME)) {
                # Order by group_by vector
                # All elements in group_by must be found in group_table
                group_table <- group_table[match(plot_order, group_table$GROUP_NAME),]
            } else if (alakazam:::checkColumns(db, plot_order)) {
                # Order by column names
                # Check for valid columns with alakazam:::checkColumns
                ORDER <- bind_rows(lapply(seq_along(1:NUM_GROUPS), function(g) {
                    this_group <- group_table[g,]
                    g_id <- as.vector(unlist(this_group["GROUPS"]))
                    this_db <- db %>% dplyr::filter(GROUPS==g_id)
                    this_db <- this_db[,c(plot_order, group),drop=F]
                    idx <- do.call("order",this_db)
                    ret <- this_db[idx[1],,drop=F]
                    ret[,plot_order, drop=F]
                } ))
                group_table <- data.frame(cbind(group_table, ORDER))
                group_table <- group_table[do.call("order",group_table[,c(plot_order,group)]),]
            } else {
                ## Could be here if the filtering step removes data
                warning("Can't apply `plot_order`, using `plot_oder=NULL` instead")
                plot_order <- NULL
            }
            
        }
        
        list("description"=group_table,
             "num_groups"=NUM_GROUPS)
        
    })
    names(group_tables) <- features
    
    ## validate similarity
    tmp_tabs <- lapply(group_tables, '[[', "description")
    last_group <- rev(group)[1]
    
    for (i in 1:length(similarity)) {
        group_tables[[i]]$similarity <- similarity[i]
        if ( ! similarity[i] %in% c(valid_similarities, 
                                    as.vector(unlist(tmp_tabs[[i]][,last_group]))
        )
        ) {
            stop(paste0(similarity[i], " is not a valid similarity value."))
        }
    }
    
    # define clone overlap matrices to be plotted
    overlap_perc_matrix = matrix(0, group_tables[[1]]$num_groups, group_tables[[length(group_tables)]]$num_groups)
    dimnames(overlap_perc_matrix) = list(group_tables[[1]]$description$GROUP_NAME, 
                                         group_tables[[length(group_tables)]]$description$GROUP_NAME )
    
    overlap_note_matrix = overlap_perc_matrix;
    
    # calculate the number of overlap
    for (ii in 1:nrow(overlap_perc_matrix)) {    
        
        group_id_i <- group_tables[[1]]$description[["GROUPS"]][ii]
        group_desc_i <- group_tables[[1]]$description[ii,]
        
        for (jj in 1:ncol(overlap_perc_matrix)) {
            
            if (ii >= jj) {
                ## upper diagonal
                feature <- features[1]
                exact_mode <- exact[1]
                similarity_method <- similarity[1]
            } else {
                feature <- features[length(features)]
                exact_mode <- exact[length(exact)]
                similarity_method <- similarity[length(similarity)]
            }
            
            group_id_j <- group_tables[[length(features)]]$description[["GROUPS"]][jj]
            group_desc_j <- group_tables[[length(features)]]$description[jj,]
            
            if ( ! similarity_method %in% valid_similarities ) {
                # if same user defined value, use min
                if (group_desc_i[last_group] == group_desc_j[last_group]) {
                    warning(paste0("Switching to similarity='min':\n",
                                   group_desc_i[['GROUP_NAME']], " - " , group_desc_j[['GROUP_NAME']]
                    )
                    )
                    similarity_method <- "min"
                }
            }         
            
            # find the clones in groups i and j
            feature_ii = unique(db[[feature]][db[["GROUPS"]] == group_id_i])
            if (na.rm) {
                feature_ii <- na.omit(feature_ii)
            }
            feature_jj = unique(db[[feature]][db[["GROUPS"]] == group_id_j])
            if (na.rm) {
                feature_jj <- na.omit(feature_jj)
            }
            
            total_ii_feature = length(feature_ii)
            total_jj_feature = length(feature_jj)
            
            # find the number of features shared between groups i and j
            if (exact_mode) {
                num_shared_features <- length(intersect(feature_ii, feature_jj))
            } else {
                # exact==FALSE
                d_mat <- matrix(FALSE, nrow=total_ii_feature, ncol=total_jj_feature)
                sapply(1:total_ii_feature, function(seq_ii) {
                    sapply(1:total_jj_feature, function(seq_jj) {
                        sequence_ii <- as.vector(unlist(feature_ii[seq_ii]))
                        sequence_jj <- as.vector(unlist(feature_jj[seq_jj]))
                        if (length(sequence_ii) > 0 & length(sequence_jj) > 0) {
                            d_mat[seq_ii,seq_jj] <<- alakazam::seqEqual(sequence_ii, sequence_jj)
                        }
                    })
                })
                
                if (similarity_method == "min") {
                    ii_num_shared_features <- sum(rowSums(d_mat)>0)
                    jj_num_shared_features <- sum(colSums(d_mat)>0)
                    
                    # Return number of shared features from the smaller set
                    # if ties, 
                    # return min
                    if (total_ii_feature < total_jj_feature) {
                        num_shared_features <- ii_num_shared_features
                    } else if (total_ii_feature > total_jj_feature) {
                        num_shared_features <- jj_num_shared_features
                    } else {
                        num_shared_features <- min(ii_num_shared_features, jj_num_shared_features)
                    }
                    
                } else {
                    num_shared_features <- sum(colSums(d_mat)>0, rowSums(d_mat)>0)
                    if (any(rowSums(d_mat)==0)) {
                        feature_ii <- feature_ii[rowSums(d_mat)==0]
                    }
                    if (any(colSums(d_mat)==0)) {
                        feature_jj <- feature_jj[colSums(d_mat)==0]
                    }
                }
            } # end exact==FALSE
            
            # Get %, different denominators for min, jaccard and user defined
            if (similarity_method == "min") {
                total_feature = min(total_ii_feature, total_jj_feature)
                # find the percentage of feature shared relative to the total number of feature
                perc_shared_feature = round ( num_shared_features / total_feature * 100 , 1)
            } else if (similarity_method=="jaccard" ) {
                ## Intersection over union
                perc_shared_feature <- round(num_shared_features/length(unique(union(feature_ii, feature_jj))) * 100, 1)
            } else { # user defined
                if ( group_desc_i[[last_group]] == similarity_method ) {
                    perc_shared_feature = round ( num_shared_features / total_ii_feature * 100 , 1)
                } else {
                    perc_shared_feature = round ( num_shared_features / total_jj_feature * 100 , 1)
                }
            }
            overlap_perc_matrix[ii, jj] <- perc_shared_feature
            
            # define the text to be displayed on the heatmap
            if (ii == jj) {
                overlap_note_matrix[ii, jj] = num_shared_features
            } else {
                overlap_note_matrix[ii, jj] = paste0 (num_shared_features, "\n(", perc_shared_feature, "%)")
            }
            if (num_shared_features == 0){
                if (print_zero) {
                    overlap_note_matrix[ii, jj] = 0
                } else {
                    overlap_note_matrix[ii, jj] = ""
                }
            }
            
        }
    }   
    
    # the diagnoal values are not used for coloring
    diag(overlap_perc_matrix) = NA
    diag(overlap_note_matrix) = NA
    
    # organize the data into ggplot format
    melt_overlap_perc <- reshape2::melt(overlap_perc_matrix, 
                                        list("Var1","Var2"), value.name="value")
    melt_overlap_note <- reshape2::melt(overlap_note_matrix,
                                        list("Var1","Var2"), value.name="value")
    
    # Get default plot_order, based on group_table order
    factor.order <- as.vector(unlist(group_tables[[1]]$description$GROUP_NAME))
    
    melt_overlap_perc$Var1 = factor(melt_overlap_perc$Var1, levels=factor.order)
    melt_overlap_perc$Var2 = factor(melt_overlap_perc$Var2, levels=factor.order)
    melt_overlap_note$Var1 = factor(melt_overlap_note$Var1, levels=factor.order)
    melt_overlap_note$Var2 = factor(melt_overlap_note$Var2, levels=factor.order) 
    
    .e = environment()
    
    ## Will rotate x-axis labels if long names found
    my_angle <- 0
    if (any(sapply(factor.order,nchar)>6)) {
        my_angle <- long_x_angle
    }
    
    if (is.null(xlab)) {
        xlab <- features[1]
    }
    
    if (is.null(ylab)) {
        ylab <- features[length(features)]
    }
    
    # plot the heatmap
    p <- ggplot(melt_overlap_perc, aes(Var1, Var2, group=Var2), environment = .e) +
        geom_tile(aes(fill = value)) + 
        theme_bw() +
        geom_text(aes(label = melt_overlap_note$value), na.rm = TRUE, size=geom_text_size) +
        geom_text(data=group_tables[[1]]$description %>% dplyr::mutate("Var1"=GROUP_NAME, "Var2"=GROUP_NAME),
                  aes(label = NUM_FEATURES), angle = 45, color="white", vjust=1.5, hjust=0.5, na.rm = TRUE, size=geom_text_size) +
        geom_text(data=group_tables[[length(features)]]$description %>% dplyr::mutate("Var1"=GROUP_NAME, "Var2"=GROUP_NAME),
                  aes(label = NUM_FEATURES), angle = 45, color="white", vjust=-0.5, hjust=0.5, na.rm = TRUE, size=geom_text_size) +     
        #annotate("rect", xmin = 0.5, xmax = 3.5, ymin = 0.5, ymax = 3.5, alpha = 0, colour = "red") + #* hard-coded for spencer project
        #annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 3.5, ymax = 6.5, alpha = 0, colour = "red") + #* hard-coded for spencer project
        #annotate("rect", xmin = 6.5, xmax = 9.5, ymin = 6.5, ymax = 9.5, alpha = 0, colour = "red") + #* hard-coded for spencer project
        annotate("segment", x = 0.5, xend = nrow(group_tables[[1]]$description)+0.5, 
                 y = 0.5, yend = nrow(group_tables[[1]]$description)+0.5, colour = "white") + 
        scale_fill_gradient(low = heatmap_colors[1], 
                            high = heatmap_colors[2],
                            na.value = heatmap_colors[3] ) +
        xlab(xlab) + ylab(ylab) +
        theme(legend.position = "none", 
              #axis.title.x = element_blank(), 
              #axis.title.y = element_blank(),
              axis.title = element_text(size=axis_title_size),
              axis.text.y = element_text(size=axis_text_size),
              axis.text.x = element_text(angle = my_angle, hjust = 1, size=axis_text_size)) +
        scale_x_discrete (breaks=levels(melt_overlap_perc$Var1),
                         labels=levels(melt_overlap_perc$Var2)) +
        scale_y_discrete(breaks=levels(melt_overlap_perc$Var2),
                         labels=levels(melt_overlap_perc$Var2))
    if (!is.null(title)) {
        if (title == "auto") {
            title <- paste0("Unique ", paste(features, collapse=" and ")," overlap across ", paste(group,collapse=", "))
            for (i in 1:length(features)) {
                new_line <- paste0("\n",
                                   features[i],": ",
                                   paste0("similarity=", similarity[i]),
                                   ", ",
                                   paste0("exact=", exact[i])
                )
                title <- paste(title, new_line, sep="")
            }
            
        }
    }
    p <- p + ggtitle(title) + 
        theme(plot.title = element_text(lineheight=.8, face="bold", size=title_text_size))
    if (!silent) print(p)
    invisible(list("p"=p,
                   "perc"=melt_overlap_perc,
                   "note"=melt_overlap_note)
    )
}
