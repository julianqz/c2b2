# Julian Q. Zhou
# https://github.com/julianqz

#' Create BCR network for B cell clones using igraph
#' 
#' @param  db         A `data.frame` containing BCR data.
#' @param  col_clone  Column name containing clone IDs.
#' @param  col_seq    Column name for sequences used for computing edges.
#' @param  threshold  Hamming distance threshold used for computing edges.
#' @param  nt_or_aa   Whether to use nucleotide or aa distances. Either `nt` or `aa`.
#' @param  bool_norm  Whether to normalize distance by length of sequence. Boolean.
#' @param  bool_mst   Whether to convert graph into minimum spanning tree. Boolean.
#' 
#' @return  An igraph graph object.
#' 
#' @details This function solely handles graph building and does not 
#'          handle aesthetics.
#' 
get_bcr_network = function(db, col_clone, col_seq, 
                           threshold, nt_or_aa=c("nt", "aa"), 
                           bool_norm=TRUE, bool_mst=TRUE) {
    
    suppressPackageStartupMessages(require(igraph))
    suppressPackageStartupMessages(require(alakazam))
    
    # related functions
    # ?graph_from_adj_list
    # ?make_empty_graph
    # ?add_edges
    # ?vertex_attr_names
    
    # checks
    stopifnot(nrow(db)>0)
    stopifnot(all( c(col_clone, col_seq) %in% colnames(db) ))
    stopifnot(threshold>=0)
    stopifnot(nt_or_aa %in% c("nt", "aa"))
    if (bool_norm) { stopifnot(threshold<=1) }
    
    # initialize an undirected graph with nrow(db) vertices with no edges
    g = make_empty_graph(n=nrow(db), directed=F) 
    
    # add edges
    # loop thru clones
    # draw edge if (normalized) hamming distance (nt/aa) btw col_seq (eg. CDR3)
    # is within (<=) threshold
    
    cat("Calculating edges based on", 
        ifelse(bool_norm, "normalized", ""),
        nt_or_aa, "Hamming distance\nbetween `", col_seq, 
        "` of clonal members, using a threshold of", threshold, "...\n")
    
    for (cl in unique(db[[col_clone]])) {
        # wrt db
        # these are the same integer IDs of vertices
        idx = which(db[[col_clone]]==cl)
        
        # if more than 1 clonal member
        if (length(idx)>1) {
            
            # distance matrix
            if (nt_or_aa=="nt") {
                mtx_ham = pairwiseDist(db[[col_seq]][idx],
                                       dist_mat=alakazam::getDNAMatrix())
            } else {
                mtx_ham = pairwiseDist(db[[col_seq]][idx],
                                       dist_mat=alakazam::getAAMatrix())
            }
            
            # normalize by length?
            if (bool_norm) {
                mtx_ham = mtx_ham / nchar(db[[col_seq]][idx][1])
            }
            
            # binarize
            mtx_ham_bin = mtx_ham <= threshold
            
            # initialize list holder for edges
            lst_edges = vector(mode="list", 
                               length=sum(mtx_ham_bin[upper.tri(mtx_ham_bin, diag=F)]))
            i_edge = 1
            
            for (i_row in 1:length(idx)) {
                for (i_col in 1:length(idx)) {
                    if (i_col > i_row) {
                        #cat(i_row, i_col, "\n")
                        if (mtx_ham_bin[i_row, i_col]) {
                            lst_edges[[i_edge]] = c(idx[i_row], idx[i_col])
                            i_edge = i_edge + 1
                        }
                    }
                }
            }
            # +1 because after adding last edge, i_edge still undergoes +1
            stopifnot( i_edge == length(lst_edges)+1 )
            
            # add edges to g
            g = add_edges(graph=g, edges=unlist(lst_edges))
            
            rm(lst_edges, mtx_ham, mtx_ham_bin)
        }
        rm(idx)
    }
    
    # minimum spanning tree?
    if (bool_mst) {
        cat("Converting graph to minimum spanning tree...\n")
        g = mst(graph=g)
    }
    
    cat("Done.\n")
    
    # sanity check
    stopifnot( length(V(g)) == nrow(db) )
    
    return(g)
}


#' adapted from ?igraph::shapes (changed marked by #*)
#' generic star vertex shape, with a parameter for number of rays
#' 
#' @exapmles
#' # no clipping, edges will be below the vertices anyway
#' add_shape("star", clip=shape_noclip,
#'           plot=mystar, 
#'           parameters=list(vertex.norays=5, vertex.frame.color="transparent"))
#' g <- make_ring(length(shapes))
#' plot(g, vertex.shape="star", vertex.color=rainbow(vcount(g)),
#'      vertex.size=seq(10,20,length=vcount(g)))
#' plot(g, vertex.shape="star", vertex.color=rainbow(vcount(g)),
#'      vertex.size=seq(10,20,length=vcount(g)),
#'      vertex.norays=rep(4:8, length=vcount(g)))
#'      
mystar <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
    }
    vertex.size  <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
    }
    norays <- params("vertex", "norays")
    if (length(norays) != 1 && !is.null(v)) {
        norays <- norays[v]
    }
    #* added by JZ
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
        vertex.frame.color <- vertex.frame.color[v]
    }
    #* JZ added vertex.frame.color and fg
    mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays, vertex.frame.color,
           FUN=function(x, y, bg, size, nor, fg) {
               symbols(x=x, y=y, bg=bg, fg=fg,
                       stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                       add=TRUE, inches=FALSE)
           })
}
