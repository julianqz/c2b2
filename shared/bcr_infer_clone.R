# infer B cell clones 

#' Define clones via hierarchical clustering
#'
#' Given a single partition of heavy chain-only bulk-seqs based on heavy chain VJL combination, 
#' or of VH:VL paired single-cell BCRs based on heavy and light chain VJL combinations,
#' perform clonal clustering with a specified distance threshold. 
#' (VJL=V gene, J gene, junction length)
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  name of the column containing heavy chain junctions. It is the distance
#'                           between these heavy chain junction sequences that is used for hierarchical
#'                           clustering. Note that light chain sequences are not used at this stage.
#' @param    VJLgroupColumn  name of the column containing information about partitioning of heavy chain
#'                           bulk-seq based on heavy chain VJL combination, or that of single-cell VH:VL 
#'                           paired BCRs based on heavy and light chain VJL combinations. 
#'                           Can be obtained by calling \link{distToNearest} with \code{keepVJLgroup=TRUE}
#'                           or \link[alakazam]{groupGenes}.
#' @param    cloneColumn     name of the column to be generated that contains clone IDs. 
#'                           Defaults to \code{"CLONE"}.
#' @param    threshold       a numeric value between 0 and 1 to be used as the clonal clustering
#'                           threshold 
#' @param    linkage         a character value specifying the type of linkage to be used for hierarchical
#'                           clustering. Any value accepted by \link[stats]{hclust} is fine. 
#'                           Defaults to \code{"single"}.
#' @param    maxmiss         The maximum number of non-ACGT characters (gaps or Ns) to permit in the 
#'                           junction sequence before excluding the record from clonal assignment. Note, 
#'                           under single linkage non-informative positions can create artifactual links 
#'                           between unrelated sequences. Use with caution. Default to \code{0}.    
#' @param    verbose         if \code{TRUE}, print a message when a sequence is excluded due to the presence
#'                           of more than (strictly \code{>}) \code{maxmiss} number of non-ACGT characters.                                                                   
#' 
#' @return   Returns a list containing two components. (i) \code{CLUSTERED}: a modified \code{db} data.frame 
#'           with a new column (named after \code{cloneColumn}) containing clone IDs as a result of hierarchical 
#'           clustering given the threshold specified by \code{threshold}. (ii) \code{EXCLUDED}: any sequence
#'           excluded from clustering for containing more than (strict \code{>}) \code{maxmiss} number of 
#'           non-ACGT characters in \code{sequenceColumn} will be stored here. \code{NULL} if no sequence is
#'           excluded.
#'
#' @details
#' 
#' The underlying operation can be expected to be the same as \code{DefineClones.py} as in 
#' \code{Change-O} with the following specifications:
#' \itemize{
#'    \item \code{--mode gene}
#'    \item \code{--act set}
#'    \item \code{--model ham}
#'    \item \code{--norm len}
#'    \item \code{--link single}
#' }
#' 
#' The advantage is that with this function in R, partitioning based on heavy chain VJL combination,
#' or heavy and light chain VJL combinations in single-cell mode, needs to be performed only once, 
#' provided that the partitioning column has been saved via, for instance, \code{keepVJLgroup} in
#' \link{distToNearest}.

defineClonePerPartition = function(db, sequenceColumn="junction", VJLgroupColumn="vjl_group", cloneColumn="clone_id", 
                                   threshold, maxmiss=0, verbose=F,
                                   linkage=c("single", "complete", "average", "ward.D", "ward.D2", 
                                             "mcquitty", "median", "centroid")) {
    require(stringi)
    require(alakazam)
    
    # all rows must be in same partition
    curGrp = unique(db[[VJLgroupColumn]])
    stopifnot( length(curGrp) == 1 )
    
    # all rows have same junciton length
    curJuncs <- db[[sequenceColumn]]
    stopifnot( length(unique(stri_length(curJuncs))) == 1 )
    curJuncsLen = stri_length(curJuncs[1])
    
    # return object
    returnLst <- vector(mode="list", length=2)
    names(returnLst) <- c("CLUSTERED", "EXCLUDED")
    
    # filter based on number of non-ATGC characters
    # changeo::filterMissing() uses <=
    miss_count <- stri_count_regex(str=db[[sequenceColumn]], pattern="[^ATGC]")
    # if TRUE, keep; if FALSE, exclude
    bool_keep <- miss_count <= maxmiss
    
    if (any(!bool_keep) & verbose) {
        cat("Excluded due to non-ATGC character(s) in", sequenceColumn, ":", sum(!bool_keep), "\n")    
    }
    
    # if nothing excluded, remains NULL
    if (sum(!bool_keep)>0) {
        returnLst[["EXCLUDED"]] <- db[!bool_keep, ]
        
        db <- db[bool_keep, ]
        curJuncs = curJuncs[bool_keep]
    }
    
    if (nrow(db)>0) {
        
        # only 1 seq: no need to cluster
        if (nrow(db)==1) {
            curIDs <- paste0(curGrp, "_1")
        } else {
            # pairwise hamming distance (square matrix)
            curHam <- pairwiseDist(seq=curJuncs) / curJuncsLen[1]
            # convert to dist class (1d condensed distance vector)
            # important not to use dist()! 
            curDist <- as.dist(m=curHam)
            # hclust using threshold
            # calls stats::hclust
            # faster altnertive: `fastcluster` package
            curClust <- fastcluster::hclust(d=curDist, method=linkage)
            #plot(curClust)
            #abline(h=threshold, lty=2, col=2)
            
            # cut
            curCut <- cutree(curClust, h=threshold)
            #rect.hclust(curClust, h=threshold)
            
            # record clone assignment
            curIDs <- paste(curGrp, curCut, sep="_")
        }
        
        stopifnot(length(curIDs) == nrow(db))
        stopifnot(!any(is.na(curIDs)))
        
        # column for clone ID
        db[[cloneColumn]] <- curIDs
        
        returnLst[["CLUSTERED"]] = db
        
    } else {
        if (verbose) { cat("Nothing left to cluster on; skipped\n") }
        # returnLst$CLUSTERED remains NULL
    }
    
    return(returnLst)
}

#TODO: add parallelization

#' Wrapper to infer B cell clones for an individual
#' 
#' Calls `defineClonePerPartition` on a `db` containing all `VJLgroupColumn` 
#' partitions belonging to an individual. 
#' 
#' See doc for `defineClonePerPartition` for parameters. 
#' The only difference is that `db` here contains all `VJLgroupColumn` partitions 
#' instead of a single `VJLgroupColumn` partition.
#' 
#' @return  A list containing `db_clust` and `db_fail`, which can be either 
#'          `data.frame` or `NULL`. 

defineCloneDb = function(db, sequenceColumn="junction", VJLgroupColumn="vjl_group", 
                         cloneColumn="clone_id", 
                         threshold, maxmiss=0, linkage="single", verbose=T) {
    # checks
    stopifnot(is.numeric(threshold) && (threshold>=0 & threshold<=1))
    
    stopifnot(all(c(sequenceColumn, VJLgroupColumn) %in% colnames(db)))
    if (cloneColumn %in% colnames(db)) {
        warning(cloneColumn, " already exists and will be overwritten.\n")
    }
    
    # loop through partitions
    
    uniqGrps = unique(db[[VJLgroupColumn]])
    
    clustPassLst = vector(mode="list", length=length(uniqGrps))
    names(clustPassLst) = uniqGrps
    
    clustFailLst = vector(mode="list", length=length(uniqGrps))
    names(clustFailLst) = uniqGrps
    
    for (i_grp in 1:length(uniqGrps)) {
        
        if (verbose) { if (i_grp%%100==0) { cat(i_grp, "\n") } }
        
        grp = uniqGrps[i_grp]
        
        # wrt db
        grp_idx = which(db[[VJLgroupColumn]]==grp)
        stopifnot(length(grp_idx)>0)
        
        # hierarchical clustering
        # returns a list with $CLUSTERED and $EXCLUDED
        curClust = defineClonePerPartition(db=db[grp_idx, ], 
                                           sequenceColumn=sequenceColumn, 
                                           VJLgroupColumn=VJLgroupColumn, 
                                           cloneColumn=cloneColumn, 
                                           maxmiss=maxmiss, linkage=linkage,
                                           threshold=threshold,
                                           verbose=F)
        
        if (!is.null(curClust[["CLUSTERED"]])) {
            clustPassLst[[grp]] = curClust[["CLUSTERED"]]
        }
        if (!is.null(curClust[["EXCLUDED"]])) {
            clustFailLst[[grp]] = curClust[["EXCLUDED"]]
        }
        
        rm(curClust, grp, grp_idx)
    }
    
    # if all seqs passed, clustFailLst will contain all NULL's
    # in which case db_fail will also be NULL
    db_fail = do.call(rbind, clustFailLst)
    if (verbose) { cat("Number of seqs failed:", 
                       ifelse(is.null(db_fail), 0, nrow(db_fail)), "\n") }
    
    # if all seqs failed, clustPassLst will contain all NULL's
    # in which case db_clust will also be NULL
    db_clust = do.call(rbind, clustPassLst)
    if (verbose) { cat("Number of seqs passed:", 
                       ifelse(is.null(db_clust), 0, nrow(db_clust)), "\n") }
    
    return(list(db_clust=db_clust, db_fail=db_fail))
}

