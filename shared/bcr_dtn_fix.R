# copied from shazam/R/DistToNearest.R, v1.0.2
# contains bug fix for btw-subj `distToNearest` pushed in commit 47bffe0 on 2020-11-10
# bug fix not in stable release v1.0.2
# changes are indicated by #* (main change on L212-213)

distToNearest_fix <- function(db, sequenceColumn="junction", vCallColumn="v_call", jCallColumn="j_call", 
                          model=c("ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "mk_rs5nf", "m1n_compat", "hs1f_compat"), 
                          normalize=c("len", "none"), symmetry=c("avg", "min"),
                          first=TRUE, VJthenLen=TRUE, nproc=1, fields=NULL, cross=NULL, 
                          mst=FALSE, subsample=NULL, progress=FALSE,
                          cellIdColumn=NULL, locusColumn="locus", 
                          onlyHeavy=TRUE, keepVJLgroup=TRUE) {
    #* 
    require(foreach)
    
    # Hack for visibility of foreach index variables
    i <- NULL
    
    # Initial checks
    model <- match.arg(model)
    normalize <- match.arg(normalize)
    symmetry <- match.arg(symmetry)
    if (!is.data.frame(db)) { stop('Must submit a data frame') }
    
    # Check base input
    check <- alakazam::checkColumns(db, c(sequenceColumn, vCallColumn, jCallColumn, fields, cross))
    if (check != TRUE) { stop(check) }
    
    # Check single-cell input
    if (!is.null(cellIdColumn)) {
        check <- alakazam::checkColumns(db, c(cellIdColumn, locusColumn))
        if (check != TRUE) { stop(check) }
    }
    
    # Cast all columns to character
    columns <- c(sequenceColumn, vCallColumn, jCallColumn, fields, cross,
                 cellIdColumn, locusColumn)
    columns <- columns[!is.null(columns) & columns %in% names(db)]
    for (cl in columns) { db[[cl]] <- as.character(db[[cl]]) }
    
    # Convert sequence columns to uppercase
    db <- shazam:::toupperColumns(db, c(sequenceColumn)) 
    
    # Single-cell mode?
    if (!is.null(cellIdColumn) & !is.null(locusColumn)) {
        singleCell <- TRUE
        
        # check locus column
        valid_loci <- c("IGH", "IGI", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
        check <- !all(unique(db[[locusColumn]]) %in% valid_loci)
        if (check) {
            stop("The locus column contains invalid loci annotations.")
        }
    } else {
        singleCell <- FALSE
    } 
    
    # Disallow multiple heavy chains per cell
    if (singleCell) {
        # check multiple heavy chains
        x <- sum(table(db[[cellIdColumn]][db[[locusColumn]] == "IGH"]) > 1)
        if (x > 0) {
            stop(paste(x, "cell(s) with multiple heavy chains found. One heavy chain per cell is expected."))
        }
        # check multiple beta chains
        x <- sum(table(db[[cellIdColumn]][db[[locusColumn]] == "TRB"]) > 1)
        if (x > 0) {
            stop(paste(x, "cell(s) with multiple beta chains found. One beta chain per cell is expected."))
        }
        # check multiple delta chains
        x <- sum(table(db[[cellIdColumn]][db[[locusColumn]] == "TRD"]) > 1)
        if (x > 0) {
            stop(paste(x, "cell(s) with multiple delta chains found. One delta chain per cell is expected."))
        }
    }
    
    # Check for invalid characters
    # heavy
    valid_seq <- sapply(db[[sequenceColumn]], shazam:::allValidChars, shazam:::getCharsInModel(model)) 
    not_valid_seq <- which(!valid_seq)
    if (length(not_valid_seq) > 0) {
        warning("Invalid sequence characters in the ", sequenceColumn, 
                " column. ", length(not_valid_seq), " sequence(s) removed")
        db <- db[valid_seq, ]
    }
    
    # junction length columns (prep for groupGenes)
    junc_len <- "JUNC_LEN"
    db[[junc_len]] <- stringi::stri_length(db[[sequenceColumn]])
    
    # create V+J grouping, or V+J+L grouping
    if (VJthenLen) {
        # 2-stage partitioning using first V+J and then L
        # V+J only first
        # creates $vj_group
        db <- alakazam::groupGenes(db, v_call=vCallColumn, j_call=jCallColumn, junc_len=NULL,
                         cell_id=cellIdColumn, locus=locusColumn, only_heavy=onlyHeavy,
                         first=first)
        # L (later)  
        group_cols <- c("vj_group", junc_len)
        
    } else {
        # 1-stage partitioning using V+J+L simultaneously
        # creates $vj_group
        # note that despite the name (VJ), this is based on V+J+L
        db <- alakazam::groupGenes(db, v_call=vCallColumn, j_call=jCallColumn, junc_len=junc_len,
                         cell_id=cellIdColumn, locus=locusColumn, only_heavy=onlyHeavy,
                         first=first)
        group_cols <- c("vj_group")
    }
    
    # groups to use
    if (!is.null(fields)) {
        group_cols <- append(group_cols,fields)
    }
    # unique groups
    # not necessary but good practice to force as df and assign colnames
    # (in case group_cols has length 1; which can happen in groupBaseline)
    uniqueGroups <- data.frame(unique(db[, group_cols]), stringsAsFactors=FALSE)
    colnames(uniqueGroups) <- group_cols
    rownames(uniqueGroups) <- NULL
    # indices
    # crucial to have simplify=FALSE 
    # (otherwise won't return a list if uniqueClones has length 1)
    uniqueGroupsIdx <- sapply(1:nrow(uniqueGroups), function(i){
        curGroup <- data.frame(uniqueGroups[i, ], stringsAsFactors=FALSE)
        colnames(curGroup) <- group_cols
        # match for each field
        curIdx <- sapply(group_cols, function(coln){
            db[[coln]]==curGroup[, coln]
        }, simplify=FALSE)
        curIdx <- do.call(rbind, curIdx)
        # intersect to get match across fields 
        curIdx <- which(colSums(curIdx)==length(group_cols))
        # sanity check
        # no NA
        stopifnot( all(!is.na(curIdx)) )
        # index within range of db
        stopifnot( max(curIdx) <= nrow(db) )
        return(curIdx)
    }, simplify=FALSE)
    
    # Create new column for distance to nearest neighbor
    db$TMP_DIST_NEAREST <- rep(NA, nrow(db))
    db$ROW_ID <- 1:nrow(db)
    
    # Create cluster of nproc size and export namespaces
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.
    if( nproc==1 ) {
        # If needed to run on a single core/cpu then, register DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        foreach::registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        doParallel::registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }
    
    # Export groups to the clusters
    if (nproc > 1) { 
        export_functions <- list("db",
                                 "uniqueGroupsIdx", 
                                 "cross",
                                 "mst",
                                 "subsample",
                                 "sequenceColumn", 
                                 "model",
                                 "normalize",
                                 "symmetry",
                                 "shazam:::nearestDist", 
                                 "shazam:::HH_S1F_Distance",
                                 "shazam:::MK_RS1NF_Distance",
                                 "shazam:::HH_S5F_Distance",
                                 "shazam:::MK_RS5NF_Distance",
                                 "shazam:::HS1F_Compat",
                                 "shazam:::M1N_Compat",
                                 "shazam:::calcTargetingDistance",
                                 "shazam:::findUniqSeq",
                                 "shazam:::pairwise5MerDist",
                                 "shazam:::nonsquare5MerDist",
                                 "singleCell",
                                 "locusColumn")
        parallel::clusterExport(cluster, export_functions, envir=environment())
    }
    
    
    
    n_groups <- length(uniqueGroupsIdx)
    if (progress) { 
        pb <- alakazam::progressBar(n_groups) 
    }
    tryCatch(list_db <- foreach(i=1:n_groups, .errorhandling='stop') %dopar% {
        # wrt db
        idx <- uniqueGroupsIdx[[i]]
        
        if (singleCell) {
            # only use IGH
            # wrt idx
            idxBool <- db[[locusColumn]][idx] == "IGH"
        } else {
            idxBool <- rep(TRUE, length(idx))
        }
        
        db_group <- db[idx, ]
        
        crossGroups <- NULL
        if (!is.null(cross)) {
            #*
            x <- dplyr::group_by(db_group, !!!rlang::syms(cross))
            crossGroups <- dplyr::group_indices(x)
        }
        
        arrSeqs <-  db[[sequenceColumn]][idx]
        
        db_group$TMP_DIST_NEAREST[idxBool] <- shazam:::nearestDist(arrSeqs[idxBool], 
                                                          model=model,
                                                          normalize=normalize,
                                                          symmetry=symmetry,
                                                          crossGroups=crossGroups[idxBool],
                                                          mst=mst,
                                                          subsample=subsample)
        # Update progress
        if (progress) { pb$tick() }
        
        return(db_group)
    }, 
    error = function(e) {
        if (nproc > 1 & grepl("Error in unserialize(socklist[[n]]) : error reading from connection", e, fixed=TRUE)) {
            warning("There is an error running the code in parallel. Try with nproc=1.")
        }
        stop(e)
    }
    )
    
    # Convert list from foreach into a db data.frame
    db <- do.call(rbind, list_db)
    db <- db[order(db$ROW_ID), ]
    
    # Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    if (!is.null(cross)) {
        db$cross_dist_nearest <- db$TMP_DIST_NEAREST
    } else {
        db$dist_nearest <- db$TMP_DIST_NEAREST
    }
    
    # prepare db for return
    if ((!VJthenLen) && keepVJLgroup) {
        db$vjl_group <- db[["vj_group"]]
    }
    db <- db[, !(names(db) %in% c(junc_len, "vj_group", "ROW_ID", "V1", "J1","TMP_DIST_NEAREST"))]
    
    return(db)
}
