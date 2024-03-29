# Julian Q. Zhou
# https://github.com/julianqz

#' Tabulate V-J gene usage
#' 
#' @param   vec_v          Character vector containing V gene/allele annotations.
#'                         Vector length should match that of `vec_j` and `vec_subj`. 
#' @param   vec_j          Character vector containing J gene/allele annotations.
#'                         Vector length should match that of `vec_v` and `vec_subj`.
#' @param   vec_subj       Character vector containing subject labels.
#'                         Vector length should match that of `vec_v` and `vec_j`.
#' @param   anno_separator Character delimiting ambiguous annotations. Passed on
#'                         to `sep` of `alakazam::getGene`.
#' 
#' @return  A nested list containing `common` and `subj_specific`, each of which
#'          is itself a list of `data.frame`. 
#'          
#'          The number of `data.frame`s matches the number of unique subjects 
#'          in `vec_subj`. Each `data.frame` contains columns `subj`, `v`, `j`, 
#'          and `count`.
#'          
#'          In `common`, each `data.frame` contains the same number of rows 
#'          (i.e. same total set of V-J combinations) for all subjects. In 
#'          other words, the V-J combinations being counted in `common` are the
#'          union of V-J combinations found in each subject. 
#'          
#'          In `subj_specific`, each `data.frame` contains a number of rows 
#'          specific to the subject being represented. 
#'          
#'          Both `common` and `subj_specific` are returned so as to provide  
#'          convenience when visualizing -- `common` allows heatmaps with common
#'          rows and columns to be drawn for all subjects, whereas `subj_specific`
#'          allows heatmaps with unique rows and columns to be drawn for each subject.
#' 
#' @details Tip: if data for different subjects are stored in separate `data.frame`s, 
#'          concatenate to create the required input vectors.
#' 
#' @example tabulate_vj(vec_v=db[["germline_v_call"]], 
#'                      vec_j=db[["germline_j_call"]],
#'                      vec_subj=db[["donor"]], ",")         
#'                      
tabulate_vj = function(vec_v, vec_j, vec_subj, 
                       anno_separator=",") {
    
    # getGene
    require(alakazam) 
    
    stopifnot(length(vec_v)==length(vec_j))
    stopifnot(length(vec_v)==length(vec_subj))
    
    # unique entries
    vec_v_uniq = unique(vec_v)
    vec_j_uniq = unique(vec_j)
    vec_subj_uniq = unique(vec_subj)
    
    # extract v/j gene for unique entries
    # note that vec_[vj]_uniq_getGene could still contain duplicates
    
    # important: do not set strip_d to True (which is the default)
    # getGene assumes any allele with a "D" in name is a duplicate, 
    #   which is decisively not the case 
    #   (e.g. IGHV3-43 != iGHV3-43D nucleotide-wise)
    
    vec_v_uniq_getGene = getGene(segment_call=vec_v_uniq,
                                 first=F, collapse=T, strip_d=F, omit_nl=F, 
                                 sep=anno_separator)
    vec_j_uniq_getGene = getGene(segment_call=vec_j_uniq,
                                 first=F, collapse=T, strip_d=F, omit_nl=F, 
                                 sep=anno_separator)
    
    # map back v/j gene to all entries
    # wrt vec_[vj]_uniq
    idx_v_uniq = match(vec_v, vec_v_uniq)
    idx_j_uniq = match(vec_j, vec_j_uniq)
    stopifnot(all.equal(vec_v, vec_v_uniq[idx_v_uniq]))
    stopifnot(all.equal(vec_j, vec_j_uniq[idx_j_uniq]))
    
    vec_v_gene = vec_v_uniq_getGene[idx_v_uniq]
    vec_j_gene = vec_j_uniq_getGene[idx_j_uniq]
    
    # unique v/j genes, across subjects
    vec_v_gene_uniq = sort(unique(vec_v_gene))
    vec_j_gene_uniq = sort(unique(vec_j_gene))
    
    # tabulate for each subject
    # dimension: [v, j, subj]
    tab = table(vec_v_gene, vec_j_gene, vec_subj)
    stopifnot(all( dim(tab) == c(length(vec_v_gene_uniq), 
                                 length(vec_j_gene_uniq),
                                 length(vec_subj_uniq)) ))
    
    # list containing df common to all subjects
    lst_df_common = vector(mode="list", length=length(vec_subj_uniq))
    names(lst_df_common) = vec_subj_uniq
    
    # list containing df unique to each subject
    lst_df_subj = vector(mode="list", length=length(vec_subj_uniq))
    names(lst_df_subj) = vec_subj_uniq
    
    # base df common to all subjects
    base_df_cols = c("subj", "v", "j", "count")
    base_df = data.frame(matrix(NA, 
                                nrow=length(vec_v_gene_uniq)*length(vec_j_gene_uniq),
                                ncol=length(base_df_cols)))
    colnames(base_df) = base_df_cols
    # as.vector linearizes matrix by column, not by row
    base_df[["v"]] = rep(vec_v_gene_uniq, times=length(vec_j_gene_uniq))
    base_df[["j"]] = rep(vec_j_gene_uniq, each=length(vec_v_gene_uniq))
    
    for (i_subj in 1:length(vec_subj_uniq)) {
        i_subj_df = base_df
        i_subj_df[["subj"]] = vec_subj_uniq[i_subj]
        i_subj_tab_common = tab[,,i_subj]
        i_subj_df[["count"]] = as.vector(i_subj_tab_common)
        
        # sanity check
        # no NA
        stopifnot(!any(is.na(i_subj_df)))
        # total counts per v gene
        stopifnot(all.equal( rowSums(i_subj_tab_common),
                             sapply(vec_v_gene_uniq, function(g){
                                 sum(i_subj_df[["count"]][i_subj_df[["v"]]==g])
                             }) ))
        # total counts per j gene
        stopifnot(all.equal( colSums(i_subj_tab_common),
                             sapply(vec_j_gene_uniq, function(g){
                                 sum(i_subj_df[["count"]][i_subj_df[["j"]]==g])
                             }) ))
        
        lst_df_common[[i_subj]] = i_subj_df
        
        # df specific to subject
        
        # remove from tab rows and cols containing only 0
        i_subj_bool_row = rowSums(i_subj_tab_common)==0
        i_subj_bool_col = colSums(i_subj_tab_common)==0
        # can't be the case that all rows contain 0 only
        stopifnot(!all(i_subj_bool_row))
        # can't be the case that all cols contain 0 only
        stopifnot(!all(i_subj_bool_col))
        i_subj_tab = i_subj_tab_common[!i_subj_bool_row, !i_subj_bool_col]
        
        # make another df based on new tab
        i_subj_df_2 = data.frame(matrix(NA, 
                                      nrow=prod(dim(i_subj_tab)),
                                      ncol=length(base_df_cols)))
        colnames(i_subj_df_2) = base_df_cols
        i_subj_df_2[["subj"]] = vec_subj_uniq[i_subj]
        # as.vector linearizes matrix by column, not by row
        i_subj_df_2[["v"]] = rep(rownames(i_subj_tab), times=ncol(i_subj_tab))
        i_subj_df_2[["j"]] = rep(colnames(i_subj_tab), each=nrow(i_subj_tab))
        i_subj_df_2[["count"]] = as.vector(i_subj_tab)
        
        # sanity check
        # no NA
        stopifnot(!any(is.na(i_subj_df_2)))
        # total counts per v gene
        stopifnot(all.equal( rowSums(i_subj_tab),
                             sapply(rownames(i_subj_tab), function(g){
                                 sum(i_subj_df_2[["count"]][i_subj_df_2[["v"]]==g])
                             }) ))
        # total counts per j gene
        stopifnot(all.equal( colSums(i_subj_tab),
                             sapply(colnames(i_subj_tab), function(g){
                                 sum(i_subj_df_2[["count"]][i_subj_df_2[["j"]]==g])
                             }) ))
        
        # two df's should sum up to the same
        stopifnot(sum(i_subj_df[["count"]]) == sum(i_subj_df_2[["count"]]))
        
        # # factorize?
        # # v/j for subject-specific df need to be done here
        # if (factorize) {
        #     i_subj_df_2[["v"]] = factor(x=i_subj_df_2[["v"]], 
        #                                 levels=sort(unique(i_subj_df_2[["v"]])))
        #     i_subj_df_2[["j"]] = factor(x=i_subj_df_2[["j"]], 
        #                                 levels=sort(unique(i_subj_df_2[["j"]])))
        # }
        
        lst_df_subj[[i_subj]] = i_subj_df_2
    }
    
    #df_common = do.call(rbind, lst_df_common)
    #df_subj = do.call(rbind, lst_df_subj)
    
    # factorize?
    # if (factorize) {
    #     df_common[["subj"]] = factor(x=df_common[["subj"]], levels=vec_subj_uniq)
    #     df_subj[["subj"]] = factor(x=df_subj[["subj"]], levels=vec_subj_uniq)
    #     
    #     df_common[["v"]] = factor(x=df_common[["v"]], levels=vec_v_gene_uniq)
    #     df_common[["j"]] = factor(x=df_common[["v"]], levels=vec_j_gene_uniq)
    # }
    
    #return(list("common"=df_common, "subj_specific"=df_subj))
    return(list("common"=lst_df_common, "subj_specific"=lst_df_subj))
}

