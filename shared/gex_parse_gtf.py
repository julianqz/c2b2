#!/usr/bin/env python3
"""
Parse (GenCode v32) Gene Tranfer Format (GTF) file
"""

# (GenCode v32) gtf format
# https://www.gencodegenes.org/pages/data_format.html

# key-value pairs in 9th column
# format: key "value";

# examples

# gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";

# gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "lncRNA"; transcript_name "DDX11L1-202"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";

import re
import numpy as np

def parse_gtf(preparsed):
    """Parse the main annotation column in (GenCode v32) gtf.
    
    Args:
      prepased:
        A string corresponding to one row in the 9th column in the gtf.
    
    Returns:
      A tuple (gene_id, gene_name, gene_type)
    """
    # pre-checks
    assert "gene_id" in preparsed
    assert "gene_name" in preparsed
    assert "gene_type" in preparsed
    
    # split by "; " with the space behind semi-colon being optional
    # \s: whitespace
    # * causes resulting RE to match 0 or more repetitions of the preceding RE
    preparsed_split = re.split(r";\s*", preparsed)
    
    # sanity check
    # must at least contain 3 elements corresponding to gene_id, gene_name, gene_type
    # +1 to arrive at 4 b/c the trailing ; in the last element in preparsed gives rise to
    # an additional empty string
    assert len(preparsed_split) >=4
    
    # remove last entry (empty string)
    preparsed_split.pop()
    
    # for elements corresponding to gene_id, gene_type, gene_name
    # split by whitespace
    # after split, if second part surrounded by quotes, remove quotes

    # gene_id
    
    # locate element corresponding to gene_id
    gene_id_idx = np.where(["gene_id" in a for a in preparsed_split])[0]
    # expect exactly 1 match
    assert len(gene_id_idx) == 1
    # retrieve index
    gene_id_idx = gene_id_idx[0]
    # split by whitespace
    gene_id_pre = re.split(r"\s", preparsed_split[gene_id_idx])
    # expect exactly 2 elements after split
    assert len(gene_id_pre) == 2
    # get 2nd post-split element
    gene_id = gene_id_pre[1]
    # if element is wrapped by quotes (as is the case in GenCode v32 gtf)
    if gene_id[0]=='"':
        # 1:(len(gene_id)-1) covers 1:(len(gene_id)-2) inclusive
        gene_id = gene_id[1:(len(gene_id)-1)]
    assert '"' not in gene_id
    
    # gene_type
    
    gene_type_idx = np.where(["gene_type" in a for a in preparsed_split])[0]
    assert len(gene_type_idx) == 1
    gene_type_idx = gene_type_idx[0]
    
    gene_type_pre = re.split(r"\s", preparsed_split[gene_type_idx])
    assert len(gene_type_pre) == 2
    gene_type = gene_type_pre[1]
    if gene_type[0]=='"':
        gene_type = gene_type[1:(len(gene_type)-1)]
    assert '"' not in gene_type
    
    # gene_name
    
    gene_name_idx = np.where(["gene_name" in a for a in preparsed_split])[0]
    assert len(gene_name_idx) == 1
    gene_name_idx = gene_name_idx[0]
    
    gene_name_pre = re.split(r"\s", preparsed_split[gene_name_idx])
    assert len(gene_name_pre) == 2
    gene_name = gene_name_pre[1]
    if gene_name[0]=='"':
        gene_name = gene_name[1:(len(gene_name)-1)]
    assert '"' not in gene_name
    
    return (gene_id, gene_name, gene_type)
