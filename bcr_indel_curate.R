
#### example sequence object ####

# sequence
# event: insertion; deletion; substitution
# event coordinates: start, end
# event length 

set.seed(394752)
seqinr::c2s(sample(x=c("A","T","G","C"), size=100, replace=T))

### del at left edge + ins at right edge

#                                                                                                                       1111
# ---                               23    33          4 4               66      77                    9    9            0111
# 321123456   7                     90    56          7 9               56      34                    5    6            9012 
#    TACTTA   GCCGAACATAGCAGATAGCAAGTAGTAGGGCAATTTATCTGTGGATACCGCAACACCAAGGTTAACCAGAATGTGGCAAGGCCGTGGTT    GTTTTAGTCTTTGTTTT  obsv
# GGGTACTTAATGGCCGAACATAGCAGATAGCAAGT      GCAATTTATCTCCCGATACCGCAACACCAA        AGAATGTGGCAAGGCCGTGGTTATGCGTTTTAGTCTTTG      germ

# ___TACTTA___GCCGAACATAGCAGATAGCAAGT+AGTAGG+GCAATTTATCT>GTG<GATACCGCAACACCAA+GGTTAACC+AGAATGTGGCAAGGCCGTGGTT____GTTTTAGTCTTTG+TTTT+

# idea: use code to generate ^^^; manually turn ___/++/>< in Word to highlights; check events for pos; check printed table() for clonal members

seq_obj_1 = vector(mode="list", length=2)
names(seq_obj_1) = c("sequence", "events")
seq_obj_1[["sequence"]] = "TACTTAGCCGAACATAGCAGATAGCAAGTAGTAGGGCAATTTATCTGTGGATACCGCAACACCAAGGTTAACCAGAATGTGGCAAGGCCGTGGTTGTTTTAGTCTTTGTTTT"

events_1 = data.frame(matrix(NA, nrow=7, ncol=5))
colnames(events_1) = c("event_type", "event_start", "event_end",
                     "event_len", "event_seq")
events_1[["event_type"]] = c("deletion", "deletion", "insertion", "substitution",
                           "insertion", "deletion", "insertion")
events_1[["event_start"]] = c(-3,6,30,47,66,95,109)
events_1[["event_end"]] = c(-1,7,35,49,73,96,112)
events_1[["event_len"]] = c(3,3,6,3,8,4,4)
events_1[["event_seq"]] = c("GGG", "ATG", "AGTAGG", "GTG", "GGTTAACC", "ATGC", "TTTT")
stopifnot(!any(is.na(events_1)))

seq_obj_1[["events"]] = events_1
rm(events_1)

### ins at left edge and del at right edge

#                                                                                                                     111111
#             1                     33    33         45555              66      77                    9    9          111111
# 123456789   0                     23    89         90123              89      67                    8    9          012345
# GGGTACTTA   GCCGAACATAGCAGATAGCAAGTAGTAGGGCAATTTATCTGTGGATACCGCAACACCAAGGTTAACCAGAATGTGGCAAGGCCGTGGTT    GTTTTAGTCTTTG      obsv
#    TACTTAATGGCCGAACATAGCAGATAGCAAGT      GCAATTTATCTCCCGATACCGCAACACCAA        AGAATGTGGCAAGGCCGTGGTTATGCGTTTTAGTCTTTGTTTT  germ

# +GGG+TACTTA___GCCGAACATAGCAGATAGCAAGT+AGTAGG+GCAATTTATCT>GTG<GATACCGCAACACCAA+GGTTAACC+AGAATGTGGCAAGGCCGTGGTT____GTTTTAGTCTTTG____

# idea: use code to generate ^^^; manually turn ___/++/>< in Word to highlights; check events for pos; check printed table() for clonal members

seq_obj_2 = vector(mode="list", length=2)
names(seq_obj_2) = c("sequence", "events")
seq_obj_2[["sequence"]] = "GGGTACTTAGCCGAACATAGCAGATAGCAAGTAGTAGGGCAATTTATCTGTGGATACCGCAACACCAAGGTTAACCAGAATGTGGCAAGGCCGTGGTTGTTTTAGTCTTTG"

events_2 = data.frame(matrix(NA, nrow=7, ncol=5))
colnames(events_2) = c("event_type", "event_start", "event_end",
                     "event_len", "event_seq")
events_2[["event_type"]] = c("insertion", "deletion", "insertion", "substitution",
                           "insertion", "deletion", "deletion")
events_2[["event_start"]] = c(1,9,33,50,69,98,112)
events_2[["event_end"]] = c(3,10,38,52,76,99,115)
events_2[["event_len"]] = c(3,3,6,3,8,4,4)
events_2[["event_seq"]] = c("GGG", "ATG", "AGTAGG", "GTG", "GGTTAACC", "ATGC", "TTTT")
stopifnot(!any(is.na(events_2)))

seq_obj_2[["events"]] = events_2
rm(events_2)


#### functions ####

# constants
NAME_SEQ = "sequence"
NAME_EVENTS = "events"
NAME_EVENT_TYPE = "event_type"
NAME_EVENT_START = "event_start"
NAME_EVENT_END = "event_end"
NAME_EVENT_LEN = "event_len"
NAME_EVENT_SEQ = "event_seq"
EVENTS_COLNAMES = c(NAME_EVENT_TYPE, NAME_EVENT_START, NAME_EVENT_END, 
                    NAME_EVENT_LEN, NAME_EVENT_SEQ)
EVENT_TYPES = c("insertion", "deletion", "substitution")


# Check the starting and end positions of the local alignment involved in an event,
# given the event type and the length of the full sequence
# Called by `check_seq_obj`
# Input:
# - seq_end: same as length of seq
# Output: 
# - TRUE/FALSE
check_event_start_end = function(event_type, event_start, event_end, seq_end) {
    
    stopifnot( event_end >= event_start )
    
    if (event_type=="substitution") {

        if (event_start >= 1 & event_end <= seq_end) {
            return(TRUE)
        } else {
            return(FALSE)
        }
        
    } else if (event_type=="deletion") {
        # at left edge
        bool_left_edge = (event_end == -1)
        # at right edge
        bool_right_edge = (event_start == seq_end+1)
        # not at edge
        bool_not_edge = (event_start > 1 & event_end < seq_end)
        
        if (sum(c(bool_left_edge, bool_right_edge, bool_not_edge))==1) {
            return(TRUE)
        } else {
            return(FALSE)
        }
        
    } else if (event_type=="insertion") {
        # at left edge
        bool_left_edge = (event_start == 1)
        # at right edge
        bool_right_edge = (event_end == seq_end)
        # not at edge
        bool_not_edge = (event_start > 1 & event_end < seq_end)
        
        if (bool_left_edge) {
            stopifnot(event_end<seq_end)
        }
        if (bool_right_edge) {
            stopifnot(event_start>1)
        }
        
        if (sum(c(bool_left_edge, bool_right_edge, bool_not_edge))==1) {
            return(TRUE)
        } else {
            return(FALSE)
        }
        
    } else {
        cat("Unexpected event_type.\n")
        return(FALSE)
    }
}

# Calculate the length of the local alignment involved in an event,
# given the event's start and end position, as well as 
#       the length of the full sequence
# Called by `check_seq_obj`
# Input: 
# - seq_end: same as length of seq
# Output:
# - If check passes, NA if deletion in the middle of full sequence; otherwise an integer
# - If check fails, an error will be raised
calculate_event_len_from_start_end = function(event_type, event_start, event_end, seq_end) {
    
    if (event_start > 1 & event_end < seq_end) {
        # not at one of two edges
        
        if (event_type %in% c("substitution", "insertion")) {
            len = event_end - event_start + 1
            return(len)
            
        } else if (event_type=="deletion") {
            # can't calculate from start/end
            return(NA)
            
        } else {
            stop("Unexpected event_type.")
        }
        
    } else {
        # at one of two edges
        len = event_end - event_start + 1
        return(len)
    }

}


# Check the validity of a sequence object ($sequence + $events)
# Input:
# - seq_obj: list with pre-defined specifications containing $sequence and $events
# - legal_alphabet: character vector defining the acceptable characters in $sequence and $events$event_seq
# Output:
# - If check passes, TRUE
check_seq_obj = function(seq_obj, 
                         legal_alphabet=c("A","T","G","C","N",".","-")) {
    require(seqinr)
    
    ### overall data structure
    stopifnot(is.list(seq_obj))
    stopifnot(length(seq_obj)==2)
    stopifnot(all.equal(names(seq_obj), c(NAME_SEQ, NAME_EVENTS)))
    
    ### sequence
    seq = seq_obj[[NAME_SEQ]]
    stopifnot(is.character(seq))
    
    # check all characters are in legal alphabet
    seq = toupper(seq)
    seq_c = s2c(seq)
    stopifnot(all(seq_c %in% legal_alphabet))
    
    # seq length (for use later with checking $events)
    seq_end = nchar(seq)
    
    ### events
    
    events = seq_obj[[NAME_EVENTS]]
    stopifnot(is.data.frame(events))
    stopifnot(all.equal(colnames(events), EVENTS_COLNAMES))
    
    # ok if $events is empty (0-row data.frame)
    if (nrow(events)>0) {
        ## all event types are ones as expected
        stopifnot(all(events[[NAME_EVENT_TYPE]] %in% EVENT_TYPES))
        
        ## range (start/end) should be acceptable
        bool_event_range_check = sapply(1:nrow(events), 
                                        function(i){ 
                                            check_event_start_end(event_type=events[[NAME_EVENT_TYPE]][i],
                                                                  event_start=events[[NAME_EVENT_START]][i],
                                                                  event_end=events[[NAME_EVENT_END]][i],
                                                                  seq_end=seq_end) 
                                        })
        stopifnot(all(bool_event_range_check))
        
        ## len should be non-negative
        stopifnot(all(events[[NAME_EVENT_LEN]]>0))
        
        ## range and len should match logically
        # calculate len from range (start/end)
        event_len_calc_init = sapply(1:nrow(events), 
                                     function(i){
                                         calculate_event_len_from_start_end(event_type=events[[NAME_EVENT_TYPE]][i],
                                                                            event_start=events[[NAME_EVENT_START]][i],
                                                                            event_end=events[[NAME_EVENT_END]][i],
                                                                            seq_end=seq_end)})
        # expect either >0 integer or NA
        stopifnot(all(event_len_calc_init>0, na.rm=T))
        # if non-NA, value should match events$event_len
        stopifnot(all(events[[NAME_EVENT_LEN]]==event_len_calc_init, na.rm=T))
        
        ## event_len and event_seq should agree
        stopifnot(all.equal(events[[NAME_EVENT_LEN]],
                            nchar(events[[NAME_EVENT_SEQ]])))
        
        ## event_seq should contain only legal alphabet
        stopifnot(all(sapply(events[[NAME_EVENT_SEQ]],
                             function(x){ all(s2c(x) %in% legal_alphabet) })))
    }
    
    return(TRUE)
}

# expect TRUE
check_seq_obj(seq_obj_1)
check_seq_obj(seq_obj_2)


# Assumes that `seq_obj` has passed `check_seq_obj`
# deletion: _ for each deleted position
# insertion: + for start and end of insertion
# substitution: > for start and < for end of substitution

parse_single_event = function(event_idx_cur, seq_obj) {
    
    events = seq_obj[[NAME_EVENTS]]
    sequence = seq_obj[[NAME_SEQ]]
    
    n_events = nrow(events)
    seq_end = nchar(sequence)

    event_type_cur = events[[NAME_EVENT_TYPE]][event_idx_cur]
    event_start_cur = events[[NAME_EVENT_START]][event_idx_cur]
    event_end_cur = events[[NAME_EVENT_END]][event_idx_cur]
    event_len_cur = events[[NAME_EVENT_LEN]][event_idx_cur]
    
    ### part 1
    
    if (event_idx_cur==1) {
        
        # if 1st event, check whether there's seq before current event point
        
        if ( (event_type_cur %in% c("insertion", "substitution") & event_start_cur==1) | 
             (event_type_cur=="deletion" & event_start_cur<0) ) {
            # no seq before current event point
            # i.e. current ins/sub is at left edge
            seq_1_pos = NA
            seq_1_str = ""
            
        } else {
            # there's seq before current event point
            seq_1_start = 1
            
            if (event_type_cur %in% c("insertion", "substitution")) {
                seq_1_end = event_start_cur-1
            } else if (event_type_cur=="deletion") {
                seq_1_end = event_start_cur
            }
            
            seq_1_pos = seq_1_start:seq_1_end
            seq_1_str = substr(sequence, seq_1_start, seq_1_end)
        }
        
    } else {
        
        event_idx_bf = event_idx_cur-1
        event_end_bf = events[[NAME_EVENT_END]][event_idx_bf]
        event_type_bf = events[[NAME_EVENT_TYPE]][event_idx_bf]

        if (event_type_bf %in% c("insertion", "substitution")) {
            seq_1_start = event_end_bf+1
        } else if (event_type_bf=="deletion") {
            if (event_end_bf!=-1) {
                seq_1_start = event_end_bf
            } else {
                # when event_idx_cur=2, event_idx_bf=1, 
                # and 1st event is a deletion at left edge
                seq_1_start = 1
            }
        }
        
        if (event_type_cur %in% c("insertion", "substitution")) {
            seq_1_end = event_start_cur-1
        } else if (event_type_cur=="deletion") {
            
            if (event_start_cur > seq_end) {
                # del at right edge
                seq_1_end = event_start_cur-1
            } else {
                # del not at right edge
                seq_1_end = event_start_cur
            }
            
        }
        
        seq_1_pos = seq_1_start:seq_1_end
        seq_1_str = substr(sequence, seq_1_start, seq_1_end)
        
    }
    
    ### part 2
    
    if (event_type_cur %in% c("insertion", "substitution")) {
        seq_2_start = event_start_cur
        seq_2_end = event_end_cur
        seq_2_pos = seq_2_start:seq_2_end
        seq_2_str_stem = substr(sequence, seq_2_start, seq_2_end)
        
        if (event_type_cur=="insertion") {
            seq_2_str = paste0("+", seq_2_str_stem, "+")
        } else if (event_type_cur=="substitution") {
            seq_2_str = paste0(">", seq_2_str_stem, "<")
        }
        
    } else if (event_type_cur=="deletion") {
        seq_2_pos = NA
        seq_2_str = paste0(rep("_", event_len_cur), collapse="")
    }
    
    ### part 3
    
    if (event_idx_cur==n_events) {
        # if last event, check whether there's seq after current event point
        if (event_end_cur>=seq_end) {
            # no seq after current event point
            seq_3_pos = NA
            seq_3_str = ""
        } else {
            # there's seq after current event point
            if (event_type_cur %in% c("insertion", "substitution")) {
                seq_3_start = event_end_cur+1
                
            } else if (event_type_cur=="deletion") {
                seq_3_start = event_end_cur
            }
            seq_3_end = seq_end
            seq_3_pos = seq_3_start:seq_3_end
            seq_3_str = substr(sequence, seq_3_start, seq_3_end)
        }
    } else {
        seq_3_pos = NA
        seq_3_str = ""
    }
    
    
    cat(event_idx_cur, "\n")
    cat("seq_1_pos:\n")
    print(seq_1_pos)
    cat("seq_1_str:", seq_1_str, "\n")
    cat("seq_2_pos:\n")
    print(seq_2_pos)
    cat("seq_2_str:", seq_2_str, "\n")
    cat("seq_3_pos:\n")
    print(seq_3_pos)
    cat("seq_3_str:", seq_3_str, "\n")
    
    lst_return = vector(mode="list", length=6)
    names(lst_return) = paste0("seq_", rep(1:3, each=2), rep(c("_pos", "_str"), times=3))
    lst_return[["seq_1_pos"]] = seq_1_pos
    lst_return[["seq_1_str"]] = seq_1_str
    lst_return[["seq_2_pos"]] = seq_2_pos
    lst_return[["seq_2_str"]] = seq_2_str
    lst_return[["seq_3_pos"]] = seq_3_pos
    lst_return[["seq_3_str"]] = seq_3_str
    
    return(lst_return)
}


parse_seq_obj = function(seq_obj) {
    
    n_events = nrow(seq_obj[[NAME_EVENTS]])
    sequence = seq_obj[[NAME_SEQ]]
    seq_end = nchar(sequence)
    
    lst_parsed = vector(mode="list", length=n_events)
    
    for (i_event in 1:n_events) {
        lst_parsed[[i_event]] = parse_single_event(event_idx_cur=i_event, 
                                                   seq_obj=seq_obj)
    }
    
    # All the seq_*_pos combined should have no duplicate and 
    # should cover all of the positions in input sequence
    
    .collapse_seq_pos = function(x) {
        return(c(x[["seq_1_pos"]], x[["seq_2_pos"]], x[["seq_3_pos"]]))
    }
    
    vec_seq_pos_w_na = unlist(lapply(lst_parsed, .collapse_seq_pos))
    vec_seq_pos_no_na = vec_seq_pos_w_na[!is.na(vec_seq_pos_w_na)]
    
    stopifnot(length(vec_seq_pos_no_na)==seq_end)
    stopifnot(!any(duplicated(vec_seq_pos_no_na)))
    stopifnot(all.equal(sort(vec_seq_pos_no_na), 1:seq_end))
    
    .collapse_seq_str = function(x) {
        return(c(x[["seq_1_str"]], x[["seq_2_str"]], x[["seq_3_str"]]))
    }
    
    # All the seq_*_str combined, and with +_<> removed, should match input sequence
    vec_seq_str_w_annotation = unlist(lapply(lst_parsed, .collapse_seq_str))
    str_annotated = paste(vec_seq_str_w_annotation, collapse="")
    # remove all occurrences of "+", "_", ">", "<"
    str_restored = gsub(pattern="[_><+]", replacement="", x=str_annotated, fixed=F)
    
    stopifnot(str_restored==sequence)
    
    lst_return = vector(mode="list", length=2)
    names(lst_return) = c("lst_parsed", "str_annotated")
    lst_return[["lst_parsed"]] = lst_parsed
    lst_return[["str_annotated"]] = str_annotated
    
    return(lst_return)
}


parse_seq_obj(seq_obj_2)
parse_seq_obj(seq_obj_2)[["str_annotated"]]

parse_seq_obj(seq_obj_1)
parse_seq_obj(seq_obj_1)[["str_annotated"]]


# TODO

# rename parse_
# add doc
# lots of additional test cases and edge cases
# gather tough cases (see ELN 2025-05-09)


#### temp ####

obsv = "TACTTAGCCGAACATAGCAGATAGCAAGTAGTAGGGCAATTTATCTGTGGATACCGCAACACCAAGGTTAACCAGAATGTGGCAAGGCCGTGGTTGTTTTAGTCTTTG"
germ = "TACTTAATGGCCGAACATAGCAGATAGCAAGTGCAATTTATCTCCCGATACCGCAACACCAAAGAATGTGGCAAGGCCGTGGTTATGCGTTTTAGTCTTTG"

#install.packages("text.alignment")

text.alignment::smith_waterman(b=germ, a=obsv, type="character", lower=F, edit_mark="-",
                               match=2, mismatch=-1, gap=-1)

#BiocManager::install("DECIPHER")

library(Biostrings)
library(DECIPHER)

seq_set = readDNAStringSet("~/Desktop/test.fasta", format="fasta")
seq_aligned = AlignSeqs(seq_set, gapOpening=-7)
seq_aligned

BrowseSeqs(seq_aligned)

seq_aligned_2 = AdjustAlignment(seq_aligned)
BrowseSeqs(seq_aligned_2)
