
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
# 321123456   7                     90    56          6 8               56      34                    5    6            9012 
#    TACTTA   GCCGAACATAGCAGATAGCAAGTAGTAGGGCAATTTATCTGTGGATACCGCAACACCAAGGTTAACCAGAATGTGGCAAGGCCGTGGTT    GTTTTAGTCTTTGTTTT  obsv
# GGGTACTTAATGGCCGAACATAGCAGATAGCAAGT      GCAATTTATCTCCCGATACCGCAACACCAA        AGAATGTGGCAAGGCCGTGGTTATGCGTTTTAGTCTTTG      germ

# ___TACTTA___GCCGAACATAGCAGATAGCAAGT+AGTAGG+GCAATTTATCT>GTG<GATACCGCAACACCAA+GGTTAACC+AGAATGTGGCAAGGCCGTGGTT____GTTTTAGTCTTTG

# idea: use code to generate ^^^; manually turn ___/++/>< in Word to highlights; check events for pos; check printed table() for clonal members

seq_obj_1 = vector(mode="list", length=2)
names(seq_obj_1) = c("sequence", "events")
seq_obj_1[["sequence"]] = "TACTTAGCCGAACATAGCAGATAGCAAGTAGTAGGGCAATTTATCTGTGGATACCGCAACACCAAGGTTAACCAGAATGTGGCAAGGCCGTGGTTGTTTTAGTCTTTGTTTT"

events_1 = data.frame(matrix(NA, nrow=7, ncol=5))
colnames(events_1) = c("event_type", "event_start", "event_end",
                     "event_len", "event_seq")
events_1[["event_type"]] = c("deletion", "deletion", "insertion", "substitution",
                           "insertion", "deletion", "insertion")
events_1[["event_start"]] = c(-3,6,30,46,66,95,109)
events_1[["event_end"]] = c(-1,7,35,48,73,96,112)
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

# ___TACTTA___GCCGAACATAGCAGATAGCAAGT+AGTAGG+GCAATTTATCT>GTG<GATACCGCAACACCAA+GGTTAACC+AGAATGTGGCAAGGCCGTGGTT____GTTTTAGTCTTTG

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
generate_indel_annotation = function(seq_obj) {
    # deletion: _ for each deleted position
    # insertion: + for start and end of insertion
    # substitution: > for start and < for end of substitution
    
    seq_orig = seq_obj[[NAME_SEQ]]
    events = seq_obj[[NAME_EVENTS]]
}



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
