#!/usr/bin/env bash
# Preprocess bulk NEB AbSeq data (phix removal and presto-abseq pipeline
# with UMI correction)
#
# Author: Julian Q Zhou
# Date:   2021-04-26
#
# Prereqs:  
# 1) raw fastq files from all samples should be in ${PROJ_ID}/data/raw/
# 2) raw fastq files named as follows: [sample_id]_R[12].fastq
# 3) sample_list_${PROJ_ID}.txt in ${PROJ_ID}/aux/
# 4) ${PROJ_ID}.yaml in ${PROJ_ID}/aux/


# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -J  Project ID."                        
    echo -e "  -T  Path to the top-level working dir." 
    echo -e "  -Y  Number of cores for parallelization." 
    echo -e "  -A  Whether to run phix removal (PR). Boolean."
    echo -e "  -B  Whether to run presto-abseq pipeline (PA) with UMI correction. Boolean."
    echo -e "  -C  [PR] Path to script for phix removal."
    echo -e "  -D  [PA] Path to script for presto-abseq pipeline with UMI correction."
    echo -e "  -E  [PR] Path to fastq2fasta.py."
    echo -e "  -F  [PA] Path to Python script removing inconsistent C primer and internal C alignments."
    echo -e "  -G  [PR] Directory containing phiX174 reference db."
    echo -e "  -H  [PA] Path to Read 1 FASTA primer sequences."
    echo -e "  -I  [PA] Path to Read 2 FASTA primer or template switch sequences."
    echo -e "  -K  [PA] Path to C-region FASTA sequences for the C-region internal to the primer."
    echo -e "  -L  [PA] Path to V-segment reference file."
    echo -e "  -M  [PA] The mate-pair coordinate format of the raw data."
    echo -e "  -N  [PA] Parameter that sets ${CS_KEEP}. Boolean."
    echo -e "  -O  [PA] Parameter that sets ${BOOL_PRE}. Boolean."
    echo -e "  -P  [PA] Parameter that sets ${BOOL_MID}. Boolean."
    echo -e "  -Q  [PA] Parameter that sets ${BOOL_POST}. Boolean."
    echo -e "  -R  [PA] Parameter that sets ${N_SUBSAMPLE}."
    echo -e "  -h  This message."
}

# Get commandline arguments
while getopts "J:T:Y:A:B:C:D:E:F:G:H:I:K:L:M:N:O:P:Q:R:h" OPT; do
    case "$OPT" in
    J)  PROJ_ID="${OPTARG}"
        ;;
    T)  PATH_ROOT=$(realpath "${OPTARG}")
        ;;
    Y)  NPROC="${OPTARG}"
        ;;
    A)  BOOL_PR="${OPTARG}"
        ;;
    B)  BOOL_PA="${OPTARG}"
        ;;
    C)  PATH_SCRIPT_PR=$(realpath "${OPTARG}")
        ;;
    D)  PATH_SCRIPT_PA=$(realpath "${OPTARG}")
        ;;
    E)  PATH_SCRIPT_Q2A=$(realpath "${OPTARG}")
        ;;
    F)  PATH_SCRIPT_C=$(realpath "${OPTARG}")
        ;;
    G)  PATH_REF_PHIX=$(realpath "${OPTARG}")
        ;;
    H)  PATH_PRIMER_R1=$(realpath "${OPTARG}")
        ;;
    I)  PATH_PRIMER_R2=$(realpath "${OPTARG}")
        ;;
    K)  PATH_IC=$(realpath "${OPTARG}")
        ;;
    L)  PATH_REF_V=$(realpath "${OPTARG}")
        ;;
    M)  COORD="${OPTARG}"
        ;;
    N)  BOOL_CS_KEEP="${OPTARG}"
        ;;
    O)  BOOL_PRE="${OPTARG}"
        ;;
    P)  BOOL_MID="${OPTARG}"
        ;;
    Q)  BOOL_POST="${OPTARG}"
        ;;
    R)  N_SUBSAMPLE="${OPTARG}"
        ;;
    h)  usage
        exit
        ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)  echo "Option -$OPTARG requires an argument" >&2
        exit 1
        ;;
    esac
done


# paths

PATH_PROJ="${PATH_ROOT}/${PROJ_ID}"
# no error if existing
mkdir -p "${PATH_PROJ}"

# raw fastqs
PATH_RAW="${PATH_PROJ}/data/raw/"
mkdir -p "${PATH_RAW}"

# csv, txt, logs
PATH_AUX="${PATH_PROJ}/aux/"
mkdir -p "${PATH_AUX}"

# output directory for phix removal
PATH_OUTPUT_PR="${PATH_PROJ}/data/phix/"
mkdir -p "${PATH_OUTPUT_PR}"

# output directory for presto-abseq pipeline with umi correction
PATH_OUTPUT_PA="${PATH_PROJ}/data/presto_umi/"
mkdir -p "${PATH_OUTPUT_PA}"

# overall log for looping thru sample list
PATH_LOG="${PATH_AUX}log_bcr_bulk_abseq_umi_$(date '+%m%d%Y_%H%M%S').log"

NAME_LIST="sample_list_${PROJ_ID}.txt"
PATH_LIST="${PATH_AUX}${NAME_LIST}" 

NAME_YAML="${PROJ_ID}.yaml"
PATH_YAML="${PATH_AUX}${NAME_YAML}"

MaskPrimers.py --version &> "${PATH_LOG}"
echo "NPROC=${NPROC}" &>> "${PATH_LOG}"
echo "Sample list: ${NAME_LIST}" &>> "${PATH_LOG}"
echo "Project yaml: ${NAME_YAML}" &>> "${PATH_LOG}"
echo "Primers to Read 1: ${PATH_PRIMER_R1}" &>> "${PATH_LOG}"
echo "Primers to Read 2: ${PATH_PRIMER_R2}" &>> "${PATH_LOG}"
echo "Internal-C fasta: ${PATH_IC}" &>> "${PATH_LOG}"
echo "V segment reference: ${PATH_REF_V}" &>> "${PATH_LOG}"
echo "BOOL_PR: ${BOOL_PRE}" &>> "${PATH_LOG}"
echo "BOOL_PRE: ${BOOL_PRE}" &>> "${PATH_LOG}"
echo "BOOL_MID: ${BOOL_MID}" &>> "${PATH_LOG}"
echo "BOOL_POST: ${BOOL_POST}" &>> "${PATH_LOG}"
echo "N_SUBSAMPLE: ${N_SUBSAMPLE}" &>> "${PATH_LOG}"
echo "CS_KEEP: ${BOOL_CS_KEEP}" &>> "${PATH_LOG}"


# phix removal

if $BOOL_PR; then

    echo "***** phix removal *****" &>> "${PATH_LOG}"

    N_LINES=$(wc -l < "${PATH_LIST}")
    echo "N_LINES: ${N_LINES}" &>> "${PATH_LOG}"

    for ((IDX=1; IDX<=${N_LINES}; IDX++)); do

        # read sample ID from file
        CUR_ID=$(sed "${IDX}q;d" "${PATH_LIST}") 

        echo "IDX: ${IDX}; CUR_ID: ${CUR_ID}" &>> "${PATH_LOG}"

        
        #* sample input fastq files
        RAW_1="${PATH_RAW}/${CUR_ID}_R1.fastq"
        RAW_2="${PATH_RAW}/${CUR_ID}_R2.fastq"

        # sample-specific log for phix removal
        PATH_LOG_PR_ID="${PATH_AUX}log_PR_${IDX}_${CUR_ID}_$(date '+%m%d%Y_%H%M%S').log"

        # R1
        echo "  - R1" &>> "${PATH_LOG}"
        echo "---------------- R1 ----------------" &> "${PATH_LOG_PR_ID}"
        "${PATH_SCRIPT_PR}" \
            -s "${RAW_1}" \
            -r "${PATH_REF_PHIX}" \
            -n "${CUR_ID}_R1" \
            -o "${PATH_OUTPUT_PR}" \
            -p "${NPROC}" \
            -t "${PATH_SCRIPT_Q2A}" \
            &>> "${PATH_LOG_PR_ID}"

        # R2
        echo "  - R2" &>> "${PATH_LOG}"
        echo "---------------- R2 ----------------" &>> "${PATH_LOG_PR_ID}"  
        "${PATH_SCRIPT_PR}" \
            -s "${RAW_2}" \
            -r "${PATH_REF_PHIX}" \
            -n "${CUR_ID}_R2" \
            -o "${PATH_OUTPUT_PR}" \
            -p "${NPROC}" \
            -t "${PATH_SCRIPT_Q2A}" \
            &>> "${PATH_LOG_PR_ID}"

        # input files to presto-abseq pipeline are now output files from phix removal
        # phix removal adds "_R2_nophix_selected.fastq" to ${CUR_ID}
        PATH_INPUT="${PATH_OUTPUT_PR}"
        SUFFIX_1="_R1_nophix_selected.fastq"
        SUFFIX_2="_R2_nophix_selected.fastq"
        
    done
    
else
    
    PATH_INPUT="${PATH_RAW}"
    SUFFIX_1="_R1.fastq"
    SUFFIX_2="_R2.fastq"

fi


# presto-abseq pipeline with umi correction

if $BOOL_PA; then

    echo "***** presto-abseq with umi correction *****" &>> "${PATH_LOG}"

    "${PATH_SCRIPT_PA}" \
        -a "${PATH_LIST}" \
        -b "${PATH_INPUT}" \
        -c "${SUFFIX_1}" \
        -d "${SUFFIX_2}" \
        -e "${N_SUBSAMPLE}" \
        -f "${PATH_PRIMER_R1}" \
        -g "${PATH_PRIMER_R2}" \
        -i "${PATH_IC}" \
        -j "${PATH_REF_V}" \
        -k "${PATH_YAML}" \
        -l "${PATH_OUTPUT_PA}" \
        -m "${COORD}" \
        -n "${NPROC}" \
        -o "${PATH_SCRIPT_C}" \
        -p "${PATH_SCRIPT_Q2A}" \
        -q "${BOOL_CS_KEEP}" \
        -r "${BOOL_PRE}" \
        -s "${BOOL_MID}" \
        -t "${BOOL_POST}" \
        &> "${PATH_LOG}"

fi

echo "Finished" &>> "${PATH_LOG}"
