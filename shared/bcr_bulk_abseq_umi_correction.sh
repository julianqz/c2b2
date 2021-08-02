#!/usr/bin/env bash
# Preprocess bulk NEB AbSeq data (phix removal and presto-abseq pipeline
# with UMI correction)
#
# Author: Julian Q Zhou
# Date:   2021-04-26
#
# Prereqs:  
# 1) sample_list_${PROJ_ID}.csv in ${PROJ_ID}/aux/
#    - 4 comma-separated columns
#    - [sample id],[directory],[R1 filename],[R2 filename]
#      Eg: b1,/home/raw,R1.fastq.gz,R2.fastq.gz
# 2) ${PROJ_ID}.yaml in ${PROJ_ID}/aux/


# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -A  Project ID."                        
    echo -e "  -B  Path to the top-level working dir." 
    echo -e "  -C  Number of cores for parallelization." 
    echo -e "  -D  Whether to run phix removal (PR). Boolean."
    echo -e "  -E  Whether PR had been run already. Boolean."
    echo -e "  -F  Whether to run presto-abseq pipeline (PA) with UMI correction. Boolean."
    echo -e "  -G  [PR] Path to script for phix removal."
    echo -e "  -H  [PA] Path to script for presto-abseq pipeline with UMI correction."
    echo -e "  -I  [PR] Path to fastq2fasta.py."
    echo -e "  -J  [PA] Path to Python script removing inconsistent C primer and internal C alignments."
    echo -e "  -K  [PR] Directory containing phiX174 reference db."
    echo -e "  -L  [PA] Path to Read 1 FASTA primer sequences."
    echo -e "  -M  [PA] Path to Read 2 FASTA primer or template switch sequences."
    echo -e "  -N  [PA] Path to C-region FASTA sequences for the C-region internal to the primer."
    echo -e "  -O  [PA] Path to V-segment reference file."
    echo -e "  -P  [PA] The mate-pair coordinate format of the raw data."
    echo -e "  -Q  [PA] Parameter that sets ${CS_KEEP}. Boolean."
    echo -e "  -R  [PA] Parameter that sets ${BOOL_PRE}. Boolean."
    echo -e "  -S  [PA] Parameter that sets ${BOOL_MID}. Boolean."
    echo -e "  -T  [PA] Parameter that sets ${BOOL_POST}. Boolean."
    echo -e "  -U  [PA] Parameter that sets ${N_SUBSAMPLE}."
    echo -e "  -h  This message."
}

# Get commandline arguments
while getopts "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:h" OPT; do
    case "$OPT" in
    A)  PROJ_ID="${OPTARG}"
        ;;
    B)  PATH_ROOT=$(realpath "${OPTARG}")
        ;;
    C)  NPROC="${OPTARG}"
        ;;
    D)  BOOL_PR="${OPTARG}"
        ;;
    E)  RAN_PR_ALREADY="${OPTARG}"
        ;;
    F)  BOOL_PA="${OPTARG}"
        ;;
    G)  PATH_SCRIPT_PR=$(realpath "${OPTARG}")
        ;;
    H)  PATH_SCRIPT_PA=$(realpath "${OPTARG}")
        ;;
    I)  PATH_SCRIPT_Q2A=$(realpath "${OPTARG}")
        ;;
    J)  PATH_SCRIPT_C=$(realpath "${OPTARG}")
        ;;
    K)  PATH_REF_PHIX=$(realpath "${OPTARG}")
        ;;
    L)  PATH_PRIMER_R1=$(realpath "${OPTARG}")
        ;;
    M)  PATH_PRIMER_R2=$(realpath "${OPTARG}")
        ;;
    N)  PATH_IC=$(realpath "${OPTARG}")
        ;;
    O)  PATH_REF_V=$(realpath "${OPTARG}")
        ;;
    P)  COORD="${OPTARG}"
        ;;
    Q)  BOOL_CS_KEEP="${OPTARG}"
        ;;
    R)  BOOL_PRE="${OPTARG}"
        ;;
    S)  BOOL_MID="${OPTARG}"
        ;;
    T)  BOOL_POST="${OPTARG}"
        ;;
    U)  N_SUBSAMPLE="${OPTARG}"
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

NAME_LIST="sample_list_${PROJ_ID}.csv"
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

        # readline from csv file
        CUR_LINE=$(sed "${IDX}q;d" "${PATH_LIST}")

        # split strings in unix
        # https://linuxhint.com/bash_split_examples/ 
        # following example 2

        # $IFS: internal field separator (default is white space)
        # -r: read backslash (\) as a character rather than escape character
        # -a: store split words into an array variable
        IFS=","
        read -a strarr <<< "${CUR_LINE}"

        # sample ID
        CUR_ID=${strarr[0]}

        # path containing input fasta files
        PATH_RAW=${strarr[1]}

        # sample input fastq files
        FN_1=${strarr[2]}
        FN_2=${strarr[3]}
        RAW_1="${PATH_RAW}/${FN_1}"
        RAW_2="${PATH_RAW}/${FN_2}"


        echo "IDX: ${IDX}; CUR_ID: ${CUR_ID}" &>> "${PATH_LOG}"
        echo " - R1: ${RAW_1}" &>> "${PATH_LOG}"
        echo " - R2: ${RAW_2}" &>> "${PATH_LOG}"

        
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
    
    if $RAN_PR_ALREADY; then
        # ran PR before and not repeating the process
        PATH_INPUT="${PATH_OUTPUT_PR}"
        SUFFIX_1="_R1_nophix_selected.fastq"
        SUFFIX_2="_R2_nophix_selected.fastq"
    else
        # skipping PR
        # the following settings trigger a special condition in the presto-abseq_UMI_correction pipe
        PATH_INPUT="SKIP_PR"
        SUFFIX_1="SKIP_PR"
        SUFFIX_2="SKIP_PR"
    fi
    
fi


# presto-abseq pipeline with umi correction

if $BOOL_PA; then

    echo "***** presto-abseq with umi correction *****" &>> "${PATH_LOG}"
    echo "PATH_INPUT: ${PATH_INPUT}" &>> "${PATH_LOG}"
    echo "SUFFIX_1: ${SUFFIX_1}" &>> "${PATH_LOG}"
    echo "SUFFIX_2: ${SUFFIX_2}" &>> "${PATH_LOG}"


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
