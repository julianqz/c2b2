#!/usr/bin/env bash

# Author: Julian Q. Zhou
# https://github.com/julianqz
# Date:   2021-04-20
#
# Preprocess bulk NEB AbSeq data (phix removal and presto-abseq pipeline)
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
    echo -e "  -J  Project ID."                        
    echo -e "  -T  Path to the top-level working dir." 
    echo -e "  -Y  Number of cores for parallelization." 
    echo -e "  -A  Whether to run phix removal (PR). Boolean."
    echo -e "  -P  If not running phix removal (PR), Whether phix-removed files already exist. Boolean."
    echo -e "  -B  Whether to run presto-abseq pipeline (PA). Boolean."
    echo -e "  -C  [PR] Path to script for phix removal."
    echo -e "  -D  [PA] Path to script for presto-abseq pipeline."
    echo -e "  -E  [PR] Path to fastq2fasta.py."
    echo -e "  -F  [PA] Path to Python script removing inconsistent C primer and internal C alignments."
    echo -e "  -G  [PR] Directory containing phiX174 reference db."
    echo -e "  -H  [PA] Path to Read 1 FASTA primer sequences."
    echo -e "  -I  [PA] Path to Read 2 FASTA primer or template switch sequences."
    echo -e "  -K  [PA] Path to C-region FASTA sequences for the C-region internal to the primer."
    echo -e "  -L  [PA] Path to V-segment reference file."
    echo -e "  -M  [PA] The mate-pair coordinate format of the raw data."
    echo -e "  -N  [PA] Parameter that sets ${CS_KEEP}. Boolean."
    echo -e "  -h  This message."
}

# Get commandline arguments
while getopts "J:T:Y:A:P:B:C:D:E:F:G:H:I:K:L:M:N:h" OPT; do
    case "$OPT" in
    J)  PROJ_ID="${OPTARG}"
        ;;
    T)  PATH_ROOT=$(realpath "${OPTARG}")
        ;;
    Y)  NPROC="${OPTARG}"
        ;;
    A)  BOOL_PR_RUN="${OPTARG}"
        ;;
    P)  BOOL_PR_EXIST="${OPTARG}"
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

# output directory for presto-abseq pipeline
PATH_OUTPUT_PA="${PATH_PROJ}/data/presto/"
mkdir -p "${PATH_OUTPUT_PA}"

# overall log for looping thru sample list
PATH_LOG="${PATH_AUX}log_bcr_bulk_abseq_$(date '+%m%d%Y_%H%M%S').log"

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
echo "CS_KEEP: ${BOOL_CS_KEEP}" &>> "${PATH_LOG}"

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


    # phix removal
    if $BOOL_PR_RUN; then

        echo "- phix removal" &>> "${PATH_LOG}"

        # sample-specific log for phix removal
        PATH_LOG_PR_ID="${PATH_AUX}log_PR_${IDX}_${CUR_ID}_$(date '+%m%d%Y_%H%M%S').log"

        # R1
        echo "   - R1" &>> "${PATH_LOG}"
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
        echo "   - R2" &>> "${PATH_LOG}"
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
        INPUT_ABSEQ_1="${PATH_OUTPUT_PR}/${CUR_ID}_R1_nophix_selected.fastq"
        INPUT_ABSEQ_2="${PATH_OUTPUT_PR}/${CUR_ID}_R2_nophix_selected.fastq"

    else

        if $BOOL_PR_EXIST; then
            # input files to presto-abseq pipeline are existing phix-removed files
             INPUT_ABSEQ_1="${PATH_OUTPUT_PR}/${CUR_ID}_R1_nophix_selected.fastq"
             INPUT_ABSEQ_2="${PATH_OUTPUT_PR}/${CUR_ID}_R2_nophix_selected.fastq"
        else
            # input files to presto-abseq pipeline are original files
            INPUT_ABSEQ_1="${RAW_1}"
            INPUT_ABSEQ_2="${RAW_2}"
        fi

    fi


    # presto-abseq pipeline
    if $BOOL_PA; then

        echo "- presto-abseq pipeline" &>> "${PATH_LOG}"

        # sample-specific log for presto-abseq pipeline
        PATH_LOG_PA_ID="${PATH_AUX}log_PA_${IDX}_${CUR_ID}_$(date '+%m%d%Y_%H%M%S').log"

        # sample-specific presto directory
        PATH_OUTPUT_PA_ID="${PATH_OUTPUT_PA}${CUR_ID}/"

        # if existing, remove first
        if [ -d "${PATH_OUTPUT_PA_ID}" ]; then
            echo "    Removed pre-exisitng folder for ${CUR_ID}" &>> "${PATH_LOG}"
            rm -r "${PATH_OUTPUT_PA_ID}"
        fi

        echo "    Created new presto-abseq folder for ${CUR_ID}" &>> "${PATH_LOG}"
        mkdir "${PATH_OUTPUT_PA_ID}"

        "${PATH_SCRIPT_PA}" \
            -1 "${INPUT_ABSEQ_1}" \
            -2 "${INPUT_ABSEQ_2}" \
            -j "${PATH_PRIMER_R1}" \
            -v "${PATH_PRIMER_R2}" \
            -c "${PATH_IC}" \
            -r "${PATH_REF_V}" \
            -y "${PATH_YAML}" \
            -n "${CUR_ID}" \
            -o "${PATH_OUTPUT_PA_ID}" \
            -x "${COORD}" \
            -p "${NPROC}" \
            -t "${PATH_SCRIPT_C}" \
            -s "${BOOL_CS_KEEP}" \
            &> "${PATH_LOG_PA_ID}"

        # rename prestor report
        # name generated from yaml is the same for all samples b/c one yaml file used for entire project
        mv "${PATH_OUTPUT_PA_ID}report/"*.html "${PATH_OUTPUT_PA_ID}report/${CUR_ID}_prestor.html"

        # convert fastq to fasta
        # creates ${CUR_ID}-final_collapse-unique_atleast-2.fasta
        cd "${PATH_OUTPUT_PA_ID}"

        echo "    Converting output fastq to fasta" &>> "${PATH_LOG}"

        "${PATH_SCRIPT_Q2A}" \
            "${PATH_OUTPUT_PA_ID}${CUR_ID}-final_collapse-unique_atleast-2.fastq" \
            &>> "${PATH_LOG_PA_ID}"

        "${PATH_SCRIPT_Q2A}" \
            "${PATH_OUTPUT_PA_ID}${CUR_ID}-final_collapse-unique.fastq" \
            &>> "${PATH_LOG_PA_ID}"

    fi

done

echo "Finished" &>> "${PATH_LOG}"
