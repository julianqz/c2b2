#!/usr/bin/env bash
# Perform (initial) V(D)J gene annotation & parse IgBLAST results
#
# Author: Julian Q Zhou
# Date:   2021-04-22
#
# Prereqs:  
# 1) In ${PROJ_ID}/aux, a csv "input_fasta_${PROJ_ID}_${RUN_ID}" 
#    in which each row indicates the sample ID and the path to its input fasta

# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -J  Project ID."          
    echo -e "  -Z  Run ID. Usually one of {bulk, mAb, 10x}, but can be anything."              
    echo -e "  -T  Path to the top-level working dir." 
    echo -e "  -Y  Number of cores for parallelization." 
    echo -e "  -A  Whether to run IgBLAST. Boolean."
    echo -e "  -B  Whether to run MakeDb.py. Boolean."
    echo -e "  -C  Whether to split by chain. Boolean."
    echo -e "  -D  Whether to perform additional QC. Boolean."
    echo -e "  -E  Whether to split by productive vs. non-productive rearrangement. Boolean."
    echo -e "  -F  Path to script to split by chain."
    echo -e "  -G  Path to script to perform additional QC."
    echo -e "  -H  Path to IGDATA."
    echo -e "  -I  Path to igblastn."
    echo -e "  -K  Path to IMGT germline reference fastas."
    echo -e "  -h  This message."
}

# Get commandline arguments
while getopts "J:Z:T:Y:A:B:C:D:E:F:G:H:I:K:h" OPT; do
    case "$OPT" in
    J)  PROJ_ID="${OPTARG}"
        ;;
    Z)  RUN_ID="${OPTARG}"
		;;
    T)  PATH_ROOT=$(realpath "${OPTARG}")
        ;;
    Y)  NPROC="${OPTARG}"
        ;;
    A)  BOOL_IG="${OPTARG}"
        ;;
    B)  BOOL_MK="${OPTARG}"
        ;;
    C)  BOOL_SC="${OPTARG}"
        ;;
    D)  BOOL_QC="${OPTARG}"
        ;;
    E)  BOOL_SP="${OPTARG}"
        ;;
    F)  PATH_SCRIPT_SC=$(realpath "${OPTARG}")
        ;;
    G)  PATH_SCRIPT_QC=$(realpath "${OPTARG}")
        ;;
    H)  PATH_IGDATA=$(realpath "${OPTARG}")
        ;;
    I)  PATH_IGBLASTN=$(realpath "${OPTARG}")
        ;;
    K)  PATH_REF=$(realpath "${OPTARG}")
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

# output directory
PATH_OUTPUT="${PATH_PROJ}/data/annotate_${RUN_ID}/"
mkdir -p "${PATH_OUTPUT}"


# overall log for looping thru sample list
PATH_LOG="${PATH_AUX}log_bcr_annotate_${RUN_ID}_$(date '+%m%d%Y_%H%M%S').log"

NAME_LIST="input_fasta_${PROJ_ID}_${RUN_ID}.csv"
PATH_LIST="${PATH_AUX}${NAME_LIST}" 


MakeDb.py --version &> "${PATH_LOG}"
echo $("${PATH_IGBLASTN}" -version | head -n 1) &>> "${PATH_LOG}"
echo "NPROC=${NPROC}" &>> "${PATH_LOG}"
echo "Input fasta list: ${NAME_LIST}" &>> "${PATH_LOG}"
echo "BOOL_IG: ${BOOL_IG}" &>> "${PATH_LOG}"
echo "BOOL_MK: ${BOOL_MK}" &>> "${PATH_LOG}"
echo "BOOL_SC: ${BOOL_SC}" &>> "${PATH_LOG}"
echo "BOOL_QC: ${BOOL_QC}" &>> "${PATH_LOG}"
echo "BOOL_SP: ${BOOL_SP}" &>> "${PATH_LOG}"


N_LINES=$(wc -l < "${PATH_LIST}")
echo "N_LINES: ${N_LINES}" &>> "${PATH_LOG}"


for ((IDX=1; IDX<=${N_LINES}; IDX++)); do

	# read sample ID from file
	CUR_ID=$(sed "${IDX}q;d" "${PATH_LIST}") 

	echo "IDX: ${IDX}; CUR_ID: ${CUR_ID}" &>> "${PATH_LOG}"

	
    #* sample input fastq files
    RAW_1="${PATH_RAW}/${CUR_ID}_R1.fastq"
    RAW_2="${PATH_RAW}/${CUR_ID}_R2.fastq"

    
    # phix removal
    if $BOOL_PR; then

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
        # input files to presto-abseq pipeline are original files
        INPUT_ABSEQ_1="${RAW_1}"
        INPUT_ABSEQ_2="${RAW_2}"
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

        # convert fastq to fasta
        # creates ${CUR_ID}-final_collapse-unique_atleast-2.fasta
        cd "${PATH_OUTPUT_PA_ID}"

        echo "    Converting output fastq to fasta" &>> "${PATH_LOG}"

        "${PATH_SCRIPT_Q2A}" \
            "${PATH_OUTPUT_PA_ID}${CUR_ID}-final_collapse-unique_atleast-2.fastq" \
            &>> "${PATH_LOG_PA_ID}"

    fi

done

echo "Finished" &>> "${PATH_LOG}"
