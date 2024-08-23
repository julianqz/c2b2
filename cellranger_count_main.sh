#!/usr/bin/env bash

# Author: Julian Q. Zhou
# https://github.com/julianqz
# Date:   2021-05-01
#
# Run `cellranger count` for a list of samples involving NO Feature Barcode
# 
# Prereqs:  
# 1) The following must be in ${PROJ_ID}/aux/
#    - a sample list: "cr_list_count_${PROJ_ID}.txt"
#      each row is semi-colon-separated
#      [sample];[comma-separated fastq id(s)]
# 2) Assumes that all fastqs are in one centralized folder

# NOTE: not for data involving Feature Barcode

# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -J  Project ID."                        
    echo -e "  -T  Path to the top-level working dir." 
    echo -e "  -F  Path to the centralized fastq dir." 
    echo -e "  -R  Path to the reference dir." 
    echo -e "  -Y  Number of cores for cellranger."    
    echo -e "  -Z  Amount of memory for cellranger."
    echo -e "  -W  Delete .bam* files. Default is true."
    echo -e "  -U  Input for --expect-cells. Optional."   
    echo -e "  -h  This message."
}

PROJ_ID_SET=false
PATH_ROOT_SET=false
PATH_FASTQ_SET=false
PATH_REF_SET=false
BOOL_DEL_BAM_SET=false

# Get commandline arguments
while getopts "J:T:F:R:Y:Z:W:U:h" OPT; do
    case "$OPT" in
    J)  PROJ_ID=$OPTARG
        PROJ_ID_SET=true
        ;;
    T)  PATH_ROOT=$(realpath $OPTARG)
        PATH_ROOT_SET=true
        ;;
    F)  PATH_FASTQ=$(realpath $OPTARG)
        PATH_FASTQ_SET=true
        ;;
    R)  PATH_REF=$(realpath $OPTARG)
        PATH_REF_SET=true
        ;;
    Y)  CR_N=$OPTARG
        ;;
    Z)  CR_M=$OPTARG
        ;;
    W)  BOOL_DEL_BAM=$OPTARG
        BOOL_DEL_BAM_SET=true
        ;;
    U)  NUM_EXP_CELLS=$OPTARG
        BOOL_NUM_EXP_CELLS_SET=true
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

# Exit if no project ID provided
if ! $PROJ_ID_SET; then
    echo "You must specify a project ID via the -J option" >&2
    exit 1
fi

# Exit if no top-level working directory provided
if ! $PATH_ROOT_SET; then
    echo "You must specify a top-level working directory that contains the project folder via the -T option" >&2
    exit 1
fi

# Exit if no centralized fastq directory provided
if ! $PATH_FASTQ_SET; then
    echo "You must specify a centralized fastq directory that contains all the fastq files via the -F option" >&2
    exit 1
fi

# Exit if no reference directory provided
if ! $PATH_REF_SET; then
    echo "You must specify a reference directory that contains the genome references via the -R option" >&2
    exit 1
fi

# Set BOOL_DEL_BAM to true if no -W specified
if ! $BOOL_DEL_BAM_SET; then
    BOOL_DEL_BAM=true
fi

# paths

PATH_PROJ="${PATH_ROOT}/${PROJ_ID}"
# no error if existing
mkdir -p "${PATH_PROJ}"

# csv, txt, logs
PATH_AUX="${PATH_PROJ}/aux/"
mkdir -p "${PATH_AUX}"

# cellranger count outputs
PATH_OUTPUT="${PATH_PROJ}/cr_count/"
mkdir -p "${PATH_OUTPUT}"

PATH_LOG="${PATH_AUX}log_cr_count_$(date '+%m%d%Y_%H%M%S').log"

NAME_LIST="cr_list_count_${PROJ_ID}.txt"
PATH_LIST="${PATH_AUX}${NAME_LIST}" 


cellranger --version &> "${PATH_LOG}"
echo "--localcores=${CR_N}; --localmem=${CR_M}" &>> "${PATH_LOG}"
echo "Sample list: ${NAME_LIST}" &>> "${PATH_LOG}"

N_LINES=$(wc -l < "${PATH_LIST}")
echo "N_LINES: ${N_LINES}" &>> "${PATH_LOG}"

if $BOOL_NUM_EXP_CELLS_SET; then
    echo "--expect-cells: ${NUM_EXP_CELLS}" &>> "${PATH_LOG}"
fi


cd "${PATH_OUTPUT}"

for ((IDX=1; IDX<=${N_LINES}; IDX++)); do

    # read current line
    CUR_LINE=$(sed "${IDX}q;d" "${PATH_LIST}") 

    # important for .txt to be semi-colon-separated
    # this allows for multiple fastq ids to be comma-separated for the same sample, if applicable
    IFS=";"
    read -a strarr <<< "${CUR_LINE}"

	# sample ID
	CUR_ID=${strarr[0]}

    # fastq id(s) (this is what cellranger calls `sample` -- confusing eh?)
    CUR_FASTQ_IDS=${strarr[1]}


	echo "IDX: ${IDX}; CUR_ID: ${CUR_ID}" &>> "${PATH_LOG}"

	# sample-specific log
	PATH_LOG_ID="${PATH_AUX}log_cr_count_${IDX}_${CUR_ID}_$(date '+%m%d%Y_%H%M%S').log"

    # default
    # --create-bam true
    # --include-introns true

    if $BOOL_NUM_EXP_CELLS_SET; then

        cellranger count \
            --id "${CUR_ID}" \
            --fastqs "${PATH_FASTQ}" \
            --sample "${CUR_FASTQ_IDS}" \
            --transcriptome "${PATH_REF}" \
            --nosecondary \
            --localcores "${CR_N}" \
            --localmem "${CR_M}" \
            --expect-cells "${NUM_EXP_CELLS}" \
            &> "${PATH_LOG_ID}"

    else

        cellranger count \
            --id "${CUR_ID}" \
            --fastqs "${PATH_FASTQ}" \
            --sample "${CUR_FASTQ_IDS}" \
            --transcriptome "${PATH_REF}" \
            --nosecondary \
            --localcores "${CR_N}" \
            --localmem "${CR_M}" \
            &> "${PATH_LOG_ID}"
            
    fi

    # remove .bam, .bambi, etc.
    
    if $BOOL_DEL_BAM; then
        rm "${PATH_OUTPUT}${CUR_ID}"/outs/*.bam*
    fi

    rm "${PATH_OUTPUT}${CUR_ID}"/_*
    rm -r "${PATH_OUTPUT}${CUR_ID}/SC_RNA_COUNTER_CS"
    rm "${PATH_OUTPUT}${CUR_ID}/${CUR_ID}.mri.tgz"

done

echo "ALL DONE" &>> "${PATH_LOG}"
