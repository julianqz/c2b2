#!/usr/bin/env bash

# Author: Julian Q. Zhou
# https://github.com/julianqz
# Date:   2024-02-08
#
# Run `cellranger count` for a list of samples involving Feature Barcode
# Written for cellranger 7.2.0+
#
# Prereqs:  
# 1) The following must be in ${PROJ_ID}/aux/
#    - a sample list: "cr_list_count_${PROJ_ID}.txt"
#      each row is semi-colon-separated
#      [sample];[name of --libraries csv];[name of --feature-ref csv]
#    - all the --libraries csv's referenced in sample list
#    - all the --feature-ref csv's referenced in sample list


# NOTE: for data involving Feature Barcode

# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -J  Project ID."                        
    echo -e "  -T  Path to the top-level working dir." 
    echo -e "  -R  Path to the reference dir." 
    echo -e "  -Y  Number of cores for cellranger."    
    echo -e "  -Z  Amount of memory for cellranger."
    echo -e "  -W  Delete .bam* files. Default is true."   
    echo -e "  -h  This message."
}

PROJ_ID_SET=false
PATH_ROOT_SET=false
PATH_REF_SET=false
BOOL_DEL_BAM_SET=false

# Get commandline arguments
while getopts "J:T:R:Y:Z:W:h" OPT; do
    case "$OPT" in
    J)  PROJ_ID=$OPTARG
        PROJ_ID_SET=true
        ;;
    T)  PATH_ROOT=$(realpath $OPTARG)
        PATH_ROOT_SET=true
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


cd "${PATH_OUTPUT}"

for ((IDX=1; IDX<=${N_LINES}; IDX++)); do

    # read current line
    CUR_LINE=$(sed "${IDX}q;d" "${PATH_LIST}") 

    # use of ; (as opposed to ,) important for version of pipeline without fbc
    # keep using ; for consistency
    IFS=";"
    read -a strarr <<< "${CUR_LINE}"

	# sample ID
	CUR_ID=${strarr[0]}

    # --libraries csv
    CUR_CSV_LIB=${strarr[1]}
    PATH_CUR_CSV_LIB="${PATH_AUX}${CUR_CSV_LIB}"

    ## --feature-ref csv
    CUR_CSV_FR=${strarr[2]}
    PATH_CUR_CSV_FR="${PATH_AUX}${CUR_CSV_FR}"

    echo "IDX: ${IDX}; CUR_ID: ${CUR_ID}" &>> "${PATH_LOG}"
    echo "    --libraries csv: ${PATH_CUR_CSV_LIB}" &>> "${PATH_LOG}"
    echo "    --feature-ref csv: ${PATH_CUR_CSV_FR}" &>> "${PATH_LOG}"

	# sample-specific log
	PATH_LOG_ID="${PATH_AUX}log_cr_count_${IDX}_${CUR_ID}_$(date '+%m%d%Y_%H%M%S').log"

    # cellranger 7.2.0 doc

    # --fastqs cannot be used when performing Feature Barcode analysis; use --libraries instead.

    # --libraries
    # Path to a libraries.csv file declaring FASTQ paths and library types of input libraries. 
    # Required for gene expression + Feature Barcode analysis. 
    # When using --libraries, --fastqs and --sample must not be used. 
    # This argument should not be used when performing gene expression-only analysis; use --fastqs instead.

    # --feature-ref
    # Required for Feature Barcode analysis. 
    # Path to a Feature Reference CSV file declaring the Feature Barcode reagents used in the experiment.

    # csv examples
    # https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-feature-bc-analysis#libraries-csv

    # --include-introns is true by default
    # --no-bam is false by default

	cellranger count \
		--id "${CUR_ID}" \
        --libraries "${PATH_CUR_CSV_LIB}" \
        --feature-ref "${PATH_CUR_CSV_FR}" \
        --transcriptome "${PATH_REF}" \
        --nosecondary \
        --no-bam "${BOOL_DEL_BAM}" \
		--localcores "${CR_N}" \
		--localmem "${CR_M}" \
		&> "${PATH_LOG_ID}"


    rm "${PATH_OUTPUT}${CUR_ID}"/_*
    #rm -r "${PATH_OUTPUT}${CUR_ID}/SC_RNA_COUNTER_CS"
    #rm "${PATH_OUTPUT}${CUR_ID}/${CUR_ID}.mri.tgz"

done

echo "ALL DONE" &>> "${PATH_LOG}"
