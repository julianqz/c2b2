#!/usr/bin/env bash
# Run `cellranger multi` for a list of samples
#
# Author: Julian Q Zhou
# Date:   2021-04-15
#
# Prereqs:  
# The following must be in ${PROJ_ID}/aux/
# - a sample list: "cr_list_multi_${PROJ_ID}.txt"
# - sample-specific config csvs: "cr_multi_*.csv"


# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -J  Project ID."                        
    echo -e "  -T  Path to the top-level working dir." 
    echo -e "  -Y  Number of cores for cellranger."    
    echo -e "  -Z  Amount of memory for cellranger."   
    echo -e "  -h  This message."
}

PROJ_ID_SET=false
PATH_ROOT_SET=false

# Get commandline arguments
while getopts "J:T:Y:Z:h" OPT; do
    case "$OPT" in
    J)  PROJ_ID=$OPTARG
        PROJ_ID_SET=true
        ;;
    T)  PATH_ROOT=$(realpath $OPTARG)
        PATH_ROOT_SET=true
        ;;
    Y)  CR_N=$OPTARG
        ;;
    Z)  CR_M=$OPTARG
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


# paths

PATH_PROJ="${PATH_ROOT}/${PROJ_ID}"
# no error if existing
mkdir -p "${PATH_PROJ}"

# csv, txt, logs
PATH_AUX="${PATH_PROJ}/aux/"
mkdir -p "${PATH_AUX}"

# cellranger multi outputs
PATH_OUTPUT="${PATH_PROJ}/cr_multi/"
mkdir -p "${PATH_OUTPUT}"

PATH_LOG="${PATH_AUX}log_cr_multi_$(date '+%m%d%Y_%H%M%S').log"

NAME_LIST="cr_list_multi_${PROJ_ID}.txt"
PATH_LIST="${PATH_AUX}${NAME_LIST}" 


cellranger --version &> "${PATH_LOG}"
echo "--localcores=${CR_N}; --localmem=${CR_M}" &>> "${PATH_LOG}"
echo "Sample list: ${NAME_LIST}" &>> "${PATH_LOG}"

N_LINES=$(wc -l < "${PATH_LIST}")
echo "N_LINES: ${N_LINES}" &>> "${PATH_LOG}"


cd "${PATH_OUTPUT}"

for ((IDX=1; IDX<=${N_LINES}; IDX++)); do

	# read sample ID from file
	CUR_ID=$(sed "${IDX}q;d" "${PATH_LIST}") 

	echo "IDX: ${IDX}; CUR_ID: ${CUR_ID}" &>> "${PATH_LOG}"

	# sample-specific log
	PATH_LOG_ID="${PATH_AUX}log_cr_multi_${IDX}_${CUR_ID}_$(date '+%m%d%Y_%H%M%S').log"

    # sample-specific config csv
    NAME_CSV_ID="cr_multi_${CUR_ID}.csv"
    PATH_CSV_ID="${PATH_AUX}${NAME_CSV_ID}"

	cellranger multi \
		--id "${CUR_ID}" \
		--csv "${PATH_CSV_ID}" \
		--localcores "${CR_N}" \
		--localmem "${CR_M}" \
		&> "${PATH_LOG_ID}"

    # clean up
    rm "${PATH_OUTPUT}${CUR_ID}/"_*
    rm "${PATH_OUTPUT}${CUR_ID}/${CUR_ID}.mri.tgz"
    rm -r "${PATH_OUTPUT}${CUR_ID}/SC_MULTI_CS"

    # outs/multi/vdj_[bt]/*.bed, *.bam, *.bambai, *.fasta.fai
    # outs/per_sample_outs/${CUR_ID}/vdj_[bt]/*.bam, *.bambai, *.fasta.fai
    # still there despite specifying --no-bam in config csv
    rm "${PATH_OUTPUT}${CUR_ID}/outs/multi/vdj_[bt]/*.bed"
    rm "${PATH_OUTPUT}${CUR_ID}/outs/multi/vdj_[bt]/*.[bf]a[im]"
    rm "${PATH_OUTPUT}${CUR_ID}/outs/per_sample_outs/${CUR_ID}/vdj_[bt]/*.[bf]a[im]"

done

echo "ALL DONE" &>> "${PATH_LOG}"
