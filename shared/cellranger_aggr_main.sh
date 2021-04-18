#!/usr/bin/env bash
# Run `cellranger aggr`
#
# Author: Julian Q Zhou
# Date:   2021-04-16
#
# Prereqs:  
# The following must be in ${PROJ_ID}/aux/
# - a project-specific config csv: "cr_aggr_${PROJ_ID}.csv"


# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -J  Project ID."                        #*
    echo -e "  -T  Path to the top-level working dir." #*
    echo -e "  -Y  Number of cores for cellranger."    #*
    echo -e "  -Z  Amount of memory for cellranger."   #*
    echo -e "  -h  This message."
}

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

PATH_PROJ="${PATH_ROOT}${PROJ_ID}/"
# no error if existing
mkdir -p "${PATH_PROJ}"

# csv, txt, logs
PATH_AUX="${PATH_PROJ}/aux/"
mkdir -p "${PATH_AUX}"

# cellranger aggr outputs
PATH_OUTPUT="${PATH_PROJ}/cr_aggr/"
mkdir -p "${PATH_OUTPUT}"

PATH_LOG="${PATH_AUX}log_cr_aggr_$(date '+%m%d%Y_%H%M%S').log"

NAME_CSV="cr_aggr_${PROJ_ID}.csv"
PATH_CSV="${PATH_AUX}${NAME_CSV}" 

AGGR_ID="${PROJ_ID}_aggr"

cellranger --version &> "${PATH_LOG}"
echo "--localcores=${CR_N}; --localmem=${CR_M}" &>> "${PATH_LOG}"
echo "Config csv: ${NAME_CSV}" &>> "${PATH_LOG}"

cd "${PATH_OUTPUT}"

cellranger aggr \
	--id "${AGGR_ID}" \
	--csv "${PATH_CSV}" \
    --normalize "mapped" \
    --nosecondary \
    --disable-ui \
	--localcores "${CR_N}" \
	--localmem "${CR_M}" \
	&> "${PATH_LOG}"

echo "ALL DONE" &>> "${PATH_LOG}"
