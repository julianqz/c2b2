#!/usr/bin/env bash
# Create clonal consensus germlines after inferring B cell clones
#
# Author: Julian Q Zhou
# Date:   2021-05-09
#
# Prereqs:  
#
# 1) A CSV at ${PATH_CSV} containing two columns (NO HEADER)
#    1st col: subject ID; 2nd col: path to input tab file
# 
# 2) Input tab files should have a clone ID columns (`--cloned` flag is on)
# 
# 3) If novel alleles were inferred and found by tigger, they should have
#    been added to ${PATH_REFS}


# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -A  Path to CSV containing subject info."
    echo -e "  -B  Path to save outputs."
    echo -e "  -C  Path to IMGT reference fasta(s)."
    echo -e "  -D  Column name of sequence."
    echo -e "  -E  Column name of V call ."
    echo -e "  -F  Column name of D call."
    echo -e "  -G  Column name of J call."
    echo -e "  -H  Column name of clone ID."
    echo -e "  -I  Output foramt. Either 'airr' or 'changeo'."
    echo -e "  -h  This message."
}

# Get commandline arguments
while getopts "A:B:C:D:E:F:G:H:I:h" OPT; do
    case "$OPT" in
    A)  PATH_CSV=$(realpath "${OPTARG}")
        ;;
    B)  PATH_WORK=$(realpath "${OPTARG}")
		;;
	C)  PATH_REFS="${OPTARG}"
		;;
    D)  COL_SEQ="${OPTARG}"
        ;;
    E)  COL_V="${OPTARG}"
        ;;
    F)  COL_D="${OPTARG}"
        ;;
    G)  COL_J="${OPTARG}"
        ;;
    H)  COL_CLONE="${OPTARG}"
        ;;
    I)  CG_FORMAT="${OPTARG}"
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


# overall log for looping thru sample csv
PATH_LOG="${PATH_WORK}/log_bcr_createGermlines_${RUN_TYPE}_$(date '+%m%d%Y_%H%M%S').log"

CreateGermlines.py --version &> "${PATH_LOG}"
echo "Input csv: ${PATH_CSV}" &>> "${PATH_LOG}"
echo "References: " ${PATH_REFS} &>> "${PATH_LOG}"
echo "--vf: ${COL_V}" &>> "${PATH_LOG}"


N_LINES=$(wc -l < "${PATH_CSV}")
echo "N_LINES: ${N_LINES}" &>> "${PATH_LOG}"


for ((IDX=1; IDX<=${N_LINES}; IDX++)); do

	# read current line
	CUR_LINE=$(sed "${IDX}q;d" "${PATH_CSV}") 

	# split strings in unix
	# https://linuxhint.com/bash_split_examples/ 
	# following example 2

	# $IFS: internal field separator (default is white space)
	# -r: read backslash (\) as a character rather than escape character
	# -a: store split words into an array variable
	IFS=","
	read -a strarr <<< "${CUR_LINE}"

	# ID
	CUR_ID=${strarr[0]}

    # sample-specific path to input db
    PATH_INPUT=${strarr[1]}

    echo "IDX: ${IDX}; CUR_ID: ${CUR_ID}" &>> "${PATH_LOG}"

    #! important: if novel alleles were inferred and found by tigger, 
    #  add novel alleles to germline .fasta before running CreateGermlines

    # if you have run the clonal assignment task prior to invoking CreateGermlines, 
    # then adding the --cloned argument is recommended, as this will generate a 
    # single germline of consensus length for each clone

    # no nproc

    # appends _germ-pass.tab, _germ-fail.tab (even if input -d is .tsv)
    CreateGermlines.py \
        -d "${PATH_INPUT}" \
        -r ${PATH_REFS} \
        -g "full dmask vonly regions" \
        --sf "${COL_SEQ}" \
        --vf "${COL_V}" \
        --df "${COL_D}" \
        --jf "${COL_J}" \
        --cf "${COL_CLONE}" \
        --cloned \
        --format "${CG_FORMAT}" \
        --failed \
        --outname "${CUR_ID}" \
        --outdir "${PATH_WORK}" \
        --log "${PATH_WORK}/log_bcr_createGermlines_${CUR_ID}.log" \
        &>> "${PATH_LOG}"

done

echo "Finished" &>> "${PATH_LOG}"
