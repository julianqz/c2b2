#!/usr/bin/env bash

# Create clonal consensus germlines after inferring B cell clones
#
# Author: Julian Q Zhou
# Date:   2021-05-09
#
# Prereqs:  
#
# 1) A CSV at ${PATH_CSV} containing 2 or 3 columns (NO HEADER)
#    1st col: subject ID
#    2nd col: path to input tab file (heavy chain)
#    3rd col: path to input tab file (light chain) -- if `-M true`
# 
# 2) Input tab files should have a clone ID columns (`--cloned` flag is on)
# 
# 3) If novel alleles were inferred and found by tigger, they should be
#    supplied via `-F`


# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -A  Path to CSV containing subject info."
    echo -e "  -B  Path to save outputs."
    echo -e "  -C  Path to IMGT reference fasta for V."
    echo -e "  -D  Path to IMGT reference fasta for D."
    echo -e "  -E  Path to IMGT reference fasta for J."
    echo -e "  -F  Path to IMGT reference fasta for novel alleles (if any). Optional."
    echo -e "  -G  Column name of sequence."
    echo -e "  -H  Column name of V call ."
    echo -e "  -I  Column name of D call."
    echo -e "  -J  Column name of J call."
    echo -e "  -K  Column name of clone ID."
    echo -e "  -L  Output foramt. Either 'airr' or 'changeo'."
    echo -e "  -M  Whether to run for both heavy and light. Either 'true' or 'false'."
    echo -e "  -h  This message."
}

PATH_REF_N_SET=false

# Get commandline arguments
while getopts "A:B:C:D:E:F:G:H:I:J:K:L:M:h" OPT; do
    case "$OPT" in
    A)  PATH_CSV=$(realpath "${OPTARG}")
        ;;
    B)  PATH_WORK=$(realpath "${OPTARG}")
		;;
	C)  PATH_REF_V=$(realpath "${OPTARG}")
		;;
    D)  PATH_REF_D=$(realpath "${OPTARG}")
        ;;
    E)  PATH_REF_J=$(realpath "${OPTARG}")
        ;;
    F)  PATH_REF_N=$(realpath "${OPTARG}")
        PATH_REF_N_SET=true
        ;;
    G)  COL_SEQ="${OPTARG}"
        ;;
    H)  COL_V="${OPTARG}"
        ;;
    I)  COL_D="${OPTARG}"
        ;;
    J)  COL_J="${OPTARG}"
        ;;
    K)  COL_CLONE="${OPTARG}"
        ;;
    L)  CG_FORMAT="${OPTARG}"
        ;;
    M)  HEAVY_LIGHT="${OPTARG}"
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
PATH_LOG="${PATH_WORK}/log_bcr_createGermlines_$(date '+%m%d%Y_%H%M%S').log"

CreateGermlines.py --version &> "${PATH_LOG}"
echo "HEAVY_LIGHT: ${HEAVY_LIGHT}" &>> "${PATH_LOG}"
echo "Input csv: ${PATH_CSV}" &>> "${PATH_LOG}"
echo "Reference, V: ${PATH_REF_V}" &>> "${PATH_LOG}"
echo "Reference, D: ${PATH_REF_D}" &>> "${PATH_LOG}"
echo "Reference, J: ${PATH_REF_J}" &>> "${PATH_LOG}"
if $PATH_REF_N_SET; then
    echo "Reference, N: ${PATH_REF_N}" &>> "${PATH_LOG}"
else
    echo "Reference, N: not supplied" &>> "${PATH_LOG}"
fi
echo "--vf: ${COL_V}" &>> "${PATH_LOG}"


if $HEAVY_LIGHT; then
    RUN_MODE=(heavy light)
else
    RUN_MODE=(heavy)
fi


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

    echo "IDX: ${IDX}; CUR_ID: ${CUR_ID}" &>> "${PATH_LOG}"

    # loop thru heavy (and light)

    for CHAIN in "${RUN_MODE[@]}"; do

        # sample-specific path to input db
        if [[ ${CHAIN} == "heavy" ]]; then
            PATH_INPUT=${strarr[1]}
        else
            PATH_INPUT=${strarr[2]}
        fi

        echo " - ${CHAIN}: ${PATH_INPUT}" &>> "${PATH_LOG}"

        #! important: if novel alleles were inferred and found by tigger, 
        #  add novel alleles to germline .fasta before running CreateGermlines

        # if you have run the clonal assignment task prior to invoking CreateGermlines, 
        # then adding the --cloned argument is recommended, as this will generate a 
        # single germline of consensus length for each clone

        # no nproc

        # [outname]_germ-pass/fail.tsv

        if ! $PATH_REF_N_SET; then

            CreateGermlines.py \
                -d "${PATH_INPUT}" \
                -r "${PATH_REF_V}" "${PATH_REF_D}" "${PATH_REF_J}" \
                -g full dmask vonly regions \
                --sf "${COL_SEQ}" \
                --vf "${COL_V}" \
                --df "${COL_D}" \
                --jf "${COL_J}" \
                --cf "${COL_CLONE}" \
                --cloned \
                --format "${CG_FORMAT}" \
                --failed \
                --outname "${CUR_ID}_${CHAIN}" \
                --outdir "${PATH_WORK}" \
                --log "${PATH_WORK}/log_bcr_createGermlines_${CUR_ID}_${CHAIN}.log" \
                &>> "${PATH_LOG}"

        else

            CreateGermlines.py \
                -d "${PATH_INPUT}" \
                -r "${PATH_REF_V}" "${PATH_REF_D}" "${PATH_REF_J}" "${PATH_REF_N}" \
                -g full dmask vonly regions \
                --sf "${COL_SEQ}" \
                --vf "${COL_V}" \
                --df "${COL_D}" \
                --jf "${COL_J}" \
                --cf "${COL_CLONE}" \
                --cloned \
                --format "${CG_FORMAT}" \
                --failed \
                --outname "${CUR_ID}_${CHAIN}" \
                --outdir "${PATH_WORK}" \
                --log "${PATH_WORK}/log_bcr_createGermlines_${CUR_ID}_${CHAIN}.log" \
                &>> "${PATH_LOG}"

        fi

    done

done

echo "Finished" &>> "${PATH_LOG}"
