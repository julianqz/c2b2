#!/usr/bin/env bash
# Perform (initial) V(D)J gene annotation & parse IgBLAST results
#
# Author: Julian Q Zhou
# Date:   2021-04-22
#
# Prereqs:  
#
# 1) In ${PROJ_ID}/aux, a "bcr_input_annotate_${PROJ_ID}_${RUN_TYPE}.csv" 
#    in which each row indicates the sample ID and the path to its input fasta
#
# 2) If annotator is imgt, in ${MK_PATH_IMGT} (`-U`), imgt output .zip or .txz files 
#    with names that follow ${sample_id}${MK_IMGT_SUFFIX}, where
#    ${sample_id} matches the sample IDs in aux/input_fasta_${PROJ_ID}_${RUN_TYPE}, and
#    ${MK_IMGT_SUFFIX} is defined via `-V`
#
# 3) If run type is "10x" (`-B`) and 10x input is to be set for MakeDb.py (`-T`), 
#    the paths to sample-specific 10x csv/tsv are parsed from a third column in 
#    "bcr_input_annotate_${PROJ_ID}_${RUN_TYPE}.csv"
#    The filename of the 10x csv/tsv for each sample is assumed to be the same: 
#    "filtered_contig_annotations.csv" (generated by cellranger)
#
# Note:
# The Boolean run controls, with the exception of BOOL_QC, are not meant to alter the workflow.
# Rather, they are meant to control whether the workflow is executed in one setting.
# In other words, setting BOOL_IG to FALSE does not mean MakdeDb.py can run without 
#    the annotation step having been previously performed.
# The only exception is BOOL_QC. If FALSE, split_db will run on MakeDb.py output, 
#    instead of QC output.


# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -A  Project ID."          
    echo -e "  -B  Run type. One of {bulk, mAb, 10x}. Note the lowercase 'x'." 
    echo -e "  -C  Annotator. One of {igblast, imgt}."             
    echo -e "  -D  Path to the top-level working dir." 
    echo -e "  -E  Number of cores for parallelization." 
    echo -e "  -F  Whether to run IgBLAST. Boolean."
    echo -e "  -G  Whether to run MakeDb.py. Boolean."
    echo -e "  -H  Whether to perform QC. Boolean."
    echo -e "  -I  Whether to perform post-QC split. Boolean."
    echo -e "  -J  [IG] Path to IGDATA."
    echo -e "  -K  [IG] Path to igblastn."
    echo -e "  -L  [IG] Path to IMGT germline reference fastas."
    echo -e "  -M  [IG] Organism."
    echo -e "  -N  [IG] Loci. One of {ig, tr}."
    echo -e "  -O  [IG] --vdb. Name of custom V reference in IgBLAST database/."
    echo -e "  -P  [IG] --ddb. Name of custom D reference in IgBLAST database/."
    echo -e "  -Q  [IG] --jdb. Name of custom J reference in IgBLAST database/."
    echo -e "  -R  [IG] --format. One of {blast, airr}."
    echo -e "  -S  [MK] Whether to set the --partial flag. Boolean."
    echo -e "  -T  [MK] Whether to set the --10x flag. Boolean.\n" \
            "           If true, path to sample-specific 10x annotation csv/tsv is parsed (see prereqs)."
    echo -e "  -U  [MK] If annotator is 'imgt', path to IMGT/HighV-QUEST files."
    echo -e "  -V  [MK] If annotator is 'imgt', common suffix in IMGT/HighV-QUEST filenames.\n" \
            "           E.g. '_imgt.txz' or '_imgt.zip' in [sample_id]_imgt.txz or [sample_id]_imgt.zip respectively."
    echo -e "  -W  [MK] --format. One of {airr, changeo}."
	echo -e "  -X  [QCSP] Path to wrapper script to perform QC & split."
    echo -e "  -Y  [QCSP] Path to helper script to perform QC & split."
    echo -e "  -Z  [QCSP] --qcMaxPercN."
    echo -e "  -1  [QCSP] --qcColPercN. If multuple values, separate by comma.\n" \
            "             E.g. 'sequence_alignment, junction' "
    echo -e "  -2  [QCSP] --qcMaxNumNonATGCN."
    echo -e "  -3  [QCSP] --qcColNoneEmpty. If multuple values, separate by comma.\n" \
            "             E.g. 'germline_alignment, junction' "
    echo -e "  -4  [QCSP] --qcColNA. If multuple values, separate by comma.\n" \
            "             E.g. 'germline_alignment, junction, PRCONS' "
    echo -e "  -h  This message."
}

# Get commandline arguments
while getopts "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:1:2:3:4:h" OPT; do
    case "$OPT" in
    A)  PROJ_ID="${OPTARG}"
        ;;
    B)  RUN_TYPE="${OPTARG}"
		;;
	C)  ANNOTATOR="${OPTARG}"
		;;
    D)  PATH_ROOT=$(realpath "${OPTARG}")
        ;;
    E)  NPROC="${OPTARG}"
        ;;
    F)  BOOL_IG="${OPTARG}"
        ;;
    G)  BOOL_MK="${OPTARG}"
        ;;
    H)  BOOL_QC="${OPTARG}"
        ;;
    I)  BOOL_SP="${OPTARG}"
        ;;
    J)  PATH_IGDATA=$(realpath "${OPTARG}")
        ;;
    K)  PATH_IGBLASTN=$(realpath "${OPTARG}")
        ;;
    L)  PATH_REFS=$(realpath "${OPTARG}")
        ;;
    M)  IG_ORGANISM="${OPTARG}"
 		;;
 	N)  IG_LOCI="${OPTARG}"
 		;;
 	O)  IG_VDB="${OPTARG}"
 		;;
 	P)  IG_DDB="${OPTARG}"
 		;;
 	Q)  IG_JDB="${OPTARG}"
 		;;
 	R)  IG_FORMAT="${OPTARG}"
 		;;
 	S)  MK_PARTIAL="${OPTARG}"
 		;;
 	T)  MK_10X="${OPTARG}"
 		;;
 	U)  MK_PATH_IMGT=$(realpath "${OPTARG}")
 		;;
 	V)  MK_IMGT_SUFFIX="${OPTARG}"
 		;;
 	W)  MK_FORMAT="${OPTARG}"
 		;;
 	X)  PATH_SCRIPT_QCSP_WRAPPER=$(realpath "${OPTARG}")
		;;
	Y)  PATH_SCRIPT_QCSP_MAIN=$(realpath "${OPTARG}")
		;;
	Z)  QC_MAX_PERC_N="${OPTARG}"
		;;
	1)  QC_COL_PERC_N="${OPTARG}"
		;;
	2)  QC_MAX_NUM_NONATGCN="${OPTARG}"
		;;
	3)  QC_COL_NONE_EMPTY="${OPTARG}"
		;;
	4)  QC_COL_NA="${OPTARG}"
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


if $MK_PARTIAL; then
	PARTIAL="--partial"
else
	PARTIAL=""
fi

# paths

PATH_PROJ="${PATH_ROOT}/${PROJ_ID}"
# no error if existing
mkdir -p "${PATH_PROJ}"

# csv, txt, logs
PATH_AUX="${PATH_PROJ}/aux/"
mkdir -p "${PATH_AUX}"

# output directory
PATH_OUTPUT="${PATH_PROJ}/data/annotate_${RUN_TYPE}/"
mkdir -p "${PATH_OUTPUT}"


# overall log for looping thru sample csv
PATH_LOG="${PATH_AUX}log_bcr_annotate_${RUN_TYPE}_$(date '+%m%d%Y_%H%M%S').log"

NAME_CSV="bcr_input_annotate_${PROJ_ID}_${RUN_TYPE}.csv"
PATH_CSV="${PATH_AUX}${NAME_CSV}" 


MakeDb.py --version &> "${PATH_LOG}"
echo $("${PATH_IGBLASTN}" -version | head -n 1) &>> "${PATH_LOG}"
echo "NPROC=${NPROC}" &>> "${PATH_LOG}"
echo "Input csv: ${NAME_CSV}" &>> "${PATH_LOG}"
echo "ANNOTATOR: ${ANNOTATOR}" &>> "${PATH_LOG}"
echo "BOOL_IG: ${BOOL_IG}" &>> "${PATH_LOG}"
echo "BOOL_MK: ${BOOL_MK}" &>> "${PATH_LOG}"
echo "BOOL_QC: ${BOOL_QC}" &>> "${PATH_LOG}"
echo "BOOL_SP: ${BOOL_SP}" &>> "${PATH_LOG}"


N_LINES=$(wc -l < "${PATH_CSV}")
echo "N_LINES: ${N_LINES}" &>> "${PATH_LOG}"


for ((IDX=1; IDX<=${N_LINES}; IDX++)); do


	# input_annotate_${PROJ_ID}_${RUN_TYPE}.csv
	
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

    # sample-specific path to input fasta file
    PATH_INPUT_IG=${strarr[1]}

    # sample-specific path to 10x csv/tsv file    
    if [[ ${RUN_TYPE} == "10x" ]] && $MK_1OX; then
    	PATH_10X_ID=${strarr[2]}

    	# CSV_10X gets passed to MakeDb.py
    	CSV_10X="--10x ${PATH_10X_ID}filtered_contig_annotations.csv"
    else
    	CSV_10X=""
    fi

    echo "IDX: ${IDX}; CUR_ID: ${CUR_ID}" &>> "${PATH_LOG}"

    # sample-specific log
    PATH_LOG_ID="${PATH_AUX}log_annotation_${RUN_TYPE}_${IDX}_${CUR_ID}_$(date '+%m%d%Y_%H%M%S').log"
    
    # initiate log outside any if/else so that anything from inside if/else can use &>>
    echo "${PROJ_ID}_${RUN_TYPE}_${CUR_ID}" &> "${PATH_LOG_ID}"
    echo "${PATH_INPUT_IG}" &>> "${PATH_LOG_ID}"
    echo "${CSV_10X}" &>> "${PATH_LOG_ID}"

    # sample-specific output directory
    PATH_OUTPUT_ID="${PATH_OUTPUT}${CUR_ID}/"
    mkdir -p "${PATH_OUTPUT_ID}"


    # Run IgBLAST
    if $BOOL_IG; then
    	
    	echo "- running IgBLAST" &>> "${PATH_LOG}"

    	# output: [outname]_igblast.fmt7

    	AssignGenes.py \
    		--outdir "${PATH_OUTPUT_ID}" \
    		--outname "${CUR_ID}" \
    		--nproc "${NPROC}" \
    		-s "${PATH_INPUT_IG}" \
    		-b "${PATH_IGDATA}" \
    		--exec "${PATH_IGBLASTN}" \
    		--organism "${IG_ORGANISM}" \
    		--loci "${IG_LOCI}" \
    		--vdb "${IG_VDB}" \
    		--ddb "${IG_DDB}" \
    		--jdb "${IG_JDB}" \
    		--format "${IG_FORMAT}" \
    		&>> "${PATH_LOG_ID}"

    fi


    # Run MakeDb.py
    if $BOOL_MK; then

    	if [[ ${ANNOTATOR} == "igblast" ]]; then

    		PATH_ALIGN="${PATH_OUTPUT_ID}${CUR_ID}_igblast.fmt7"

    	elif [[ ${ANNOTATOR} == "imgt" ]]; then

    		PATH_ALIGN="${MK_PATH_IMGT}/${CUR_ID}${MK_IMGT_SUFFIX}"

    	fi

    	# assumes that igblast/imgt output exists
    	if [ -s "${PATH_ALIGN}" ]; then

    		echo "- running MakeDb.py; annotator is ${ANNOTATOR}" &>> "${PATH_LOG}"

    		# do not put "" around ${PATH_REFS} (otherwise * will be interpreted as is)
			# output: [outname]_db-pass.tab

			# if `false` passed to `MK_PARTIAL`, `PARTIAL` is set to empty
			# if `MK_10X` is set, `CSV_10X` is set to `--10x ${?}`; otherwise, empty

        	MakeDb.py "${ANNOTATOR}" \
        		--outdir "${PATH_OUTPUT_ID}" \
        		--outname "${CUR_ID}" \
        		--log "${PATH_OUTPUT_ID}log_makedb_${ANNOTATOR}_${CUR_ID}.log" \
        		--failed \
        		--extended \
        		"${PARTIAL}" \
        		--format "${MK_FORMAT}" \
        		-i "${PATH_ALIGN}" \
        		-r ${PATH_REFS} \
        		-s "${PATH_INPUT_IG}" \
        		"${CSV_10X}" \
        		&>> "${PATH_LOG_ID}"

        else
        	echo "${PATH_ALIGN} does not exist." &>> "${PATH_LOG_ID}"
        fi

    fi

    PATH_MK="${PATH_OUTPUT_ID}${CUR_ID}_db-pass.tab"


    # QC
    if $BOOL_QC; then

    	echo "- performing QC" &>> "${PATH_LOG}"
    	
    	# [outname]_qc.tsv

    	"${PATH_SCRIPT_QCSP_WRAPPER}" \
    		--helper "${PATH_SCRIPT_QCSP_MAIN}" \
    		--qc "TRUE" \
    		--qcSeq "${BOOL_QC_SEQ}" \
    		--qcCell "${BOOL_QC_CELL}" \
    		--qcDb "${PATH_MK}" \
    		--qcOutname "${CUR_ID}" \
    		--qcOutdir "${PATH_OUTPUT_ID}" \
    		--qcMaxPercN "${QC_MAX_PERC_N}" \
    		--qcColPercN "${QC_COL_PERC_N}" \
    		--qcMaxNumNonATGCN "${QC_MAX_NUM_NONATGCN}" \
    		--qcColNoneEmpty "${QC_COL_NONE_EMPTY}" \
    		--qcColNA "${QC_COL_NA}" \
    		--sp "FALSE" \
    		&>> "${PATH_LOG_ID}"

    	# set input name for split db
    	PATH_INPUT_SP="${CUR_ID}_qc.tsv"

    else
    	# set input name for split db
    	PATH_INPUT_SP="${PATH_MK}"
    fi


    # split
    if $BOOL_SP; then

    	echo "- splitting db" &>> "${PATH_LOG}"

    	# output: [outname]_[heavy|light]_[pr|npr].tsv

    	"${PATH_SCRIPT_QCSP_WRAPPER}" \
    		--helper "${PATH_SCRIPT_QCSP_MAIN}" \
    		--qc "FALSE" \
    		--sp "TRUE" \
    		--spDb "${PATH_INPUT_SP}" \
    		--spOutname "${CUR_ID} " \
    		--spOutdir "${PATH_OUTPUT_ID}" \
    		&>> "${PATH_LOG_ID}"
   
    fi

done

echo "Finished" &>> "${PATH_LOG}"
