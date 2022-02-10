#!/usr/bin/env bash

# Author: Julian Q. Zhou
# https://github.com/julianqz
# Date:   2021-04-22
#
# Perform (initial) V(D)J gene annotation & parse IgBLAST results
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
#    the paths to the directory containing the sample-specific 10x csv/tsv are 
#    parsed from a third column in "bcr_input_annotate_${PROJ_ID}_${RUN_TYPE}.csv"
#    Note that the path is to the directory, and should not end in .fasta.
#    The filename of the 10x csv/tsv for each sample is assumed to be the same: 
#    "filtered_contig_annotations.csv" (generated by cellranger)
#
# 4) It is implicitly assumed that concat_no_dup_[VDJ].fasta exist in the 
#    IMGT reference folder (in imgt_select/), which get passed to `-r` of MakeDb.py.
#    This should be the case as long as this script is used with julianqz/wu_cimm:ref_[ver].
#    Tried putting the regex in PATH_REFS and passing PATH_REFS to CMD in the bsub script
#    ==> expansion appears to get messed up
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
    echo -e "  -B  Run type. One of {bulk, mab, nested, 10x}. Note the lowercase 'x'." 
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
    echo -e "  -M  [IG] -organism. {human, mouse, etc.}"
    echo -e "  -N  [IG] -ig_seqtype. One of {Ig, TCR}."
    echo -e "  -O  [IG] -germline_db_V. Name of custom V reference in IgBLAST database/."
    echo -e "  -P  [IG] -germline_db_D. Name of custom D reference in IgBLAST database/."
    echo -e "  -Q  [IG] -germline_db_J. Name of custom J reference in IgBLAST database/."
    echo -e "  -r  [IG] Whether to annotate C region. One of {true, false}."
    echo -e "  -s  [IG] -c_region_db. Name of custom C reference in IgBLAST database/."
    echo -e "  -R  [IG] Determines -outfmt. One of {blast, airr}."
    echo -e "  -T  [MK] Whether to set the --10x flag. Boolean.\n" \
            "           If true, path to sample-specific 10x annotation csv/tsv is parsed (see prereqs)."
    echo -e "  -U  [MK] If annotator is 'imgt', path to IMGT/HighV-QUEST files."
    echo -e "  -V  [MK] If annotator is 'imgt', common suffix in IMGT/HighV-QUEST filenames.\n" \
            "           E.g. '_imgt.txz' or '_imgt.zip' in [sample_id]_imgt.txz or [sample_id]_imgt.zip respectively."
    echo -e "  -W  [MK] --format. One of {airr, changeo}."
	echo -e "  -X  [QCSP] Path to wrapper script to perform QC & split."
    echo -e "  -Y  [QCSP] Path to helper script to perform QC & split."
    echo -e "  -Z  [QCSP] --qcSeq. Boolean for R."
    echo -e "  -1  [QCSP] --qcCell. Boolean for R."
    echo -e "  -2  [QCSP] --qcSequential. Boolean for R."
    echo -e "  -3  [QCSP] --qcColV."
    echo -e "  -4  [QCSP] --qcColD."
    echo -e "  -5  [QCSP] --qcColJ."
    echo -e "  -6  [QCSP] --qcColC."
    echo -e "  -7  [QCSP] --qcColObsv."
    echo -e "  -8  [QCSP] --qcColGerm."
    echo -e "  -9  [QCSP] --qcMaxN."
    echo -e "  -a  [QCSP] --qcColN. If multuple values, separate by comma.\n" \
            "             E.g. 'sequence_alignment, cdr3' "
    echo -e "  -b  [QCSP] --qcLastPosN. If multiple values, separate by comma.\n" \
            "             E.g. '12, 312' "
    echo -e "  -c  [QCSP] --qcAsPercN. Boolean for R."        
    echo -e "  -d  [QCSP] --qcMaxNonATGC."
    echo -e "  -e  [QCSP] --qcAsPercNonATGC. Boolean for R."
    echo -e "  -f  [QCSP] --qcColNoneEmpty. If multuple values, separate by comma.\n" \
            "             E.g. 'germline_alignment, cdr3' "
    echo -e "  -g  [QCSP] --qcColNA. If multuple values, separate by comma.\n" \
            "             E.g. 'germline_alignment, cdr3, PRCONS' "
	echo -e "  -i  [QCSP] --qcColLenMod3."
    echo -e "  -j  [QCSP] --qcColLocus."
    echo -e "  -k  [QCSP] --qcColCell."
    echo -e "  -m  [QCSP] --qcColUMI."
    echo -e "  -n  [QCSP] --qcLogicNumHL."
	echo -e "  -o  [QCSP] --spColV."
	echo -e "  -p  [QCSP] --spColProd."            
	echo -e "  -q  [QCSP] --spValProd."
    echo -e "  -h  This message."
}

# Get commandline arguments
while getopts "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:r:s:R:T:U:V:W:X:Y:Z:1:2:3:4:5:6:7:8:9:a:b:c:d:e:f:g:i:j:k:m:n:o:p:q:h" OPT; do
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
    L)  PATH_REFS="${OPTARG}"
        ;;
    M)  IG_ORGANISM="${OPTARG}"
		;;
	N)  IG_SEQTYPE="${OPTARG}"
 		;;
 	O)  IG_VDB="${OPTARG}"
 		;;
 	P)  IG_DDB="${OPTARG}"
 		;;
 	Q)  IG_JDB="${OPTARG}"
 		;;
    r)  BOOL_ALIGN_C="${OPTARG}"
        ;;
    s)  IG_CDB="${OPTARG}"
        ;;
 	R)  IG_FORMAT="${OPTARG}"
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
	Z)  BOOL_QC_SEQ="${OPTARG}"
		;;
	1)  BOOL_QC_CELL="${OPTARG}"
		;;
    2)  BOOL_QC_SEQUENTIAL="${OPTARG}"
        ;;
	3)  QC_COL_V="${OPTARG}"
		;;
	4)  QC_COL_D="${OPTARG}"
		;;
	5)  QC_COL_J="${OPTARG}"
		;;
	6)  QC_COL_C="${OPTARG}"
		;;
	7)  QC_COL_OBSV="${OPTARG}"
		;;
	8)  QC_COL_GERM="${OPTARG}"
		;;	
	9)  QC_MAX_N="${OPTARG}"
		;;
	a)  QC_COL_N="${OPTARG}"
		;;
    b)  QC_LAST_POS_N="${OPTARG}"
        ;;
    c)  BOOL_AS_PERC_N="${OPTARG}"
        ;;
	d)  QC_MAX_NONATGC="${OPTARG}"
		;;
    e)  BOOL_AS_PERC_NONATGC="${OPTARG}"
        ;;
	f)  QC_COL_NONE_EMPTY="${OPTARG}"
		;;
	g)  QC_COL_NA="${OPTARG}"
		;;
	i)  QC_COL_LEN_MOD3="${OPTARG}"
		;;
    j)  QC_COL_LOCUS="${OPTARG}"
        ;;
    k)  QC_COL_CELL="${OPTARG}"
        ;;
    m)  QC_COL_UMI="${OPTARG}"
        ;;
    n)  QC_LOGIC_NUM_HL="${OPTARG}"
        ;;
	o)  SP_COL_V="${OPTARG}"
		;;
	p)  SP_COL_PROD="${OPTARG}"
		;;
	q)  SP_VAL_PROD="${OPTARG}"
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
echo "BOOL_ALIGN_C: ${BOOL_ALIGN_C}" &>> "${PATH_LOG}"
echo "BOOL_MK: ${BOOL_MK}" &>> "${PATH_LOG}"
echo "BOOL_QC: ${BOOL_QC}" &>> "${PATH_LOG}"
echo "BOOL_SP: ${BOOL_SP}" &>> "${PATH_LOG}"

echo "${IG_ORGANISM}" &>> "${PATH_LOG}"

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

    	# CSV_10X gets passed to `--10x` of MakeDb.py
    	CSV_10X="${PATH_10X_ID}/filtered_contig_annotations.csv"
    fi

    echo "IDX: ${IDX}; CUR_ID: ${CUR_ID}" &>> "${PATH_LOG}"

    # sample-specific log
    PATH_LOG_ID="${PATH_AUX}log_bcr_annotate_${RUN_TYPE}_${IDX}_${CUR_ID}_$(date '+%m%d%Y_%H%M%S').log"
    
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


        if [[ ${IG_FORMAT} == "blast" ]]; then
            OUTFMT="7 std qseq sseq btop"
        elif [[ ${IG_FORMAT} == "airr" ]]; then
            OUTFMT="19"
        fi

        echo "IG_FORMAT: ${IG_FORMAT}" &>> "${PATH_LOG_ID}"
        echo "-outfmt: ${OUTFMT}" &>> "${PATH_LOG_ID}"

        export IGDATA="${PATH_IGDATA}"
        export BLASTDB="${PATH_IGDATA}/database"

        # output: [outname]_igblast.fmt7

        if $BOOL_ALIGN_C; then

            "${PATH_IGBLASTN}" \
                -query "${PATH_INPUT_IG}" \
                -out "${PATH_OUTPUT_ID}${CUR_ID}_igblast.fmt7" \
                -num_threads "${NPROC}" \
                -ig_seqtype "${IG_SEQTYPE}" \
                -organism "${IG_ORGANISM}" \
                -auxiliary_data "${PATH_IGDATA}/optional_file/${IG_ORGANISM}_gl.aux" \
                -germline_db_V "${IG_VDB}" \
                -germline_db_D "${IG_DDB}" \
                -germline_db_J "${IG_JDB}" \
                -c_region_db "${IG_CDB}" \
                -outfmt "${OUTFMT}" \
                -domain_system "imgt" \
                &>> "${PATH_LOG_ID}"

        else

            "${PATH_IGBLASTN}" \
                -query "${PATH_INPUT_IG}" \
                -out "${PATH_OUTPUT_ID}${CUR_ID}_igblast.fmt7" \
                -num_threads "${NPROC}" \
                -ig_seqtype "${IG_SEQTYPE}" \
                -organism "${IG_ORGANISM}" \
                -auxiliary_data "${PATH_IGDATA}" \
                -germline_db_V "${IG_VDB}" \
                -germline_db_D "${IG_DDB}" \
                -germline_db_J "${IG_JDB}" \
                -outfmt "${OUTFMT}" \
                -domain_system "imgt" \
                &>> "${PATH_LOG_ID}"

        fi


    	#AssignGenes.py igblast \
    	#	--outdir "${PATH_OUTPUT_ID}" \
    	#	--outname "${CUR_ID}" \
    	#	--nproc "${NPROC}" \
    	#	-s "${PATH_INPUT_IG}" \
    	#	-b "${PATH_IGDATA}" \
    	#	--exec "${PATH_IGBLASTN}" \
    	#	--organism "${IG_ORGANISM}" \
    	#	--loci "${IG_LOCI}" \
    	#	--vdb "${IG_VDB}" \
    	#	--ddb "${IG_DDB}" \
    	#	--jdb "${IG_JDB}" \
    	#	--format "${IG_FORMAT}" \
    	#	&>> "${PATH_LOG_ID}"

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

    		# do not put "" around * (otherwise * will be interpreted as is)
			# output: [outname]_db-pass.tsv

			# if `false` passed to `MK_PARTIAL`, `PARTIAL` is set to empty
			# if `MK_10X` is set, `CSV_10X` is set to `--10x ${?}`; otherwise, empty

            # use --extended so that fields such as "cdr3" are created

			if $MK_10X; then

	        	MakeDb.py "${ANNOTATOR}" \
	        		--outdir "${PATH_OUTPUT_ID}" \
	        		--outname "${CUR_ID}" \
	        		--log "${PATH_OUTPUT_ID}log_makedb_${ANNOTATOR}_${CUR_ID}.log" \
	        		--failed \
	        		--extended \
	        		--format "${MK_FORMAT}" \
	        		-i "${PATH_ALIGN}" \
	        		-r "${PATH_REFS}concat_no_dup_"*.fasta \
	        		-s "${PATH_INPUT_IG}" \
	        		--10x "${CSV_10X}" \
	        		&>> "${PATH_LOG_ID}"

	        else

	        	MakeDb.py "${ANNOTATOR}" \
	        		--outdir "${PATH_OUTPUT_ID}" \
	        		--outname "${CUR_ID}" \
	        		--log "${PATH_OUTPUT_ID}log_makedb_${ANNOTATOR}_${CUR_ID}.log" \
	        		--failed \
	        		--extended \
	        		--format "${MK_FORMAT}" \
	        		-i "${PATH_ALIGN}" \
	        		-r "${PATH_REFS}concat_no_dup_"*.fasta \
	        		-s "${PATH_INPUT_IG}" \
	        		&>> "${PATH_LOG_ID}"

        	fi

        else
        	echo "${PATH_ALIGN} does not exist." &>> "${PATH_LOG_ID}"
        fi

    fi

    PATH_MK="${PATH_OUTPUT_ID}${CUR_ID}_db-pass.tsv"


    # QC
    if $BOOL_QC; then

    	echo "- performing QC" &>> "${PATH_LOG}"
    	
    	# [outname]_qc.tsv

    	"${PATH_SCRIPT_QCSP_WRAPPER}" \
    		--helper "${PATH_SCRIPT_QCSP_MAIN}" \
    		--qc "TRUE" \
    		--qcSeq "${BOOL_QC_SEQ}" \
    		--qcCell "${BOOL_QC_CELL}" \
            --qcSequential "${BOOL_QC_SEQUENTIAL}" \
    		--qcDb "${PATH_MK}" \
    		--qcOutname "${CUR_ID}" \
    		--qcOutdir "${PATH_OUTPUT_ID}" \
    		--qcColV "${QC_COL_V}" \
    		--qcColD "${QC_COL_D}" \
    		--qcColJ "${QC_COL_J}" \
    		--qcColC "${QC_COL_C}" \
    		--qcColObsv "${QC_COL_OBSV}" \
    		--qcColGerm "${QC_COL_GERM}" \
    		--qcMaxN "${QC_MAX_N}" \
    		--qcColN "${QC_COL_N}" \
            --qcLastPosN "${QC_LAST_POS_N}" \
            --qcAsPercN "${BOOL_AS_PERC_N}" \
    		--qcMaxNonATGC "${QC_MAX_NONATGC}" \
            --qcAsPercNonATGC "${BOOL_AS_PERC_NONATGC}" \
    		--qcColNoneEmpty "${QC_COL_NONE_EMPTY}" \
    		--qcColNA "${QC_COL_NA}" \
    		--qcColLenMod3 "${QC_COL_LEN_MOD3}" \
            --qcColLocus "${QC_COL_LOCUS}" \
            --qcColCell "${QC_COL_CELL}" \
            --qcColUMI "${QC_COL_UMI}" \
            --qcLogicNumHL "${QC_LOGIC_NUM_HL}" \
    		--sp "FALSE" \
    		&>> "${PATH_LOG_ID}"

    	# set input name for split db
    	PATH_INPUT_SP="${PATH_OUTPUT_ID}${CUR_ID}_qc-pass.tsv"

    else
    	# set input name for split db
    	PATH_INPUT_SP="${PATH_MK}"
    fi


    # split
    if $BOOL_SP; then

        if [ -s "${PATH_INPUT_SP}" ]; then

        	echo "- splitting db" &>> "${PATH_LOG}"

        	# output: [outname]_[heavy|light]_[pr|npr].tsv

        	"${PATH_SCRIPT_QCSP_WRAPPER}" \
        		--helper "${PATH_SCRIPT_QCSP_MAIN}" \
        		--qc "FALSE" \
        		--sp "TRUE" \
        		--spDb "${PATH_INPUT_SP}" \
        		--spOutname "${CUR_ID}_qc-pass" \
        		--spOutdir "${PATH_OUTPUT_ID}" \
        		--spColV "${SP_COL_V}" \
        		--spColProd "${SP_COL_PROD}" \
        		--spValProd "${SP_VAL_PROD}" \
        		&>> "${PATH_LOG_ID}"

        else
            echo "- splitting db skipped (${PATH_INPUT_SP} does not exist)" &>> "${PATH_LOG}"
        fi
   
    fi

done

echo "Finished" &>> "${PATH_LOG}"
