#!/bin/bash
usage="$(basename "$0") [-h] [-t Threads] [-i Directory of pre_preak_call.sh output] [-c Compressed] -- runs JAMM ATAC-seq pipeline from FastQ files.

where:
	-h	show this help text
	-t	number of threads to use
	-i	Directory of pre_peak_call.sh output for FastQ
	-c	Input files, compressed?
	-f	filename of the error log
"

while getopts ':ht:i:cf:' option;do
    case "${option}" in
	h) echo "$usage"
	   exit
	   ;;
	t) THREADS=${OPTARG}
	   ;;
	i) INPUT=${OPTARG}
	   ;;
	c) COMPRESSED=true
	   ;;
	f) ERROR_LOG=${OPTARG}
	   ;;
	:) printf "missing argument for -%s\n" "$OPTARG" >&2
	   echo "$usage" >&2
	   exit 1
	   ;;
	\?) printf "illegal option: -%s\n" "$OPTARG" >&2
	    echo "$usage" >&2
	    exit 1
	    ;;
    esac
done
shift $((OPTIND - 1))

if [[ -z $THREADS ]] || [[ -z $INPUT ]]
then
    printf "Need threads, input, outdir, or all of the above"
    echo "$usage" >&2
    exit 1
fi

pipeline="orolab"

OUTDIR=$(basename $INPUT)
PREPEAK=$(dirname $INPUT)

echo "Creating ${pipeline}/"
test ! -r ${pipeline}/ && mkdir ${pipeline}

LOG_ERR="${OUTDIR}.error.log"
LOG_OUT="${OUTDIR}.output.log"

if [ -f ${pipeline}/${LOG_ERR} ];then
    echo "Removing previous error log"
    mv ${pipeline}/${LOG_ERR} ${pipeline}/${LOG_ERR}.bak
    # rm $LOG_ERR
fi

if [ -f ${pipeline}/${LOG_OUT} ];then
    echo "Removing previous output log"
    mv ${pipeline}/${LOG_OUT} ${pipeline}/${LOG_OUT}.bak
    # rm $LOG_OUT
fi

touch ${pipeline}/$LOG_ERR
touch ${pipeline}/$LOG_OUT

# Locations of files that can be tweaked in pipeline.
index=/home/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/assembly/indices/Bdistachyon_314_v3.0.hardmasked
scripts=$HOME/ENCODE/Sarit_Reads/Scripts/PyATAC
gtf=$HOME/Brachypodium_distachyon/Phytozome/v3.1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf
fasta=$HOME/Brachypodium_distachyon/Phytozome/v3.1/assembly/Bdistachyon_314_v3.0.hardmasked.fa

# Call Peaks
test ! -r ${pipeline}/${OUTDIR}/callpeaks && mkdir -p ${pipeline}/${OUTDIR}/callpeaks
find ${INPUT} -type f -name "*preShift75.bed" | xargs -I {} bash -c 'name=$(basename "{}");ln -sr -T {} '"${pipeline}/${OUTDIR}"'/${name}'
prevdir=$(pwd)
cd ${pipeline}/${OUTDIR}/callpeaks
# macs2 callpeak -t ../${pipeline}/${OUTDIR}/${OUTDIR}_fragment_preShift75.bed -f BED -g 258400000 -n ${pipeline}/${OUTDIR}/${OUTDIR} --nomodel --shiftsize 75
(time macs2 callpeak -t ../${OUTDIR}_fragment_preShift75.bed -f BED -g 258400000 -n ${OUTDIR} --nomodel --extsize 75 | tee -a ../../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../../$LOG_ERR

# Peak annotation by HOMER
# (time annotatePeaks.pl ${OUTDIR}_peaks.narrowPeak $fasta -gtf $gtf > ${OUTDIR}_peaks_anno.txt | tee -a ../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../$LOG_ERR
(time annotatePeaks.pl ${OUTDIR}_peaks.narrowPeak $fasta -gtf $gtf > ${OUTDIR}_peaks_anno.txt | tee -a ../../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../../$LOG_ERR
cd $prevdir

