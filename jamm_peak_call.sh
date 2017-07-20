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

pipeline="JAMM"

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
scripts=$HOME/Projects/PyATAC
gtf=$HOME/Brachypodium_distachyon/Phytozome/v3.1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf
fasta=$HOME/Brachypodium_distachyon/Phytozome/v3.1/assembly/Bdistachyon_314_v3.0.hardmasked.fa
JAMM=$HOME/usr/local/JAMM/JAMM.sh
noscaff=$HOME/Projects/PyATAC/Bd21.chrom.sizes.noscaff

# Call Peaks
# No_Pre filtering rmdups.bam
NOPRE=${pipeline}/No_Pre
# test ! -r $NOPRE &&  mkdir $NOPRE
echo "Creating ${NOPRE}/${OUTDIR}"
test ! -r ${NOPRE}/${OUTDIR} && mkdir -p ${NOPRE}/${OUTDIR}
find ${INPUT} -type f -name "*.rmdup.bam*" | xargs -I {} bash -c 'name=$(basename "{}");output='"$NOPRE/$OUTDIR"';output=$output/$name;ln -rs -T {} $output'
test ! -r ${NOPRE}/${OUTDIR}/input && mkdir ${NOPRE}/${OUTDIR}/input
bam2bed < ${NOPRE}/${OUTDIR}/${OUTDIR}.sorted.rmdup.bam > ${NOPRE}/${OUTDIR}/input/${OUTDIR}.sorted.rmdup.bed
prevdir=$(pwd)
cd ${NOPRE}/${OUTDIR}
(time bash $JAMM -s input/ -g $noscaff -o results -f1 -dy -p $THREADS | tee -a ../../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../../$LOG_ERR

# Peak annotation by HOMER
(time annotatePeaks.pl results/peaks/all.peaks.narrowPeak $fasta -gtf $gtf > results/peaks/${OUTDIR}_all_peaks_anno.txt | tee -a ../../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../../$LOG_ERR
(time annotatePeaks.pl results/peaks/filtered.peaks.narrowPeak $fasta -gtf $gtf > results/peaks/${OUTDIR}_filtered_peaks_anno.txt | tee -a ../../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../../$LOG_ERR

# Pre_filtering preshift75.bed from ${pipeline}
cd $prevdir
PRE=${pipeline}/Pre
# test ! -r $PRE && mkdir $PRE
echo "Creating ${PRE}/${OUTDIR}"
test ! -r ${PRE}/${OUTDIR} && mkdir -p ${PRE}/${OUTDIR}
test ! -r ${PRE}/${OUTDIR}/input && mkdir ${PRE}/${OUTDIR}/input
find ${INPUT} -type f -name "*preShift75.bed" | xargs -I {} bash -c 'name=$(basename "{}");ln -sr -T {} '"${PRE}/${OUTDIR}"'/input/${name}'
prevdir=$(pwd)
cd ${PRE}/${OUTDIR}
# -f1 -d y is for single-end ATAC-seq for JAMM. The JAMM documentation addresses this information
bash $JAMM -s input/ -g $noscaff -o results -f1 -d y -p $THREADS

# Peak annotation by HOMER
(time annotatePeaks.pl results/peaks/all.peaks.narrowPeak $fasta -gtf $gtf > results/peaks/${OUTDIR}_all_peaks_anno.txt | tee -a ../../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../../$LOG_ERR
(time annotatePeaks.pl results/peaks/filtered.peaks.narrowPeak $fasta -gtf $gtf > results/peaks/${OUTDIR}_filtered_peaks_anno.txt | tee -a ../../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../../$LOG_ERR
