#!/bin/bash
usage="$(basename "$0") [-h] [-t Threads] [-i FastQ input file(s)] [-o Output directory/filename] [-c Compressed] -- runs JAMM ATAC-seq pipeline from FastQ files.

where:
	-h	show this help text
	-t	number of threads to use
	-i	FastQ input file(s)
	-o	Output directory
	-c	Input files, compressed?
	-f	filename of the error log
"

while getopts ':ht:i:o:cf:' option;do
    case "${option}" in
	h) echo "$usage"
	   exit
	   ;;
	t) THREADS=${OPTARG}
	   ;;
	i) INPUT=${OPTARG}
	   ;;
	o) OUTDIR=${OPTARG}
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

if [[ -z $THREADS ]] || [[ -z $INPUT ]] || [[ -z $OUTDIR ]]
then
    printf "Need threads, input, outdir, or all of the above"
    echo "$usage" >&2
    exit 1
fi

pipeline="orolab"
echo "Creating ${pipeline}/${OUTDIR}"
test ! -r ${pipeline}/${OUTDIR} &&  mkdir ${pipeline}/${OUTDIR}

LOG_ERR="${OUTDIR}.error.log"
LOG_OUT="${OUTDIR}.output.log"

if [ -f $LOG_ERR ];then
    echo "Removing previous error log"
    mv $LOG_ERR ${LOG_ERR}.bak
    # rm $LOG_ERR
fi

if [ -f $LOG_OUT ];then
    echo "Removing previous output log"
    mv $LOG_OUT ${LOG_OUT}.bak
    # rm $LOG_OUT
fi

touch ${pipeline}/${OUTDIR}/$LOG_ERR
touch ${pipeline}/${OUTDIR}/$LOG_OUT

# Locations of files that can be tweaked in pipeline.
index=/home/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/assembly/indices/Bdistachyon_314_v3.0.hardmasked
scripts=$HOME/ENCODE/Sarit_Reads/Scripts/PyATAC
gtf=$HOME/Brachypodium_distachyon/Phytozome/v3.1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf
fasta=$HOME/Brachypodium_distachyon/Phytozome/v3.1/assembly/Bdistachyon_314_v3.0.hardmasked.fa

# Alignment to reference
echo "bowtie alignment" `date`
if [ "COMPRESSED" = true ]; then
    (time bowtie --chunkmbs 256 -p $THREADS -S -m 1 -X 2000 -t $index <(gzip -dc $INPUT) | samtools view - -bS | samtools sort -o ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.bam -T ${2}.pre | tee -a ${pipeline}/${OUTDIR}/${LOG_OUT}) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/${LOG_ERR}
else
    (time bowtie --chunkmbs 256 -p $THREADS -S -m 1 -X 2000 -t $index $INPUT | samtools view - -bS | samtools sort -o ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.bam -T ${2}.pre | tee -a ${pipeline}/${OUTDIR}/${LOG_OUT}) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/${LOG_ERR}
fi
(time samtools index ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.bam | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/${LOG_ERR}

# Output only mapped reads.
echo "sam to bam (only output mapped)" `date`
(time samtools view -b -F 4 -o ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.mapped.bam ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.bam | tee -a ${pipeline}/${OUTDIR}/${LOG_OUT}) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/${LOG_ERR}
echo "No non sorted bam exists, don't delete"

# Remove duplicates.
echo "remove duplicates"
date
(time samtools rmdup ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.mapped.bam ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.rmdup.bam | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/$LOG_ERR
(time samtools index ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.rmdup.bam | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/$LOG_ERR

# Make Fragment Bed
echo "Making fragment bed file"
(time perl $scripts/getFragmentBed.pl ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.rmdup.bam ${pipeline}/${OUTDIR}/${OUTDIR}_fragment.bed | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/$LOG_ERR
# Make UCSC Tracks
# Bd21 is ostensibly the database I'm accessing but it just needs the command for sanity checking (nothing?)
echo "make bedgraph" `date`
(time sort -k1,1 -k2,2n ${pipeline}/${OUTDIR}/${OUTDIR}_fragment.bed | bedItemOverlapCount Bd21 chromSize=$scripts/Bd21.chrom.sizes stdin | sort -k1,1 -k2,2n > ${pipeline}/${OUTDIR}/${OUTDIR}.bedGraph | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/${LOG_ERR}

# Make bigWig file.
(time bedGraphToBigWig ${pipeline}/${OUTDIR}/${OUTDIR}.bedGraph $scripts/Bd21.chrom.sizes ${pipeline}/${OUTDIR}/${OUTDIR}.bw | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/$LOG_ERR

# Call Peaks
# Before calling peaks (macs2) pre-shift the reads.
(time perl $scripts/preShift.pl ${pipeline}/${OUTDIR}/${OUTDIR}_fragment.bed 75 ${pipeline}/${OUTDIR}/${OUTDIR}_fragment_preShift75.bed | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/$LOG_ERR
test ! -r ${pipeline}/${OUTDIR}/callpeaks && mkdir ${pipeline}/${OUTDIR}/callpeaks
cd ${pipeline}/${OUTDIR}/callpeaks
# macs2 callpeak -t ../${pipeline}/${OUTDIR}/${OUTDIR}_fragment_preShift75.bed -f BED -g 258400000 -n ${pipeline}/${OUTDIR}/${OUTDIR} --nomodel --shiftsize 75
(time macs2 callpeak -t ../${OUTDIR}_fragment_preShift75.bed -f BED -g 258400000 -n $base --nomodel --extsize 75 | tee -a ../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../$LOG_ERR

# Peak annotation by HOMER
# (time annotatePeaks.pl ${OUTDIR}_peaks.narrowPeak $fasta -gtf $gtf > ${OUTDIR}_peaks_anno.txt | tee -a ../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../$LOG_ERR
(time annotatePeaks.pl ${OUTDIR}_peaks.narrowPeak $fasta -gtf $gtf > ${OUTDIR}_peaks_anno.txt | tee -a ../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../$LOG_ERR

