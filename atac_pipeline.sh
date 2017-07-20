#!/bin/bash

if [[ $# -lt 3 ]]; then
    echo "$0: <Number of Threads> <BAM input file> <Output Dir/Filename> <Control> <Compressed>"
    exit 1
fi

echo "$1 is threads, $2 is input reads, $3 is output filename, $4 indicates use control or not, $5 indicates compressed or not"

control=$HOME/ENCODE/ATAC-seq/ATAC-pos_control-gDNA-0.1-gDNA-2/ATAC-pos_control-gDNA-0.1-gDNA-2.sorted.rmdup.bam

# base=$(basename "$3")
# outdir=$(dirname "$3")
base=${3}
outdir=${3}
echo "Creating orolab/${outdir}"
test ! -r orolab/${outdir} && mkdir orolab/${outdir}

LOG_ERR="$base.error.log"
LOG_OUT="$base.output.log"

if [ -f $LOG_ERR ];then
    rm $LOG_ERR
fi

if [ -f $LOG_OUT ];then
    rm $LOG_OUT
fi

touch orolab/${outdir}/$LOG_ERR
touch orolab/${outdir}/$LOG_OUT

## Checking if compressed and making named pipe
# if [ "$5" = true ]; then
#     pipeName=${2}.pipe
#     mkfifo ${pipeName}
#     gzip -dc ${2} > $pipeName &
#     INPUT=${pipeName}
# else
#     INPUT=${2}
# fi

# Locations of files that can be tweaked in pipeline.
index=/home/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/assembly/indices/Bdistachyon_314_v3.0.hardmasked
scripts=$HOME/ENCODE/Sarit_Reads/Scripts/PyATAC
gtf=$HOME/Brachypodium_distachyon/Phytozome/v3.1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf
fasta=$HOME/Brachypodium_distachyon/Phytozome/v3.1/assembly/Bdistachyon_314_v3.0.hardmasked.fa

# Alignment to reference
echo "bowtie alignment" `date`
if [ "$5" = true ]; then
    (time bowtie --chunkmbs 256 -p $1 -S -m 1 -X 2000 -t $index <(gzip -dc $2) | samtools view - -bS | samtools sort -o orolab/${outdir}/${3}.sorted.bam -T ${2}.pre | tee -a orolab/${outdir}/${LOG_OUT}) 3>&1 1>&2 2>&3 | tee -a orolab/${outdir}/${LOG_ERR}
else
    (time bowtie --chunkmbs 256 -p $1 -S -m 1 -X 2000 -t $index $2 | samtools view - -bS | samtools sort -o orolab/${outdir}/${3}.sorted.bam -T ${2}.pre | tee -a orolab/${outdir}/${LOG_OUT}) 3>&1 1>&2 2>&3 | tee -a orolab/${outdir}/${LOG_ERR}
fi
(time samtools index orolab/${outdir}/${3}.sorted.bam | tee -a orolab/${outdir}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a orolab/${outdir}/${LOG_ERR}

# Output only mapped reads.
echo "sam to bam (only output mapped)" `date`
(time samtools view -b -F 4 -o orolab/${outdir}/${3}.sorted.mapped.bam orolab/${outdir}/${3}.sorted.bam | tee -a orolab/${outdir}/${LOG_OUT}) 3>&1 1>&2 2>&3 | tee -a orolab/${outdir}/${LOG_ERR}
echo "No non sorted bam exists, don't delete"

# Remove duplicates.
echo "remove duplicates"
date
(time samtools rmdup orolab/${outdir}/${3}.sorted.mapped.bam orolab/${outdir}/${3}.sorted.rmdup.bam | tee -a orolab/${outdir}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a orolab/${outdir}/$LOG_ERR
(time samtools index orolab/${outdir}/${3}.sorted.rmdup.bam | tee -a orolab/${outdir}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a orolab/${outdir}/$LOG_ERR

# Make Fragment Bed
echo "Making fragment bed file"
(time perl $scripts/getFragmentBed.pl orolab/${outdir}/${3}.sorted.rmdup.bam orolab/${outdir}/${3}_fragment.bed | tee -a orolab/${outdir}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a orolab/${outdir}/$LOG_ERR
# Make UCSC Tracks
# Bd21 is ostensibly the database I'm accessing but it just needs the command for sanity checking (nothing?)
echo "make bedgraph" `date`
(time sort -k1,1 -k2,2n orolab/${outdir}/${3}_fragment.bed | bedItemOverlapCount Bd21 chromSize=$scripts/Bd21.chrom.sizes stdin | sort -k1,1 -k2,2n > orolab/${outdir}/${3}.bedGraph | tee -a orolab/${outdir}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a orolab/${outdir}/${LOG_ERR}

# Make bigWig file.
(time bedGraphToBigWig orolab/${outdir}/${3}.bedGraph $scripts/Bd21.chrom.sizes orolab/${outdir}/${3}.bw | tee -a orolab/${outdir}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a orolab/${outdir}/$LOG_ERR

# Call Peaks
# Before calling peaks (macs2) pre-shift the reads.
(time perl $scripts/preShift.pl orolab/${outdir}/${3}_fragment.bed 75 orolab/${outdir}/${3}_fragment_preShift75.bed | tee -a orolab/${outdir}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a orolab/${outdir}/$LOG_ERR
test ! -r orolab/${outdir}/callpeaks && mkdir orolab/${outdir}/callpeaks
cd orolab/${outdir}/callpeaks
# macs2 callpeak -t ../orolab/${outdir}/${3}_fragment_preShift75.bed -f BED -g 258400000 -n orolab/${outdir}/${3} --nomodel --shiftsize 75
if [ "$4" = true ]; then
    (time macs2 callpeak -t ../${base}_fragment_preShift75.bed -c ${control} -f BED -g 258400000 -n $base --nomodel --extsize 75 | tee -a ../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../$LOG_ERR
else
    (time macs2 callpeak -t ../${base}_fragment_preShift75.bed -f BED -g 258400000 -n $base --nomodel --extsize 75 | tee -a ../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../$LOG_ERR
fi

# Peak annotation by HOMER
# (time annotatePeaks.pl ${base}_peaks.narrowPeak $fasta -gtf $gtf > ${base}_peaks_anno.txt | tee -a ../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../$LOG_ERR
(time annotatePeaks.pl ${OUTDIR}_peaks.narrowPeak $fasta -gtf $gtf > ${base}_peaks_anno.txt | tee -a ../$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ../$LOG_ERR

