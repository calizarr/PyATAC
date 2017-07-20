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

pipeline="Pre_Peak"

OUTDIR=$(basename $OUTDIR)
OUTDIR=$(echo $OUTDIR | sed 's/\.fastq//g' | sed 's/\.gz//g')
echo "Creating ${pipeline}/${OUTDIR}"
test ! -r ${pipeline}/${OUTDIR} &&  mkdir -p ${pipeline}/${OUTDIR}

LOG_ERR="${OUTDIR}.error.log"
LOG_OUT="${OUTDIR}.output.log"

if [ -f ${pipeline}/${OUTDIR}/${LOG_ERR} ];then
    echo "Removing previous error log"
    mv ${pipeline}/${OUTDIR}/${LOG_ERR} ${pipeline}/${OUTDIR}/${LOG_ERR}.bak
    # rm $LOG_ERR
fi

if [ -f ${pipeline}/${OUTDIR}/${LOG_OUT} ];then
    echo "Removing previous output log"
    mv ${pipeline}/${OUTDIR}/${LOG_OUT} ${pipeline}/${OUTDIR}/${LOG_OUT}.bak
    # rm $LOG_OUT
fi

touch ${pipeline}/${OUTDIR}/$LOG_ERR
touch ${pipeline}/${OUTDIR}/$LOG_OUT

# Locations of files that can be tweaked in pipeline.
index=/home/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/assembly/indices/Bdistachyon_314_v3.0.hardmasked
scripts=$HOME/Projects/PyATAC
gtf=$HOME/Brachypodium_distachyon/Phytozome/v3.1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf
fasta=$HOME/Brachypodium_distachyon/Phytozome/v3.1/assembly/Bdistachyon_314_v3.0.hardmasked.fa
noscaff=$HOME/Projects/PyATAC/Bd21.chrom.sizes.noscaff

# Alignment to reference genome
echo "bowtie alignment" `date`
if [ -f ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.bam ]
then
    echo "File: ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.bam has already been created. Moving on"
else
    if [ "$COMPRESSED" = true ]; then
	(time bowtie --chunkmbs 512 -p $THREADS -S -m 1 -X 2000 -t $index <(gzip -dc $INPUT) | samtools view - -bS | samtools sort -o ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.bam -T ${INPUT}.pre | tee -a ${pipeline}/${OUTDIR}/${LOG_OUT}) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/${LOG_ERR}
    else
	(time bowtie --chunkmbs 512 -p $THREADS -S -m 1 -X 2000 -t $index $INPUT | samtools view - -bS | samtools sort -o ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.bam -T ${INPUT}.pre | tee -a ${pipeline}/${OUTDIR}/${LOG_OUT}) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/${LOG_ERR}
    fi
    (time samtools index ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.bam | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/${LOG_ERR}
fi

# Output only mapped reads.
if [ -f ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.mapped.bam ]
then
    echo "Mapped file: ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.mapped.bam already exists. Moving on..."
else
    echo "sam to bam (only output mapped)" `date`
    (time samtools view -b -F 4 -o ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.mapped.bam ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.bam | tee -a ${pipeline}/${OUTDIR}/${LOG_OUT}) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/${LOG_ERR}
    echo "No non sorted bam exists, don't delete"
fi

# Remove duplicates.
if [ -f ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.rmdup.bam ]
then
    echo "Removed duplicates file: ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.rmdup.bam already exists. Moving on..."
else
    echo "remove duplicates" `date`
    (time samtools rmdup ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.mapped.bam ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.rmdup.bam | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/$LOG_ERR
    (time samtools index ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.rmdup.bam | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/$LOG_ERR
fi

# Make Fragment Bed
if [ -f ${pipeline}/${OUTDIR}/${OUTDIR}_fragment.bed ]
then
    echo "Fragmend Bed file: ${pipeline}/${OUTDIR}/${OUTDIR}_fragment.bed already exists. Moving on..."
else
    echo "Making fragment bed file"
    (time perl $scripts/getFragmentBed.pl ${pipeline}/${OUTDIR}/${OUTDIR}.sorted.rmdup.bam ${pipeline}/${OUTDIR}/${OUTDIR}_fragment.bed | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/$LOG_ERR
fi

if [ -f ${pipeline}/${OUTDIR}/${OUTDIR}.bedGraph ]
then
    echo "bedGraph file: ${pipeline}/${OUTDIR}/${OUTDIR}.bedGraph already exists. Moving on..."
else
    # Make UCSC Tracks
    # Bd21 is ostensibly the database I'm accessing but it just needs the command for sanity checking (nothing?)
    echo "make bedgraph" `date`
    (time sort -k1,1 -k2,2n ${pipeline}/${OUTDIR}/${OUTDIR}_fragment.bed | bedItemOverlapCount Bd21 chromSize=$scripts/Bd21.chrom.sizes stdin | sort -k1,1 -k2,2n > ${pipeline}/${OUTDIR}/${OUTDIR}.bedGraph | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/${LOG_ERR}
fi

# Make bigWig file.
if [ -f ${pipeline}/${OUTDIR}/${OUTDIR}.bw ]
then
    echo "bigWig file: ${pipeline}/${OUTDIR}/${OUTDIR}.bw already exists. Moving on..."
else
    (time bedGraphToBigWig ${pipeline}/${OUTDIR}/${OUTDIR}.bedGraph $scripts/Bd21.chrom.sizes ${pipeline}/${OUTDIR}/${OUTDIR}.bw | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/$LOG_ERR
fi


# Before calling peaks (macs2) pre-shift the reads.
if [ -f ${pipeline}/${OUTDIR}/${OUTDIR}_fragment_preShift75.bed ]
then
    echo "Shifted 75 fragment bed file: ${pipeline}/${OUTDIR}/${OUTDIR}_fragment_preShift75.bed already exists. Moving on..."
else
    (time perl $scripts/preShift.pl ${pipeline}/${OUTDIR}/${OUTDIR}_fragment.bed 75 ${pipeline}/${OUTDIR}/${OUTDIR}_fragment_preShift75.bed | tee -a ${pipeline}/${OUTDIR}/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a ${pipeline}/${OUTDIR}/$LOG_ERR
fi
