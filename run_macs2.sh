#!/usr/bin/bash
source /shares/tmockler_share/clizarraga/usr/virtualenvs/MACS2/bin/activate

# A POSIX Variable. Reset in case getopts has been used in shell.
OPTIND=1

while getopts "h?t:c:n:o:mf:" opt; do
    case "$opt" in
        h|\?)
            macs2 callpeak -h
            echo "This script needs only: -t -c -n -o (--outdir) -m (--nomodel) -f (--extsize)"
            exit 0
            ;;
        t)
            treatment_file=$OPTARG
            echo "Treatment file is $treatment_file"
            ;;
        c)
            control_file=$OPTARG
            echo "Control file is $control_file"
            ;;
        n)
            prefix_name=$OPTARG
            echo "Prefix for filenames is $prefix_name"
            ;;
        o)
            outdir=$OPTARG
            echo "Output directory is $outdir"
            ;;
        m)
            nomodel="--nomodel"
            echo "No model selected."
            ;;
        f)
            frag=$OPTARG
            echo "Fragment size is $frag"
    esac
done

# Broad model calling
macs2 callpeak -t $treatment_file -c $control_file -g 2.584e8 -B --broad --broad-cutoff 0.05 -n $prefix_name --outdir $outdir -f BAM $nomodel --extsize $frag
