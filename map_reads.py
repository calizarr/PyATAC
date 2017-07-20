#!/usr/bin/python
from __future__ import print_function
import argparse
import subprocess
import os
import sys

def options():
    parser = argparse.ArgumentParser(description="Map ATAC-seq reads with bowtie")
    parser.add_argument("nThreads", help="Number of threads to use for bowtie")
    parser.add_argument("input_file", help="Read file to be mapped with bowtie")
    parser.add_argument("out_file_dir", help="Directory for output file (directory for the sam files)")
    parser.add_argument("basepre", help="Base prefix for output file (Everything before .sam)")
    parser.add_argument("-i", "--indices",
                        default="/shares/tmockler_share/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/assembly/indices/Bdistachyon_314_v3.0.hardmasked",
                        help="Directory of bowtie indices (default: Brachy)")
    args = parser.parse_args()
    return args


def main():
    args = options()
    # Getting the absolute path to file if relative path given.
    args.input_file = os.path.abspath(args.input_file)
    # Getting the directory name
    dirname = os.path.dirname(args.input_file)
    # Getting just the filename
    output_file = os.path.join(dirname, args.out_file_dir, args.basepre)
    print(output_file)
    indices = args.indices
    # Map reads to reference, convert to bam and sort via pipes.
    subprocess.call('echo "bowtie alignment" `date`')
    cmd = "bowtie --chunkmbs 256 -p {0} -S -m 1 -X 2000 -t {1} {2} | samtools view - -bS | samtools sort - -o {3}.sorted.bam -O bam -T {2}.pre" \
          .format(args.nThreads, indices, args.input_file, output_file)
    print("Running cmd: ")
    print(cmd)
    subprocess.call(cmd, shell=True)
    # Indexing bam file.
    cmd = "samtools index {0}.sorted.bam".format(output_file)
    print("Running cmd: ")
    print(cmd)
    subprocess.call(cmd, shell=True)
    # Only mapped reads.
    cmd = "samtools view -b -F 4 -o {0}.sorted.mapped.bam {0}.sorted.bam".format(output_file)
    print("Running cmd: ")
    print(cmd)
    subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
