#!/usr/bin/bash

set -eu

input_bed=$1 # CHROM, POS1, POS2 of variants
output_pdf=$2

b37_10mb_bin_bed=data/b37.10mb.5mb.count.bed
bin_count_bed=$(basename ${input_bed}).count.bed
DRAW_PLOT=src/draw_genomic_dist_plot.R
fai_path=data/human_g1k_v37.fasta.fai

bedtools intersect -a ${b37_10mb_bin_bed} -b ${input_bed} -c -loj > ${bin_count_bed}
${DRAW_PLOT} ${bin_count_bed} ${output_pdf} ${fai_path}
