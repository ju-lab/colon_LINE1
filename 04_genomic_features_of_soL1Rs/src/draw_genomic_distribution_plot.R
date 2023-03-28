#!/usr/bin/Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(grid))

args <- commandArgs(trailingOnly=TRUE)
### args1 = input file (chrom, start, end, length, count)
### args2 = output pdf file name
### args3 = path to fai file

# args[1] = "/home/users/changhyunnam/Projects/08_LINE1/08_Whole/99_scripts/bin/L1_normal_1198.10mb.5mb.count.bed"
# args[2] = "/home/users/changhyunnam/Projects/08_LINE1/03_Somatic/11_EMseq/03_VCF/04_Somatic_CpG/01_reference_genome/distribution/L1_normal_1198.10mb.5mb.count.bed.count.bed.pdf"
# args[3] = "/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta.fai"

data <- read_tsv(args[1],col_names=c("chrom","start","end","length","count"), col_types=cols(chrom="c")) %>% filter(length > 10**6)
outfile <- args[2]
fai  <- read_tsv(args[3],col_names=c("chrom","length","C","D","E"), col_types=cols(chrom="c",length="d")) %>% select(chrom,length)

chrom_sizes = structure(fai$length,names=fai$chrom)
chrom_sizes = chrom_sizes[str_detect(names(chrom_sizes), "^(chr)?[0-9XY]*$")] # only get autosomal + X Y chromosomes, ignore small contigs

chrom_sizes_acc = c()
chrom_sizes_acc[1] = 0
for(i in c(2:length(chrom_sizes))){
  chrom_sizes_acc[i] = sum(chrom_sizes[1:i - 1])
}
chrom_sizes_acc = structure(chrom_sizes_acc, names = names(chrom_sizes))

data$pos_genome = data$start + chrom_sizes_acc[as.character(data$chrom)]

xrange = max(data$pos_genome)
x = data$pos_genome / xrange
y = data$count

yunit = 10**round(log10(max(y)))
yrange = ceiling(max(y)/yunit)*yunit
y = y / yrange

pdf(outfile, width=21, height=7, useDingbats = F)

pushViewport(viewport(x=0, y=0, width=1, height=1, just=c('left', 'bottom')))
grid.rect(gp=gpar(lty="solid"))

pushViewport(viewport(x=0, y=0.95, height=0.05, just = c('left', 'bottom')))
grid.text(paste0("Genomic distribution"), x=0.05, y=0.5, just=c('left','bottom'))

popViewport(2)

pushViewport(viewport(x=0, y=0, height=0.95, just=c('left', 'bottom')))
pushViewport(viewport(x=0.5, y=0.5, height = 0.8, width = 0.9))
grid.rect(gp=gpar(lty="solid"))

xaxis_actual_pos = c(0)
for(i in c(1:length(chrom_sizes))){xaxis_actual_pos[i + 1] = sum(chrom_sizes[1:i])}
xaxis_plot_pos = xaxis_actual_pos / xrange
grid.xaxis(at = xaxis_plot_pos, label = FALSE)

yaxis_plot_pos = seq(0, 1, by = yunit/yrange)
grid.yaxis(at = yaxis_plot_pos, label = yaxis_plot_pos*yrange)
for(i in c(1:length(chrom_sizes))){
  grid.text(names(chrom_sizes)[i], x=(xaxis_plot_pos[i+1] + xaxis_plot_pos[i])/2, y=-0.05)
}
grid.text('Count per 10Mb', x=-0.04, rot=90)
grid.text('Chromosome', y=-0.1)

for(i in c(1:length(chrom_sizes))){
  if(i %% 2 == 1){bgcolor = 8}
  else{bgcolor = FALSE}
  grid.rect(x = xaxis_plot_pos[i], y = 0, just = c('left','bottom'), height = 1, width = xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE, alpha = 0.2))
}

line = grid.segments(x0 = x, x1 = x, y0 = 0, y1 = y)
grid.draw(line)

popViewport(2)

invisible(dev.off())
