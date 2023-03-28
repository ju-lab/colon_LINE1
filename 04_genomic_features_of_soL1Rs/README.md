# Genomic Features of soL1Rs

This directory contains data and scripts to describe the genomic features of soL1R insertion sites

## Genomic distribution of soL1Rs

If you have a list of soL1Rs (or any variants) with their genomic coordinates, you can draw the distribution ofthe  variants across the genomic regions in 10 Mb sliding windows with 5-Mb-sized steps. 

### Prerequisite
R 3.6+ with R packages (tidyverse, grid)

### Command
```bash
./draw_genomic_distribution_plot.sh ${input_bed} ${output_pdf}
```

If you run this command with input_bed file with 3 columns for chromosome, position1, position2, the distribution of variants is saved in the file named output_pdf.
