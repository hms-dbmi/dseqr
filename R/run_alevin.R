salmon_idx <- system.file('indices', 'homo_sapiens', package = 'drugseqr')

system2('salmon alevin -l ISR -1 10X_FID12518_Normal_3hg_S1_L002_R1_001.fastq.gz -2 10X_FID12518_Normal_3hg_S1_L002_R2_001.fastq.gz --chromium  -i /home/alex/R/x86_64-pc-linux-gnu-library/3.5/drugseqr/indices/homo_sapiens -p 8 -o alevin_output --tgMap /home/alex/Documents/Batcave/zaklab/drugseqr/inst/extdata/tgmap.tsv')
