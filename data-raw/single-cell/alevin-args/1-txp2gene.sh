# install bioawk
sudo apt install bison
git clone git://github.com/lh3/bioawk.git && cd bioawk && make && mv awk bioawk && sudo cp bioawk /usr/local/bin/

cd ~/Documents/Batcave/zaklab/drugseqr/data-raw/single-cell/txp2hgnc

# transcript to ensembl gene tsv
wget -nv ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
bioawk -c gff '$feature=="transcript" {print $attribute}' <(gunzip -c gencode.v29.annotation.gtf.gz) | awk -F ' ' '{print substr($4,2,length($4)-3) "\t" substr($2,2,length($2)-3)}' - > txp2gene.tsv
