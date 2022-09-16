## Subsample FASTQ files - and only for '24' samples to get a small dataset
fqdir_raw=data/fastq/concat
fqdir_sub=results/fqsub
bash mcic-scripts/utils/fqsub_dir.sh -i "$fqdir_raw" -o "$fqdir_sub" -x "24"

## Get Kraken std db - https://benlangmead.github.io/aws-indexes/k2
cd /fs/project/PAS0471/jelmer/refdata/kraken/
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220607.tar.gz
tar -xzvf k2_standard_20220607.tar.gz


# REFERENCE GENOMES ------------------------------------------------------------
## Copy psan genomes from Linda
mkdir -p data/ref/psan65/ data/ref/psan68/
cp ../../../linda/results/genome_assembly/Psan65_genome/Genome_P65_45X_final.fa data/ref/psan65/
cp ../../../linda/results/genome_assembly/Psan65_genome/BRAKER2_P65.nomask_augustus.hints.gff3 data/ref/psan65/

## Copy Gmax genome from Nghi's project
cp -r ../2021-10_nghi/data/ref/jgi/phytozome/Gmax/Wm82.a4.v1/ data/ref

## Download Gmax and P. sojae from NCBI
mkdir -p data/ref/psojae data/ref/gmax_ncbi
wget -P data/ref/psojae https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/755/GCF_000149755.1_P.sojae_V3.0/GCF_000149755.1_P.sojae_V3.0_genomic.fna.gz
wget -P data/ref/psojae https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/755/GCF_000149755.1_P.sojae_V3.0/GCF_000149755.1_P.sojae_V3.0_genomic.gff.gz
wget -P data/ref/gmax_ncbi https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz
wget -P data/ref/gmax_ncbi https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.gff.gz
for zfile in data/ref/psojae/*gz data/ref/gmax_ncbi/*gz; do gunzip "$zfile"; done
for fna in data/ref/psojae/*fna data/ref/gmax_ncbi/*fna; do mv -v "$fna" "${fna/.fna/.fa}"; done
