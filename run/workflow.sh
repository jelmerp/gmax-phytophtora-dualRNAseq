# SETUP ------------------------------------------------------------------------
## Dirs and files
fqdir=data/fastq/220802_Meulia_GSL-FCC-2946/     # Dir with raw FASTQ files

fasta_host=data/ref/gmax_ncbi/GCF_000004515.6_Glycine_max_v4.0_genomic.fa       # Host reference genome FASTA file
gff_host=data/ref/gmax_ncbi/GCF_000004515.6_Glycine_max_v4.0_genomic.gff        # Host reference genome GFF file
fasta_psoj=data/ref/psojae/GCF_000149755.1_P.sojae_V3.0_genomic.fa              # Pathogen #1 reference genome FASTA file
gff_psoj=data/ref/psojae/GCF_000149755.1_P.sojae_V3.0_genomic.gff               # Pathogen #1 reference genome GFF file
fasta_psan=data/ref/psan65/Genome_P65_45X_final.fa                                     # Pathogen #2 reference genome FASTA file
gff_psan=data/ref/psan65/BRAKER2_P65.nomask_augustus.hints.gff3                        # Pathogen #2 reference genome GFF file

## Settings
minq=5                          # Min. qual for TrimGalore trimming
minlen=36                       # Min. read length for TrimGalore trimming
kraken_conf=0.5                 # Kraken confidence required for assignment (Kraken's `--confidence` flag)
gff_feature_type="exon"         # The relevant name in the 3rd column of the GFF/GTF (e.g. 'gene' or 'exon')
gff_feature_id="gene_id"        # The key used in the last column of the GFF/GTF for identifying the gene-level parent (e.g. 'Name' or 'gene_id') 


# PREPROCESSING ----------------------------------------------------------------
## Combine FASTQ files for the two lanes
sbatch mcic-scripts/utils/fqconcat.sh -i "$fqdir" -o data/fastq/concat

## QC with FastQC -- separately for non-concat and lane-concat files
shopt -s globstar
for fq in "$fqdir"/**/*fastq.gz; do
    sbatch mcic-scripts/qc/fastqc.sh -i "$fq" -o results/fastqc/raw_nonconcat
done
sbatch mcic-scripts/qc/multiqc.sh -i results/fastqc/raw_nonconcat -o results/multiqc/fastqc_raw_nonconcat

for fq in data/fastq/concat/*fastq.gz; do
    sbatch mcic-scripts/qc/fastqc.sh -i "$fq" -o results/fastqc/raw_concat
done
sbatch mcic-scripts/qc/multiqc.sh -i results/fastqc/raw_concat -o results/multiqc/fastqc_raw_concat

## Remove adapters & very low qual bases with TrimGalore
for R1 in data/fastq/lane_concat/*_R1_*fastq.gz; do
    sbatch mcic-scripts/trim/trimgalore.sh \
        -i "$R1" -o results/trim -O results/fastqc/trim -q "$minq" -l "$minlen"
done
grep "Reads with adapters" slurm-trimgalore*    # Check how many reads contained adapters
grep "Total written" slurm-trimgalore*          # Check the % of bp that was retained

## Remove rRNA reads with SortMeRNA
repo_dir=results/sortmerna/sortmerna_repo && mkdir -p "$repo_dir"
git clone https://github.com/biocore/sortmerna "$repo_dir"
for R1 in results/trim/*_R1_*fastq.gz; do
    sbatch mcic-scripts/rnaseq/sortmerna.sh -i "$R1" -o results/sortmerna -r "$repo_dir"
done
grep "Number of reads mapped" slurm-sortmerna*   # Check how many reads were mapped to the rRNA database

## Check for contamination with Kraken2 
kraken_db=/fs/project/PAS0471/jelmer/refdata/kraken/std
for R1 in results/sortmerna/unmapped/*_R1_*fastq.gz; do 
    sbatch mcic-scripts/metagenomics/kraken-run.sh \
        -i "$R1" -o results/kraken -d "$kraken_db" -W -w -c "$kraken_conf"
done

## Visualize Kraken output with Krona
for kraken_out in results/kraken/*_main.txt; do
    sampleID=$(basename "$kraken_out" _main.txt)
    sbatch mcic-scripts/metagenomics/krona.sh -i "$kraken_out" -o results/krona/"$sampleID".html
done

## Run FastQC again on final ("filtered") FASTQ files
for fq in results/kraken/unclassified/*fastq.gz; do
    sbatch mcic-scripts/qc/fastqc.sh -i "$fq" -o results/fastqc/post_kraken
done
sbatch mcic-scripts/qc/multiqc.sh -i results/fastqc/post_kraken -o results/multiqc/fastqc_post_kraken


# MAPPING & COUNTING WITH TWO REF. GENOMES -------------------------------------
## "Fix" GFF files for use with the workflow
mkdir -p results/gff_fix
gff_host_ed=results/gff_fix/$(basename "$gff_host")
gff_psoj_ed=results/gff_fix/$(basename "$gff_psoj")
gff_psan_ed=results/gff_fix/$(basename "$gff_psan")

bash scripts/gff_fix_gmax_ncbi.sh "$gff_host" "$gff_host_ed"
bash scripts/gff_fix_psoj.sh "$gff_psoj" "$gff_psoj_ed"
bash scripts/gff_fix_psan.sh "$gff_psan" "$gff_psan_ed"

## Run Nextflow workflow with G.max + P.sojae
fqd=results/fqsub  #fqd=results/kraken/unclassified
sbatch mcic-scripts/rnaseq/nf_dual.sh -i "$fqd" \
    -o results/nfdual_psoj -c workflow/conf/nf_dual.config \
    -f "$fasta_host" -g "$gff_host_ed" -F "$fasta_psoj" -G "$gff_psoj_ed"

## Run Nextflow workflow with G.max + P.santomeana
fqd=results/fqsub  #fqd=results/kraken/unclassified
sbatch mcic-scripts/rnaseq/nf_dual.sh -i "$fqd" \
    -o results/nfdual_psan -c workflow/conf/nf_dual.config \
    -f "$fasta_host" -g "$gff_host_ed" -F "$fasta_psan" -G "$gff_psan_ed"


# MAPPING & COUNTING WITH ONE REF. GENOME --------------------------------------
## Index the reference genome with STAR
sbatch mcic-scripts/rnaseq/star_index.sh -i "$fasta_host" -a "$gff_host" -o results/map/star_idx

## Align the reads to the genome with STAR
for R1 in results/kraken/unclassified/*_R1*fastq.gz; do
    sbatch mcic-scripts/rnaseq/star_align.sh \
        -i "$R1" -a "$gff_host" -r results/map/star_idx -o results/map
done
sbatch mcic-scripts/qc/multiqc.sh -i results/map/star_logs -o results/multiqc/star

## QC the mapping results with Qualimap
#i If the original annotation is in GFF format, first convert to GTF since Qualimap needs a GTF file
ref_gtf=${gff_host/.gff/.gtf}
sbatch mcic-scripts/convert/gff2gtf.sh -i "$gff_host" -o "$ref_gtf"

for bam in results/map/*bam; do
    sbatch mcic-scripts/rnaseq/qualimap.sh -i "$bam" -a "$ref_gtf" -o results/qualimap
done
sbatch mcic-scripts/qc/multiqc.sh -i results/qualimap -o results/multiqc/qualimap

## Produce a gene count table with Featurecounts
sbatch mcic-scripts/rnaseq/featurecounts.sh -i results/map -a "$ref_gtf" \
    -o results/fcounts/counts.txt -t "$gff_feature_type" -g "$gff_feature_id"
sbatch mcic-scripts/qc/multiqc.sh -i results/fcounts -o results/multiqc/featurecounts
