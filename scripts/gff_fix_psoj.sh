#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --job-name=gff_fix
#SBATCH --output=slurm-gff_fix-%j.out

## Bash strict mode
set -euo pipefail

## Hardcoded variables
gff_org=$1
gff_ed=$2

## Test
[[ ! "$#" -eq 2 ]] && echo "## ERROR: Please pass 2 arguments to the script; you passed $#" && exit 1
[[ ! -f "$gff_org" ]] && echo "## ERROR: Input file ($gff_org) does not exist" && exit 1

## Output dir
mkdir -p "$(dirname "$gff_ed")"

## Report
echo "## Running gff_fix_psoj.sh for $gff_org..."

# MAIN -------------------------------------------------------------------------
# add gene_id reference:
sed -E '/\tgene\t/s/\tID=([^;]+);/\tID=\1;gene_id=\1;/' "$gff_org" |
    sed -E '/\texon\t/s/locus_tag=([^;]+);/locus_tag=\1;gene_id=gene-\1;/' \
    > "$gff_ed"

## Wrap up
echo "## Listing output file:"
ls -lh "$gff_ed"
echo
