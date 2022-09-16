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
echo "## Running gff_fix_psan.sh for $gff_org..."

# MAIN -------------------------------------------------------------------------
# add gene_id reference:
sed -E '/\tgene\t/s/\tID=g([0-9]+);$/\tID=g\1;gene_id=g\1;/' "$gff_org" |
    sed -E '/\texon\t/s/Parent=([^;]+);/Parent=\1;gene_id=\1;/' |
    sed -E '/\texon\t/s/gene_id=g([0-9]+)\.t[0-9]+/gene_id=g\1/' \
    > "$gff_ed"

## Wrap up
echo "## Listing output file:"
ls -lh "$gff_ed"
echo
