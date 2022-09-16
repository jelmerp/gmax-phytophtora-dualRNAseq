#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --job-name=gff_fix
#SBATCH --output=slurm-gff_fix-%j.out

## Bash strict mode
set -euo pipefail

## Command-line args
gff_org=$1
gff_ed=$2

## Test
[[ ! "$#" -eq 2 ]] && echo "## ERROR: Please pass 2 arguments to the script; you passed $#" && exit 1
[[ ! -f "$gff_org" ]] && echo "## ERROR: Input file ($gff_org) does not exist" && exit 1

## Make output dir
mkdir -p "$(dirname "$gff_ed")"

## Report
echo "## Running gff_fix_gmax.sh for $gff_org..."

# MAIN -------------------------------------------------------------------------
# First, Need to add dummy 'gene_name' and 'gene_type' attributes
# for compatibility with workflow's `extract_annotations_from_gff.py` script
sed -E '/^\w.*/s/$/;gene_name=x;gene_type=x/' "$gff_org" |
    # Need to change the 'Name' ID to 'gene_id' for genes,
    # because this is erroneously hardcoded in the workflow's `extract_annotations_from_gff.py` script
    sed -E '/\tgene\t/s/;Name=/;gene_id=/' |
    # Finally, add 'gene_id' references to exons
    sed -E '/\texon\t/s/ID=([^;]+);/ID=\1;gene_id=\1;/' |
    sed -E 's/gene_id=Glyma\.([[:alnum:]]+)\.[0-9]+\.Wm82.a4.v1[^;]*;/gene_id=Glyma.\1.Wm82.a4.v1;/' \
    > "$gff_ed"

## Wrap up
echo "## Listing the output file:"
ls -lh "$gff_ed"
echo
