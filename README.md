# liftover_helper

## Requirements
Python3

## Packages:
pyvcf

## Setup Project

### Download repo

git clone https://github.com/color/liftover_helper.git

cd liftover_helper

python3 -m venv ./venv

### run make
make

### run tests to make sure it works
python -m unittest tests/test_liftover.py

### Running from commandline

python scripts/liftover.py convert_grch38_ref_mismatch_sites_to_grch37 input_file output_path

python scripts/liftover.py convert_hg19_vcf_to_grch37_vcf input_file output_file

## Description and Notes
For update_grch38_ref_to_grch37_for_record_if_needed():

1. If the genotype is malformed or if there are more than 2 alleles, the record is not handled and written as is in the output file. Picard's LiftoverVcf will put them in reject file.

2. TODO: For records that are updated, update or drop the following fields

  - AD
  - AF


