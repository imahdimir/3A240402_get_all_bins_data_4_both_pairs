#!/usr/bin/bash

# cwd == project root

# filters 1000 snps for each info score in .3:.01:1 range
#   from all over the genome, 70 * 1K = 70K snps
python py/filter_snps.py

# Obtains the bgen files for the filtered snps
bash sh/filter_bgen.sh
