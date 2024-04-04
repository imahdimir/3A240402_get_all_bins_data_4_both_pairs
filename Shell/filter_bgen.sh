#!/usr/bin/bash

sample_fp="/disk/genetics2/ukb/orig/UKBv3/sample/ukb11425_imp_chr1_22_v3_s487395.sample"
fs_po_ids="/disk/genetics4/ws_dirs/mahdimir/db/C1-P-G/A1-Git-SF/get-all-bins-data-4-both-pairs-240402-SF/inp/fs_po_ids.txt"
plink_out="/disk/genetics4/ws_dirs/mahdimir/db/C1-P-G/A1-Git-SF/get-all-bins-data-4-both-pairs-240402-SF/med/plink_out"

bg_pat="/disk/genetics/data/ukb/private/v3/raw/imputed/ukb_imp_chr%d_v3.bgen"
snps_pat="/disk/genetics4/ws_dirs/mahdimir/db/C1-P-G/A1-Git-SF/get-all-bins-data-4-both-pairs-240402-SF/med/flt-snps/c%d.txt"

for i in {1..22}
do
  bg_fp=$(printf $bg_pat $i)
  snps_fp=$(printf $snps_pat $i)

  echo $bg_fp
  echo $snps_fp

  plink2 --bgen "$bg_fp" ref-first --export bgen-1.2 --sample "$sample_fp" --keep "$fs_po_ids" --extract "$snps_fp" --out "$plink_out/c$i"
done
