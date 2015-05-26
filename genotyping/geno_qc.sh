#!/bin/bash
# Zhenyu Zhang
#
# Input: cell_file_list.txt in the current directory
#			first line: cel_files
#			other lines: full path to individual CEL file per line
#
# Output: results.txt in the current directory

apt-geno-qc \
  --cdf-file /glusterfs/users/ZHANGZ18/meta/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
  --qca-file /glusterfs/users/ZHANGZ18/meta/GenomeWideSNP_6/GenomeWideSNP_6.r2.qca \
  --qcc-file /glusterfs/users/ZHANGZ18/meta/GenomeWideSNP_6/GenomeWideSNP_6.r2.qcc \
  --chrX-probes /glusterfs/users/ZHANGZ18/meta/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
  --chrY-probes /glusterfs/users/ZHANGZ18/meta/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
  --cel-files ./cel_file_list.txt \
  --out-file ./results.txt