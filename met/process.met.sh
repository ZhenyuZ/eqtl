#!/bin/bash
# Process TCGA methylation data downloaded from firehose to be simple matrix form of just beta values

while read -r line;
do
  awk -F"\t" 'NR!=2 {for(i=2;i<=NF;i=i+4) printf $i"\t"; printf "\n"}' $line.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > $line.raw.met.txt;
done < "../diseases.txt"

