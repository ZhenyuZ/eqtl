#!/bin/bash
## script to run IMPUTE2 by chromosome chunks
## need to replace ##NAME##, ##CHR##, ##START##, ##END##

refbucket=s3://bioinformatics_scratch/zhenyu/shapeit/ref
reflegend=1000GP_Phase3_chr##CHR##.legend.gz
refhap=1000GP_Phase3_chr##CHR##.hap.gz
refmap=genetic_map_chr##CHR##_combined_b37.txt

hapbucket=s3://bioinformatics_scratch/zhenyu/shapeit/phase
hap=##CHR##.phased.haps

outdir=~/SCRATCH/impute2
rm -rf $outdir
mkdir -p $outdir
cd $outdir

outfile=$outdir/chunk.##NAME##.impute2
outbucket=s3://bioinformatics_scratch/zhenyu/impute2

s3cmd get $refbucket/$reflegend
s3cmd get $refbucket/$refhap
s3cmd get $refbucket/$refmap
s3cmd get $hapbucket/$hap

impute2 \
   -m $refmap \
   -known_haps_g $hap \
   -h $refhap \
   -l $reflegend \
   -Ne 20000 \
   -k_hap 1000 \
   -int ##START## ##END## \
   -o $outfile \
   -filt_rules_l 'EUR<0.01' \
   -seed 367946

s3cmd put chunk.##NAME##.* $outbucket/chunk.##NAME##/
