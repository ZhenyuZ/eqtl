# multiple bash scripts to run shapeit for pre-phasing

# split genotypes by chromosomes
for i in `seq 1 22`
do
  plink --noweb --bfile euro.strandfixed --chr $i --make-bed --out chr$i & 
done

# shapeit check strand
for i in `seq 1 22`
do 
    shapeit -check -B chr$i --input-ref /home/ubuntu/SCRATCH/1000GP_Phase3/1000GP_Phase3_chr$i.hap.gz /home/ubuntu/SCRATCH/1000GP_Phase3/1000GP_Phase3_chr$i.legend.gz /home/ubuntu/SCRATCH/1000GP_Phase3/1000GP_Phase3.sample --output-log chr$i & 
done

# extract slip strand and missing snp information
for i in `seq 1 22`
do
    cat chr$i.snp.strand | grep "Strand" | cut -f4 | uniq > chr$i.strand_flip
    cat chr$i.snp.strand | grep "Missing" | cut -f4 | uniq > chr$i.missing
done
 
# flip snp and exclude missing ones
for i in `seq 1 22`
do
    plink --noweb --bfile chr$i --flip chr$i.strand_flip --exclude chr$i.missing --recode --out chr$i.flipped --make-bed &
done

# shapeit check strand again
for i in `seq 1 22`
do 
    shapeit -check -B chr$i.flipped --input-ref /home/ubuntu/SCRATCH/1000GP_Phase3/1000GP_Phase3_chr$i.hap.gz /home/ubuntu/SCRATCH/1000GP_Phase3/1000GP_Phase3_chr$i.legend.gz /home/ubuntu/SCRATCH/1000GP_Phase3/1000GP_Phase3.sample --output-log chr$i.again & 
done

# extract different snps
for i in `seq 1 22`
do
    cat chr$i.again.snp.strand | grep "rs" - | cut -f4 | uniq > chr$i.different
done

# exclude different snps
for i in `seq 1 22`
do
    plink --noweb --bfile chr$i.flipped --exclude chr$i.different --recode --out chr$i.flipped.again --make-bed &
done



# shapeit check strand again
for i in `seq 1 22`
do 
    shapeit -check -B ../chr$i.flipped.again --input-ref /home/ubuntu/SCRATCH/1000GP_Phase3/1000GP_Phase3_chr$i.hap.gz /home/ubuntu/SCRATCH/1000GP_Phase3/1000GP_Phase3_chr$i.legend.gz /home/ubuntu/SCRATCH/1000GP_Phase3/1000GP_Phase3.sample --output-log chr$i.again & 
done

# collect missing and different snps
for i in `seq 1 22`
do
    cat chr$i.missing >> all.missing.snps
    cat chr$i.different >> all.different.snps
done

# recover genotypes
cat all.missing.snps all.different.snps > all.recover.snps
plink --noweb --bfile euro.strandfixed --extract all.recover.snps --make-bed --out euro.recovered

# main bash script to run shapeit
#!/bin/bash
## shapeit prephasing on clusterfs main bash script, need to replace ##par1## by 1:22

mkdir -p ~/SCRATCH/shapeit
cd ~/SCRATCH/shapeit
mkdir -p ref
cd ~/SCRATCH/shapeit/ref
s3cmd get s3://bioinformatics_scratch/zhenyu/shapeit/ref/1000GP_Phase3.sample
s3cmd get s3://bioinformatics_scratch/zhenyu/shapeit/ref/1000GP_Phase3_chr##par1##.legend.gz
s3cmd get s3://bioinformatics_scratch/zhenyu/shapeit/ref/1000GP_Phase3_chr##par1##.hap.gz
s3cmd get s3://bioinformatics_scratch/zhenyu/shapeit/ref/genetic_map_chr##par1##_combined_b37.txt
cd ~/SCRATCH/shapeit
mkdir -p output
cd ~/SCRATCH/shapeit/output
export http_proxy=http://cloud-proxy:3128; export https_proxy=http://cloud-proxy:3128;
sudo -E apt-get install libboost-iostreams1.55.0
sudo -E apt-get install libboost-program-options1.55.0
shapeit \
        --input-bed    ~/shapeit/data/chr##par1##.flipped.again --thread 30 \
        --input-map    ~/SCRATCH/shapeit/ref/genetic_map_chr##par1##_combined_b37.txt \
        --output-max   ~/SCRATCH/shapeit/output/##par1##.phased \
        --output-graph ~/SCRATCH/shapeit/output/##par1##.hgraph \
        --output-log   ~/SCRATCH/shapeit/output/##par1##.log \
        --input-ref    ~/SCRATCH/shapeit/ref/1000GP_Phase3_chr##par1##.hap.gz \
                       ~/SCRATCH/shapeit/ref/1000GP_Phase3_chr##par1##.legend.gz \
                       ~/SCRATCH/shapeit/ref/1000GP_Phase3.sample

s3cmd put ~/SCRATCH/shapeit/output/* s3://bioinformatics_scratch/zhenyu/shapeit/chr##par1##/


# generate scripts for each chromosome
for i in `seq 1 22`; do sed "s/##par1##/$i/g" getHaplotypes.sh > ./sh/get.chr$i.sh ; done


