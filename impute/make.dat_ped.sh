#!/bin/bash
# prepare ped file for imputation
plinkAll="/glusterfs/users/ZHANGZ18/euro/luad3/geno3/luad.clean"

for i in `seq 1 22` 
do
	dir="chr"$i
	mkdir $dir
	chr=$(echo $dir"/chr")
	bimFile=$(echo $dir"/chr.bim")
	pedFile=$(echo $dir"/chr.ped")
	datFile=$(echo $dir"/chr.dat")
	tempFile=$(echo $dir"/temp")
	plink --noweb --bfile $plinkAll --chr $i --make-bed --out $chr
	awk '{print "M\t"$2}' $bimFile > $datFile
	plink --noweb --bfile $chr	--recode --out $chr
	cut -d" " -f1-4 $pedFile > $(echo $tempFile"1")
	cut -d" " -f5 $pedFile > $(echo $tempFile"2")
	sed 's/1/M/g' $(echo $tempFile"2") | sed 's/2/F/g' > $(echo $tempFile"3")
	cut -d" " -f7- $pedFile > $(echo $tempFile"4")
	paste -d" " $(echo $tempFile"1") $(echo $tempFile"3") $(echo $tempFile"4") > $tempFile
	mv $tempFile $pedFile
	rm $(echo $tempFile"*")
done


