#!/bin/bash
	plink --noweb --bfile euro --chr $i  --make-bed --out chr$i
	awk '{print "M\t"$2}' chr$i.bim > chr$i.dat
	plink --noweb --bfile chr$i  --recode --out chr$i
	cut -d" " -f1-4 chr$i.ped > temp1.chr$i
	cut -d" " -f5 chr$i.ped > temp2.chr$i
	sed 's/1/M/g' temp2.chr$i | sed 's/2/F/g' > temp3.chr$i
	cut -d" " -f7- chr$i.ped > temp4.chr$i
	paste -d" " temp1.chr$i temp3.chr$i temp4.chr$i > chr$i.ped

