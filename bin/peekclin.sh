#! /bin/bash
cat $1 | sed 's/\[Not_Available\]/NA/g' > temp.clin
head -n+1 temp.clin > temp.head
tail -n+2 temp.clin > temp.body
NF=$(awk '{print NF}' temp.head)
for i in $(seq 2 $NF); do echo; echo "Column #"$i; cut -f$i temp.head; cut -f$i temp.body | sort | uniq -c ; done
echo
echo "Total samples:"
cat temp.body | wc -l
echo "Samples with missing values:"
grep -n $"\tNA" temp.body | wc -l

