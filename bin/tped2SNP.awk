
{printf "%s", $2; for(i=1;i<(NF/2-1);i++){if($(2*i+3)=="0" || $(2*i+4)=="0") printf "\tNA"; else printf "\t%s", $(2*i+3)+$(2*i+4)-2}; printf "\n"}
