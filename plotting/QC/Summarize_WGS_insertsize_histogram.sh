#Summarize the insertsize histogram created by GATK collect insertsize

paste <(ls ..|grep txt|sed 's/\.insert_size_metrics\.txt//g') <(ls ..|grep txt|sed 's/\.insert_size_metrics\.txt//g'|sed 's/-/\t/g'|cut -f2)|while read name group;do tail -n +12  ../${name}.insert_size_metrics.txt|awk -v OFS="\t" -v aaa=${name} -v bbb=${group} '{print aaa,bbb,$0}'|awk -v OFS="\t" '{if($3!="") print $0}'>>../summary_insert_size.txt;done
