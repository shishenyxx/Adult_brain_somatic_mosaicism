#Generate random samples from gnomAD

for i in {1..10000};do shuf -n 4000 ../gnomAD_strelka_intersect_no_chr_no_any_repeats.bed >> ${i}.random;done

for i in {1..10000}; do cat random/${i}.random|sed 's/:/\t/g' |awk -v OFS="\t" '{print $1, $2-1,$3+1,$4,$5}'>bed/${i}.bed;done
