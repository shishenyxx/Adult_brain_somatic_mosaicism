###codes for permutation
##L_T
#in shell
for i in {1..10000};do shuf seed_L_T.txt|sed ":a;N;s/\n//g;ta"|awk -v OFS="\t" '{print substr($1,1,2)"g"substr($1,3,9)}'|awk -v OFS="\t" '{print $1,substr($1,6,2),substr($1,9,2),substr($1,11,2)}'>>permutation_L_T.txt;done

paste <(paste <(my_join.pl -F 1 -f 1 -a <(cat permutation_L_T.txt|cut -f2) -b distance.txt |cut -f3) <(my_join.pl -F 1 -f 1 -a <(cat permutation_L_T.txt|cut -f3) -b distance.txt |cut -f3)) <(my_join.pl -F 1 -f 1 -a <(cat permutation_L_T.txt|cut -f4) -b distance.txt |cut -f3)|awk '{print $1+$2+$3}'>numbers_permutation_L_T.txt

#in R
raw<-read.table(file="numbers_permutation_L_T.txt",header=F)
ggplot(raw,aes(V1))+
	geom_histogram()+
	geom_vline(xintercept=3,col="red")+
	theme_classic()

##R_PF
#in shell
for i in {1..10000};do shuf seed_R_PF.txt|sed ":a;N;s/\n//g;ta"|awk -v OFS="\t" '{print $1,substr($1,1,2),substr($1,4,2),substr($1,6,2),substr($1,10,2)}'>>permutation_R_PF.txt;done

paste <(paste <(my_join.pl -F 1 -f 1 -a <(cat permutation_R_PF.txt|cut -f2) -b distance.txt |cut -f3) <(my_join.pl -F 1 -f 1 -a <(cat permutation_R_PF.txt|cut -f3) -b distance.txt |cut -f3)) <(paste <(my_join.pl -F 1 -f 1 -a <(cat permutation_R_PF.txt|cut -f4) -b distance.txt |cut -f3) <(my_join.pl -F 1 -f 1 -a <(cat permutation_R_PF.txt|cut -f5) -b distance.txt |cut -f3))|awk '{print $1+$2+$3+$4}'>numbers_permutation_R_PF.txt

raw<-read.table(file="numbers_permutation_R_PF.txt",header=F)

#in R
ggplot(raw,aes(V1))+
	geom_histogram()+
	geom_vline(xintercept=4,col="red")+
	theme_classic()

##L_PF
#in shell
for i in {1..10000};do shuf seed_L_PF.txt|sed ":a;N;s/\n//g;ta">>permutation_L_PF.txt;done
my_join.pl -F 1 -f 1 -a <(cat permutation_L_PF.txt|sort|uniq|awk '{print substr($1,4,2)}') -b distance.txt>numbers_permutation_L_PF.txt

For PF_L: p = 48/120 = 0.4

#in R
raw<-read.table(file="numbers_permutation_L_PF.txt",header=F)
ggplot(raw,aes(V3))+
	geom_histogram()+
	geom_vline(xintercept=1,col="red")+
	xlim(0,30)+
	theme_classic()
	
##all 3 lobes	
#in shell
paste numbers_permutation_L_T.txt numbers_permutation_R_PF.txt <(my_join.pl -F 1 -f 1 -a <(cat permutation_L_PF.txt|awk '{print substr($1,4,2)}') -b distance.txt|cut -f3)|awk '{print $1+$2+$3+$4}'>numbers_permutation_all_three_lobes.txt

#in R
raw<-read.table(file="numbers_permutation_all_three_lobes.txt",header=F)
ggplot(raw,aes(V1))+
	geom_histogram()+
	geom_density()+
	geom_vline(xintercept=8,col="red")+
	xlim(0,30)+
	theme_classic()
