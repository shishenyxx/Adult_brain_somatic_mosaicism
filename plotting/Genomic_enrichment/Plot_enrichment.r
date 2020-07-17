#Codes for plotting the genomic enrichment permutation result in R

library(ggplot2)

raw<-read.csv(file="4DBSM_plot.csv",header=T)
names(raw)
##[1] "GROUP"    "CATEGORY" "FRACTION" "Lower_CI" "Upper_CI" "Color"
ggplot(raw,aes(x=factor(GROUP,levels=c("All Variants","Brain and Organs","Brain only","Cortex only","Brain Lateralized","Brain Single Sample")),y=FRACTION,ymin=Lower_CI,ymax=Upper_CI))+
	geom_linerange(col=raw$Color,size=3)+
	geom_point(size=2,aes(fill=GROUP),color="black",pch=21)+
	scale_fill_manual(values = c("#808080","#CCCCCC","#86BAEC","#0B0E50","#D6F8F8","#1F4EC3"))+
	facet_wrap(~CATEGORY,nrow=1)+
	theme_classic()+
	theme(axis.text.x=element_blank(),strip.background = element_blank(),strip.text = element_text(angle=90))
