#Plot the cumulative depth of coverage based on bedtools coverage

library(ggplot2)

raw<-read.delim(file="genome.hist.cumulative.txt",header=T)
names(raw)
##[1] "contig" "cov"    "n"      "leng"   "prop"   "sample" "size"   "side"  
##[9] "Organs"

ggplot(raw,aes(x=cov,y=cumulative_prop,col=Organs,group=sample))+
	geom_line()+
	#scale_colour_manual(values = c("#ff5b00","brown","#107010"))+
	#facet_wrap(~cohort,ncol=2)+
	ylim(0,1)+
	theme_classic()
