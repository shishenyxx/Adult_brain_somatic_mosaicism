##Plot the estimated starting population of variants shared by both hemispheres and the organs

library(ggplot2)

raw<-read.csv(file="left_right_violins.csv",header=T)

ggplot(raw,aes(x=Group,y=log10(Estimated_number)))+
	geom_violin()+
	geom_jitter()+
	geom_segment(aes(x = 0, y = 0, xend = 0, yend = 5),color="black")+
	theme_classic()+
	geom_hline(yintercept=log10(211),linetype=2,col="grey")+
	geom_hline(yintercept=log10(86),linetype=2,col="grey")+
	theme(axis.line=element_blank())
