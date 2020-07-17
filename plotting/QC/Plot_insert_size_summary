#Codes for plotting insert size

library(ggplot2)

raw<-read.delim(file="summary_insert_size.txt",header=T)
names(raw)
##[1] "name"    "group"   "size"    "count"   "color"   "L_R"     "Lrg_Sml"

ggplot(raw,aes(x=size,y=count,col=color,group=group))+
	geom_line()+
	#scale_colour_manual(values = c("#ff5b00","brown","#107010"))+
	#facet_wrap(~cohort,ncol=2)+
	xlim(0,1500)+
	theme_classic()
