##Supplement data for AF measured by snMPAS
library(ggplot2)

raw<-read.csv(file="20200618_4dbsm_mos_259.csv",header=T)
raw2<-raw[raw$TISSUE_CELLS=="cells",]
bks=c(0,0.005,0.1,0.20,0.3,0.4,0.5,0.6,1)

pdf(file="4DBSM_Plot_single_cells.pdf",width=16,height=9)
for(i in 1:259){
temp<-raw2[raw2$CHR_POS_REF_ALT==names(table(raw2$CHR_POS_REF_ALT))[i],]
plot <-
ggplot(temp,aes(x=ID,y=sqrt(MAF),ymin=sqrt(LOWER_CI),ymax=sqrt(UPPER_CI),fill=Cell.Type,col=Cell.Type))+
	geom_point(alpha=0.5)+
	geom_pointrange()+
	geom_hline(yintercept=sqrt(0.005),linetype=2,col="black")+
	geom_hline(yintercept=sqrt(0.6),linetype=2,col="grey")+
	geom_hline(yintercept=sqrt(0.5),linetype=2,col="grey")+
	geom_hline(yintercept=sqrt(0.4),linetype=2,col="grey")+
	geom_hline(yintercept=sqrt(1),linetype=2,col="grey")+
	geom_hline(yintercept=sqrt(0),linetype=2,col="grey")+
	scale_y_continuous(breaks=sqrt(bks),labels=bks)+
	scale_color_manual(values = c("pink","#2A3280","#A73326","#666666","#CCCCCC"))+
	theme_classic()+
	labs(title=names(table(raw2$CHR_POS_REF_ALT))[i],x="")+
	theme(axis.line=element_blank(),axis.text.x=element_text(face="bold",size=5,angle=90,vjust = 0.5))+
	geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1),color="black")
print(plot)
}
dev.off()
