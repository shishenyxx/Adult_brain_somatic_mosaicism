##filtering for the raw AF data calculated from snMPAS 

library(ggplot2)

raw<-read.delim(file="single_cell_with_categories.tsv",header=T)
raw2<-raw[raw$CATEGORY!="MOSAIC"&((raw$REF_COUNT+raw$ALT_COUNT)>30)&raw$ID!="JGG-Cells",]
raw5<-raw[raw2$REF_COUNT+raw2$ALT_COUNT>30,]
raw6<-raw5[raw$CATEGORY=="HET",]
raw7<-raw5[raw$CATEGORY=="REF_HOMO",]

bks=c(0,0.0005,0.01,0.1,0.2,0.3,1)
pdf(file="2020_07_03_AF_Distribution_lower_homo_single_cell.pdf",width=4,hight=8)
ggplot(raw7,aes(x=sqrt(LOWER_CI),fill=CATEGORY))+
	geom_histogram(bins=201,fill="grey")+
	scale_x_continuous(breaks=sqrt(bks),labels=bks)+
	geom_vline(xintercept=sqrt(5.425525e-04),linetype=2)+
#	geom_vline(xintercept=sqrt(0.5),linetype=2)+
#	geom_vline(xintercept=sqrt(0.4),linetype=2,col="pink")+
#	geom_vline(xintercept=sqrt(0.6),linetype=2,col="pink")+
#	geom_vline(xintercept=sqrt(0.35),linetype=2,col="blue")+
#	geom_vline(xintercept=sqrt(0.65),linetype=2,col="blue")+
	labs(title="Ref homo AF distribution (sqr-t)")+
#	xlim(0,1)+
	ylim(0,2000)+
	theme_classic()
dev.off()

bks2=c(0,0.1,0.4,0.5,0.6,1)
pdf(file="2020_07_03_AF_Distribution_upper_het_single_cell.pdf",width=4,height=8)
ggplot(raw6,aes(x=sqrt(UPPER_CI),fill=CATEGORY))+
	geom_histogram(bins=201)+
	scale_x_continuous(breaks=sqrt(bks2),labels=bks2)+
#	geom_vline(xintercept=sqrt(0.0014),linetype=2)+
	geom_vline(xintercept=sqrt(0.5),linetype=2)+
	geom_vline(xintercept=sqrt(0.4),linetype=2,col="blue")+
	geom_vline(xintercept=sqrt(0.6),linetype=2,col="blue")+
	geom_vline(xintercept=sqrt(0.35),linetype=2,col="pink")+
	geom_vline(xintercept=sqrt(0.65),linetype=2,col="pink")+
	labs(title="Het and homo AF distribution (sqr-t), after removing the 7 variants")+
#	xlim(0,1)+
	ylim(0,2000)+
	theme_classic()
dev.off()
