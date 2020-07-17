##Quality control plots for reference homozygous and heterozygous control variants in MPAS analysis

library(ggplot2)

raw<-read.delim(file="2020_06_18_4DBSM_tissue_validation_with_categories_cell_ctrl_ORSml_fixed.tsv",header=T)
raw2<-raw[raw$CATEGORY!="MOSAIC"&raw$REF_COUNT+raw$ALT_COUNT>30,]
raw5<-raw[raw$REF_COUNT+raw$ALT_COUNT>30,]
raw7<-raw5[raw5$REF_COUNT+raw$ALT_COUNT>30,]

bks=c(0,0.0014,0.005,0.1,0.2,0.3)

pdf(file="2020_06_19_depth_distribution_homo_lowerCI.pdf",width=4,hight=8)
ggplot(raw5[raw5$CATEGORY=="REF_HOMO",],aes(x=sqrt(LOWER_CI),fill=CATEGORY))+
	geom_histogram(bins=201,fill="grey")+
	scale_x_continuous(breaks=sqrt(bks),labels=bks)+
	geom_vline(xintercept=sqrt(0.0014),linetype=2)+
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

bks2=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1)

pdf(file="2020_06_19_depth_distribution_het_upperCI.pdf",width=4,height=8)
ggplot(raw5[raw5$CATEGORY=="HET",],aes(x=sqrt(UPPER_CI),fill=CATEGORY))+
	geom_histogram(bins=201)+
	scale_x_continuous(breaks=sqrt(bks2),labels=bks2)+
#	geom_vline(xintercept=sqrt(0.0014),linetype=2)+
	geom_vline(xintercept=sqrt(0.5),linetype=2)+
	geom_vline(xintercept=sqrt(0.4),linetype=2,col="pink")+
	geom_vline(xintercept=sqrt(0.6),linetype=2,col="pink")+
	geom_vline(xintercept=sqrt(0.35),linetype=2,col="blue")+
	geom_vline(xintercept=sqrt(0.65),linetype=2,col="blue")+
	labs(title="Het and homo AF distribution (sqr-t), after removing the 7 variants")+
#	xlim(0,1)+
	ylim(0,2000)+
	theme_classic()
dev.off()
