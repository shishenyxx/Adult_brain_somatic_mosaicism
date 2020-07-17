##plotting for the correlation between AF measured in WGS and MPAS

library(ggplot2)

raw<-read.delim(file="2020_06_18_4DBSM_positive_variants_with_categories_WGS_and_AmpliSeq.txt",header=T)

names(raw)
## [1] "ID"                  "CHROM.POS.REF.ALT"   "IDCHROM.POS.REF.ALT"
## [4] "FLAG_NOISE"          "CHROM"               "POS"                
## [7] "REF"                 "ALT"                 "DEPTH"              
##[10] "REF_COUNT"           "ALT_COUNT"           "WGS_MAF"            
##[13] "WGS_LOWER_CI"        "WGS_UPPER_CI"        "CATEGORY"           
##[16] "AmpliSeq_DEPTH"      "AmpliSeq_REF_COUNT"  "AmpliSeq_ALT_COUNT" 
##[19] "AmpliSeq_MAF"        "AmpliSeq_LOWER_CI"   "AmpliSeq_UPPER_CI"  
##[22] "SET_ONE_TISSUE"      "SET_ONE_REGION"      "SET_CORTEX_ONLY"    
##[25] "SET_BRAIN_ONLY"      "SET_KIDNEY_ONLY"     "SET_LEFT_ONLY"      
##[28] "SET_RIGHT_ONLY"      "SET_IN_CORTEX"       "SET_IN_BRAIN"       
##[31] "SET_IN_CEREBELLUM"   "SET_IN_HEART"        "SET_IN_LIVER"       
##[34] "SET_IN_KIDNEY"       "SET_LEFT_ONLY_CTX"   "SET_RIGHT_ONLY_CTX" 
##[37] "CAT_LABEL"           "LEFT_RIGHT"          "Ctx_PF"             
##[40] "Ctx_F"               "Ctx_O"               "Ctx_T"              
##[43] "Ctx_P"               "no_Region"   

raw<-raw[raw$FLAG_NOISE=="FALSE",]

bks2=c(0,0.0014,0.05,0.1,0.20,0.3,0.4,0.5,0.6)

pdf(file="2020_06_18_MAF_in_WGS_and_AmpliSeq_no_noise.pdf",width=10,height=7)

ggplot(raw,aes(x=sqrt(WGS_MAF),xmax=sqrt(WGS_UPPER_CI),xmin=sqrt(WGS_LOWER_CI),y=sqrt(AmpliSeq_MAF),ymax=sqrt(AmpliSeq_UPPER_CI),ymin=sqrt(AmpliSeq_LOWER_CI),col=CAT_LABEL,fill=CAT_LABEL))+
	geom_point(size=0.5)+
	geom_pointrange(size=0.5)+
	geom_abline(slope=1,intercept=0,linetype=2)+
	scale_y_continuous(breaks=sqrt(bks2),labels=bks2)+
	scale_x_continuous(breaks=sqrt(bks2),labels=bks2)+
#	annotate("text", x = sqrt(0.2), y = sqrt(0.7), label = bquote(*(R^2)*" = 0.66"))+
	annotate("text", x = sqrt(0.2), y = sqrt(0.7), label = bquote(~ R^2~"=0.68")  )+
	labs(title="Category Labels",x=bquote("sqr-t ("~MAF[WGS]~")"),y=bquote("sqr-t ("~MAF[AmpliSeq]~")"))+
	theme_classic()
 
dev.off()
