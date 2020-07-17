#Rank plot of variants restricted in one hemisphere

library(ggplot2)

raw2<-raw[raw$Group=="Lower",]

quantile(raw2$Estimated_number,probs=c(0,0.025,0.05,0.95,0.975,1))
##     0%    2.5%      5%     95%   97.5%    100% 
##   86.0   219.0   296.6 47569.8 57176.2 66441.0
   
ggplot(raw2,aes(x=Rank,y=log10(Estimated_number)))+
	geom_point()+
	geom_segment(aes(x = 0, y = 0, xend = 0, yend = 5),color="black")+
	theme_classic()+
	geom_hline(yintercept=log10(86),linetype=2,col="grey")+
	theme(axis.line=element_blank())
