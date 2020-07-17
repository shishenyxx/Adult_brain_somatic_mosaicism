##estimation of the starting population of variants shared between left right and tissues

library(FSA)
library(ggplot2)

pdf(file="try_to_estimate_cells_hypergenomic_distribution_step5.pdf")

input$flag<-NA
input$temp<-NA

for(i in 1:nrow(input))
{
raw_input$idd<-NA
raw_input$idd<-seq(5,5000,by=5)
raw_input$hyperLB <- NA
raw_input$hyperUB <- NA
input$flag[i]=1
for(j in 1:1000){
	temp<-hyperCI(2*ceiling(raw_input$idd[j]*input$MEAN_LEFT_RIGHT[i]*2),raw_input$idd[j],ceiling(raw_input$idd[j]*input$MEAN_LEFT_RIGHT[i]*2))
#	raw_input$hyperLB[j] <- temp[ceiling(nrow(temp)/2),1]/raw_input$idd[j]
	raw_input$hyperLB[j] <- max(0,2*ceiling(raw_input$idd[j]*input$MEAN_LEFT_RIGHT[i]*2)/temp[2])
	raw_input$hyperUB[j] <- 2*ceiling(raw_input$idd[j]*input$MEAN_LEFT_RIGHT[i]*2)/temp[1]
	if(!(raw_input$hyperLB[j]>min(input$LEFT_AVG_CMAF[i]*2,input$RIGHT_AVG_CMAF[i]*2)&raw_input$hyperUB[j]<max(input$LEFT_AVG_CMAF[i]*2,input$RIGHT_AVG_CMAF[i]*2))&input$flag[i]==1){
	input$temp[i]=raw_input$idd[j]*2
	input$flag[i]=input$flag[i]*1
	}
	else{
	input$flag[i]=input$flag[i]*0
	}
}
plot <-
ggplot(as.data.frame(raw_input),aes(x=idd,y=input$MEAN_LEFT_RIGHT[i]*2,ymin=hyperLB,ymax=hyperUB))+
	geom_pointrange()+
	geom_hline(yintercept=input$LEFT_AVG_CMAF[i]*2,linetype=2,col="red")+
	geom_hline(yintercept=input$RIGHT_AVG_CMAF[i]*2,linetype=2,col="blue")+
	labs(x="Starting population",y="Cellular fractions",title=paste(input$CHR_POS_REF_ALT[i]," ",input$temp[i]))+
	theme_classic()
	
print(plot)
}
dev.off()
