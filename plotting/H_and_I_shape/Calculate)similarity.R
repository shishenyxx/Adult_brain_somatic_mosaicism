#After spliting the whole plain x∈(-1,1) and y∈(-1,1) in to 25 grids by assigning each coordinate into group =IF(A2<-0.2,IF(A2<-0.6,1,2),IF(A2>0.2,IF(A2>0.6,5,4),3))+IF(B2<-0.2,IF(B2<-0.6,20,15),IF(B2>0.2,IF(B2>0.6,0,5),10))

#Define number of variants in each grid.

library(lsa)

H_shape<-c(1,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0,1)
Anti_H<-c(1,1,1,1,1,0,1,1,1,0,0,0,1,0,0,0,1,1,1,0,1,1,1,1,1)
ID01_sorted<-c(15,0,0,0,10,8,1,0,0,10,8,2,12,6,9,8,2,4,1,9,13,0,1,2,11)
ID02_03_04<-c(22,0,3,1,95,4,6,6,9,8,4,11,37,12,9,5,4,3,3,9,66,2,0,2,130)
ID01_bulk<-c(36,0,1,0,17,5,2,1,2,5,0,3,9,3,5,3,5,3,3,4,46,1,0,0,33)
ID02_sorted<-c(5,0,2,2,9,2,3,5,0,1,1,2,14,4,2,3,0,0,2,2,7,0,1,0,10)
ID_02<-c(13,0,2,0,14,2,2,3,3,1,1,4,15,3,2,2,1,3,1,1,12,2,0,0,8)
ID_03<-c(3,1,1,1,59,1,4,0,3,2,1,4,9,6,2,0,1,0,0,4,3,0,0,2,119)
ID_04<-c(6,0,0,0,22,1,0,3,3,5,2,3,13,3,5,3,2,0,2,4,51,0,0,0,3)

cosine(H_shape,ID01_sorted)
cosine(Anti_H,ID01_sorted)

cosine(H_shape,ID01_bulk)
cosine(Anti_H,ID01_bulk)

cosine(H_shape,ID02_03_04)
cosine(Anti_H,ID02_03_04)

cosine(H_shape,ID02_sorted)
cosine(Anti_H,ID02_sorted)

cosine(H_shape,ID_02)
cosine(Anti_H,ID_02)

cosine(H_shape,ID_03)
cosine(Anti_H,ID_03)

cosine(H_shape,ID_04)
cosine(Anti_H,ID_04)


#calculate the same R values and cosine similarities for all the 1000 simulations

raw<-read.csv(file="ID_01_4tissue_bulk_simulation.txt",header=F)

output_new<-NA
output_new$corr_H<-NA
output_new$corr_anti_H<-NA
output_new$R<-NA

for(i in 1:10000)
{
output_new$corr_H[i]<-cosine(as.numeric(raw[i,]),H_shape)
output_new$corr_anti_H[i]<-cosine(as.numeric(raw[i,]),Anti_H)
output_new$R[i]<-output_new$corr_H[i]/output_new$corr_anti_H[i]
}

write.csv(output_new, file="ID_01_4tissue_bulk_summary_sim.csv")



raw<-read.csv(file="ID_01_sorted_simulation.txt",header=F)

output_new<-NA
output_new$corr_H<-NA
output_new$corr_anti_H<-NA
output_new$R<-NA

for(i in 1:10000)
{
output_new$corr_H[i]<-cosine(as.numeric(raw[i,]),H_shape)
output_new$corr_anti_H[i]<-cosine(as.numeric(raw[i,]),Anti_H)
output_new$R[i]<-output_new$corr_H[i]/output_new$corr_anti_H[i]
}

write.csv(output_new, file="ID_01_sorted_summary_sim.csv")



raw<-read.csv(file="ID_02_03_04_simulation.txt",header=F)

output_new<-NA
output_new$corr_H<-NA
output_new$corr_anti_H<-NA
output_new$R<-NA

for(i in 1:10000)
{
output_new$corr_H[i]<-cosine(as.numeric(raw[i,]),H_shape)
output_new$corr_anti_H[i]<-cosine(as.numeric(raw[i,]),Anti_H)
output_new$R[i]<-output_new$corr_H[i]/output_new$corr_anti_H[i]
}

write.csv(output_new, file="ID_02_03_04_summary_sim.csv")



raw<-read.csv(file="ID_02_bulk_simulation.txt",header=F)

output_new<-NA
output_new$corr_H<-NA
output_new$corr_anti_H<-NA
output_new$R<-NA

for(i in 1:10000)
{
output_new$corr_H[i]<-cosine(as.numeric(raw[i,]),H_shape)
output_new$corr_anti_H[i]<-cosine(as.numeric(raw[i,]),Anti_H)
output_new$R[i]<-output_new$corr_H[i]/output_new$corr_anti_H[i]
}

write.csv(output_new, file="ID_02_bulk_summary_sim.csv")



raw<-read.csv(file="ID_02_sorted_simulation.txt",header=F)

output_new<-NA
output_new$corr_H<-NA
output_new$corr_anti_H<-NA
output_new$R<-NA

for(i in 1:10000)
{
output_new$corr_H[i]<-cosine(as.numeric(raw[i,]),H_shape)
output_new$corr_anti_H[i]<-cosine(as.numeric(raw[i,]),Anti_H)
output_new$R[i]<-output_new$corr_H[i]/output_new$corr_anti_H[i]
}

write.csv(output_new, file="ID_02_sorted_summary_sim.csv")



raw<-read.csv(file="ID_03_bulk_simulation.txt",header=F)

output_new<-NA
output_new$corr_H<-NA
output_new$corr_anti_H<-NA
output_new$R<-NA

for(i in 1:10000)
{
output_new$corr_H[i]<-cosine(as.numeric(raw[i,]),H_shape)
output_new$corr_anti_H[i]<-cosine(as.numeric(raw[i,]),Anti_H)
output_new$R[i]<-output_new$corr_H[i]/output_new$corr_anti_H[i]
}

write.csv(output_new, file="ID_03_bulk_summary_sim.csv")



raw<-read.csv(file="ID_04_bulk_bulk_simulation.txt",header=F)

output_new<-NA
output_new$corr_H<-NA
output_new$corr_anti_H<-NA
output_new$R<-NA

for(i in 1:10000)
{
output_new$corr_H[i]<-cosine(as.numeric(raw[i,]),H_shape)
output_new$corr_anti_H[i]<-cosine(as.numeric(raw[i,]),Anti_H)
output_new$R[i]<-output_new$corr_H[i]/output_new$corr_anti_H[i]
}

write.csv(output_new, file="ID_04_bulk_summary_sim.csv")
