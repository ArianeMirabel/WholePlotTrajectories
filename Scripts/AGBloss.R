load("DB/AGBlostDB")

init<-AllYears[which(AllYears[,"campagne"]==1984),]
logg<-AllYears[which(AllYears[,"campagne"]==1988),]

#### get data damage for Tortue and Manaré, using Guyafor
#############################################################################

#  This calculation uses the methods of Camilia but with the allometric equation fo Rejou-Mechain 2017 for AGB
# and the tarif de cubage de l'Est for vol.ext
# The methods calculate the ABG before logging and for all inventory after logging 
# then dAGB = AGB initiale - AGB minimum i nthe 4 year after logging
# AGB damage = dAGB - AGB extracted
# uses a C content in AGB of 0.5
# and average wd bellow

### open libray and set values
library(data.table)
library(BIOMASS)

wd=0.736 # average Wodd density for Vol.ext
 #coord=c(5.18,52.53) paracou latitude and longitude
AGBi<-unlist(lapply(1:12,function(plo){
Plo<-init[which(init[,"n_parcelle"]==plo),]
return(sum(computeAGB(D=Plo[,"circonf"],WD=rep(wd,nrow(Plo)),coord=c(mean(Plo[,"Lon"]),mean(Plo[,"Lat"])))))
}))
names(AGBi)<-1:12

AGBp<-unlist(lapply(1:12,function(plo){
  Plo<-logg[which(logg[,"n_parcelle"]==plo),]
  return(sum(computeAGB(D=Plo[,"circonf"],WD=rep(wd,nrow(Plo)),coord=c(mean(Plo[,"Lon"]),mean(Plo[,"Lat"])))))
}))
names(AGBp)<-1:12

AGBloss<-((AGBi-AGBp)/AGBi)*100
AGBloss<-AGBloss[c(1,6,11,2,7,9,3,5,10,4,8,12)]
AGBloss<-cbind(AGBloss,names(AGBloss),rep(c(1,2,3,4),each=3),rep(c("darkolivegreen2","deepskyblue2","darkorange1","red2"),each=3))
colnames(AGBloss)<-c("AGB","plot","treat","col")

plot(c(1,4),c(min(as.numeric(AGBloss[,"AGB"])),max(as.numeric(AGBloss[,"AGB"]))))
apply(AGBloss,1,function(li){
  points(as.numeric(li[3]),as.numeric(li[1]),col=li[4],pch=19)
})

AGBloss<-as.matrix(apply(AGBloss[,1:2],2,as.numeric))
save(AGBloss,file="DB/LostAGB")
