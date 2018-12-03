load("DB/AGBlostDB")
load("DB/AllCirc")

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
Carr<-unlist(lapply(unique(Plo[,"n_carre"]),function(car){
  carr<-Plo[which(Plo[,"n_carre"]==car),]
return(sum(computeAGB(D=carr[,"circ_corr"],WD=rep(wd,nrow(carr)),coord=c(mean(carr[,"Lon"]),mean(carr[,"Lat"])))))
}))
names(Carr)<-paste(plo,unique(Plo[,"n_carre"]),sep=".")
return(Carr)
}))

AGBp<-unlist(lapply(1:12,function(plo){
  Plo<-logg[which(logg[,"n_parcelle"]==plo),]
  Carr<-unlist(lapply(unique(Plo[,"n_carre"]),function(car){
    carr<-Plo[which(Plo[,"n_carre"]==car),]
    return(sum(computeAGB(D=carr[,"circ_corr"],WD=rep(wd,nrow(carr)),coord=c(mean(carr[,"Lon"]),mean(carr[,"Lat"])))))
  }))
  names(Carr)<-paste(plo,unique(Plo[,"n_carre"]),sep=".")
  return(Carr)
  #return(cbind(rep(plo,length(unique(Plo[,"n_carre"]))),unique(Plo[,"n_carre"]),Carr))
}))
#colnames(AGBp)<-c("Plot","carre","AGB")

AGBloss<-((AGBi-AGBp)/AGBi)*100
c(1,6,11,2,7,9,3,5,10,4,8,12)
AGBloss<-cbind(rep(1:12,rep(4,12)),rep(1:4,12),AGBloss)
colnames(AGBloss)<-c("plot","carre","AGBloss")

save(AGBloss,file="DB/LostAGBcarre")

