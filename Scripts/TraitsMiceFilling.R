library(mice)

Traits_filling<-function(Traits1,Traits2){
  Traits1<-Traits1[which(apply(Traits1,1,function(li){return(any(!is.na(li)))})),]
  Traits2<-Traits2[which(apply(Traits2,1,function(li){return(any(!is.na(li)))})),]
  
## Traits correlation?
Traits1["name"]<-as.factor(paste(Traits1$Genus,"_",Traits1$species,sep=""))
Traits2["name"]<-sub(" ","_",Traits2[,"Name"])

traits1<-Traits1[,c("Family","Genus","species","name","bar_code","thickness","SPAD","toughness","dry_mass","traits_surf_area","ind_surf_area",
                    "sapwood_dens","moisture","bark_thick")]
colnames(traits1)<-c("Family","Genus","species","name","bar_code","L_thickness","L_chloro","L_toughness","L_DryMass","LA","SLA","WD",
                     "moisture","Bark_thick")
traits2<-Traits2[,c("name","Hauteur","Masse")]
colnames(traits2)<-c("name","Hmax","S_mass")

traits<-merge(traits2,traits1,by="name")

##SLA calc for empty values with LA available
traits[which(!is.na(traits[,"LA"])),"SLA"]<-
  traits[which(!is.na(traits[,"LA"])),"LA"]/traits[which(!is.na(traits[,"LA"])),"L_DryMass"]
Seltraits<-c("Hmax","S_mass","L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","moisture","Bark_thick")
traits<-traits[,c("Family","Genus","name","bar_code",Seltraits)]

##########################  Only some missing values: MICE regression to keep all available information
### Non entire missing lines
traitsPartial<-traits[which(apply(traits[,Seltraits],1,function(li){return(any(!is.na(li)))})),]
traitsPartial<-traitsPartial[which(apply(traitsPartial[,Seltraits],1,function(li){return(any(is.na(li)))})),]
traitsComp<-traits[which(apply(traits[,Seltraits],1,function(li){return(!any(is.na(li)))})),]
traitsEmpty<-traits[which(apply(traits[,Seltraits],1,function(li){return(!any(!is.na(li)))})),]

traitsPartial_gen<-lapply(unique(traitsPartial[,"Genus"]),function(gen){return(traitsPartial[which(traitsPartial[,"Genus"]==gen),])})

#sort(unique(unlit(lapply(traitsPartial_gen,nrow))))
#For some genus, only 1 species: we group by family the genuses where there is not enough samples: at least 20 individuals
traitsPartial_fam<-do.call(rbind,traitsPartial_gen[which(unlist(lapply(traitsPartial_gen,nrow))<=20)])

traitsPartial_gen<-traitsPartial_gen[which(unlist(lapply(traitsPartial_gen,nrow))>20)]
traitsPartial_gen<-lapply(traitsPartial_gen,function(sub){
  ret<-complete(mice(sub[,Seltraits],printFlag=F))
  return(cbind(sub[,c("Family","Genus","name","bar_code")],ret))})

traitsPartial_fam<-rbind(traitsPartial_fam,
                  do.call(rbind,traitsPartial_gen[which(unlist(lapply(traitsPartial_gen,function(sub){anyNA(sub[,"WD"])})))]))
traitsPartial_gen<-do.call(rbind,
                  traitsPartial_gen[which(!unlist(lapply(traitsPartial_gen,function(sub){anyNA(sub[,"WD"])})))])

traitsPartial_fam<-lapply(unique(traitsPartial_fam[,"Family"]),function(fam){
  return(traitsPartial_fam[which(traitsPartial_fam[,"Family"]==fam & traitsPartial_fam[,"Family"]!="Arecaceae"),])})
traitsPartial_all<-do.call(rbind,traitsPartial_fam[which(unlist(lapply(traitsPartial_fam,nrow))<=20)])
traitsPartial_all<-rbind(traitsPartial_all,traitsPartial[which(traitsPartial[,"Family"]=="Arecaceae"),])
traitsPartial_fam<-traitsPartial_fam[which(unlist(lapply(traitsPartial_fam,nrow))>20)]

traitsPartial_fam<-lapply(traitsPartial_fam,function(sub){
  ret<-complete(mice(sub,printFlag=F))
  return(cbind(sub[,c("Family","Genus","name","bar_code")],ret))})

idScarce<-as.character(traitsPartial_all[,"bar_code"],
   do.call(rbind,traitsPartial_fam[which(unlist(lapply(traitsPartial_fam,function(sub){anyNA(sub[,"WD"])})))])[,"bar_code"])
comScarce<-cbind(traits[,c("Family","Genus","name","bar_code")],
                 complete(mice(traits[,Seltraits],printFlag=F)))
comScarce<-comScarce[which(comScarce[,"bar_code"]%in%idScarce),]
traitsPartial_fam<-do.call(rbind,
                    traitsPartial_fam[which(!unlist(lapply(traitsPartial_fam,function(sub){anyNA(sub[,"WD"])})))])

traitsPartial_filled<-rbind(traitsPartial_gen,traitsPartial_fam,comScarce)

## No necessary predictor selection: less than 15 vriables
#### it is not possible to use multilevel model, because some Families have only one sample

##### For completely missing lines: sample the traits values from the Melaine's algorithm
if(nrow(traitsEmpty)!=0){
traitsEmpty_filled<-do.call(rbind,lapply(1:nrow(traitsEmpty),function(li){
  sp<-traitsEmpty[li,]
  Tosample<-traitsComp[which(traitsComp[,"name"]==as.character(sp["name"])),8:15]
  if(nrow(Tosample)!=0){sp[8:15]<-Tosample[sample(rownames(Tosample),1),]}
  
  if(nrow(Tosample)==0){
    
    Tosample<-traitsComp[which(traitsComp[,"Genus"]==as.character(sp["Genus"])),8:15]
    if(nrow(Tosample)!=0){sp[8:15]<-Tosample[sample(rownames(Tosample),1),]}
    
    if(nrow(Tosample)==0){
      Tosample<-traitsComp[which(traitsComp[,"Family"]==as.character(sp["Family"])),8:15]
      
      if(nrow(Tosample)!=0){sp[8:15]<-Tosample[sample(rownames(Tosample),1),]}
    if(nrow(Tosample)==0){sp[8:15]<-traitsComp[sample(rownames(traitsComp),1),8:15]}
  }}
  return(sp)
}))
traits_filled<-rbind(traitsEmpty_filled,traitsPartial_filled,traitsComp)
}

if(nrow(traitsEmpty)==0){traits_filled<-rbind(traitsPartial_filled,traitsComp)}

traits_filled<-traits_filled[order(sort(traits_filled[,"name"])),which(colnames(traits_filled)!="moisture")]
traits_filled<-traits_filled[which(!is.na(traits_filled[,"Hmax"])),]
return(traits_filled)
}
