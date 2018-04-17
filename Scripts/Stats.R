load("DB/Turnover_toInit")
treats<-cbind(c(1,6,11,2,7,9,3,5,10,4,8,12),rep(0:3,each=3))
treats<-treats[order(treats[,1]),]
colnames(treats)<-c("plot","treat"); rownames(treats)<-treats[,"plot"]

Maxs<-cbind(apply(Turn,1,max),treats[,"treat"])
colnames(Maxs)<-c("Max","treat")

cor(Maxs[,"Max"],Maxs[,"treat"],method="spearman")

##################################################"

Traits<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")
Traits["name"]<-as.factor(paste(Traits$Genus,"_",Traits$species,sep=""))
traits<-Traits1[,c("Family","Genus","species","name","bar_code","thickness","SPAD","toughness","dry_mass","traits_surf_area","ind_surf_area",
                   "sapwood_dens","moisture","bark_thick")]
colnames(traits)<-c("Family","Genus","species","name","bar_code","L_thickness","L_chloro","L_toughness","L_DryMass","LA","SLA","WD",
                    "moisture","Bark_thick")
##SLA calc for empty values with LA available
traits[which(!is.na(traits[,"LA"])),"SLA"]<-
  traits[which(!is.na(traits[,"LA"])),"LA"]/traits[which(!is.na(traits[,"LA"])),"L_DryMass"]

traitsName<-c("L_thickness","L_chloro","L_toughness","L_DryMass","LA","SLA","WD","moisture","Bark_thick")

allSd<-apply(traits[,traitsName],2,function(x){return(sd(x, na.rm=T))})

SpSd<-aggregate(traits[,traitsName],list(traits$name),function(x){return(sd(x, na.rm=T))})
SpSd<-apply(SpSd[,traitsName],2,function(x){return(sd(x, na.rm=T))})

GenSd<-aggregate(traits[,traitsName],list(traits$Genus),function(x){return(sd(x, na.rm=T))})
GenSd<-apply(GenSd[,traitsName],2,function(x){return(sd(x, na.rm=T))})

Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="")
Traits2<-Traits2[which(apply(Traits2,1,function(li){return(any(!is.na(li)))})),]
Traits2["name"]<-sub(" ","_",Traits2[,"Name"])
traits2<-Traits2[,c("name","Hauteur","Masse")]
colnames(traits2)<-c("name","Hmax","S_mass")
traits2[which(traits2[,"name"]=="Sterculia_speciosa"),"Hmax"]<-43

traits<-merge(traits2,traits,by="name")[,c("Hmax","L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Bark_thick")]
### To be continued
