load("DB/ReplacementTraj_ForGraphs")

CompTaxo<-CompleteTaxo

treats<-cbind(c(1,6,11,2,7,9,3,5,10,4,8,12),rep(0:3,each=3))
treats<-treats[order(treats[,1]),]
colnames(treats)<-c("plot","treat"); rownames(treats)<-treats[,"plot"]

for(q in 1:3){     
  Toplot<-lapply(CompTaxo,function(tr){return(tr[,,,q])})
  Toplot<-lapply(Toplot,function(toplot){return(toplot[,which(colnames(toplot)>=1989),])})
  Toplot<-lapply(Toplot,function(tr){
    ret<-lapply(1:dim(tr)[3],function(rep){return(apply(tr[,,rep],2,function(col){col<-col-tr[,1,rep]}))})#
    ret<-array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),length(ret)),
               dimnames=list(rownames(ret[[1]]),as.numeric(colnames(ret[[1]]))-1984,1:length(ret)))
    return(apply(apply(ret,c(1,2),median),1,max))
    })
  Ret<-as.data.frame(unlist(Toplot))
  Ret$plot<-rownames(Ret)
  Ret<-Ret[order(as.numeric(Ret$plot)),]
  Ret<-merge(Ret,treats,by="plot")
  assign(c("Richness","Shannon","Simpson")[q],Ret)
}

colnames(Richness)<-c("plot","Max","treat")
cor(Richness[,"Max"],Richness[,"treat"],method="spearman")
colnames(Shannon)<-c("plot","Max","treat")
cor(Shannon[,"Max"],Shannon[,"treat"],method="spearman")
colnames(Simpson)<-c("plot","Max","treat")
cor(Simpson[,"Max"],Simpson[,"treat"],method="spearman")


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
