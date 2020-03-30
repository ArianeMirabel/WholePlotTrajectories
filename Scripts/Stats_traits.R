load("DB/Paracou_R_Subdivided_ok")
Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="")
TraitsName<-c("Hmax","L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Bark_thick")
InventorySp<-do.call(rbind,lapply(LivingStand_all,function(yr){
  +   ret<-do.call(rbind,yr)[,c("Famille","Genre","name")]
  +   ret<-ret[which(!grepl("Indet.",ret[,"name"])),]
  +   return(ret[which(!duplicated(ret)),])}))
InventorySp<-InventorySp[which(!duplicated(InventorySp)),]

Traits1<-Traits1[which(apply(Traits1,1,function(li){return(any(!is.na(li)))})),]
Traits2<-Traits2[which(apply(Traits2,1,function(li){return(any(!is.na(li)))})),]

#Ajout colonne nom genre_espece pour avoir un seul identifiant par espèce
Traits1["name"]<-as.factor(paste(Traits1$Genus,"_",Traits1$species,sep=""))
Traits2["name"]<-sub(" ","_",Traits2[,"Name"])

#Renommer les colonnes,purement pratique
traits1<-Traits1[,c("Family","Genus","species","name","bar_code","thickness","SPAD","toughness","dry_mass","traits_surf_area","ind_surf_area",
                    "sapwood_dens","moisture","bark_thick")]
colnames(traits1)<-c("Family","Genus","species","name","bar_code","L_thickness","L_chloro","L_toughness","L_DryMass","LA","SLA","WD",
                     "moisture","Bark_thick")
traits2<-Traits2[,c("name","Hauteur","Masse")]
colnames(traits2)<-c("name","Hmax","S_mass")
traits2[which(traits2[,"name"]=="Sterculia_speciosa"),"Hmax"]<-43    # Anecdotique, pb d'estimation de la hauteur pour cette espèce

# Réunir les deux tables de traits en 1
traits<-merge(traits2,traits1,by="name")

traits[which(!is.na(traits[,"LA"])),"SLA"]<-
  traits[which(!is.na(traits[,"LA"])),"LA"]/traits[which(!is.na(traits[,"LA"])),"L_DryMass"]

traits<-subset(traits, traits$name %in% InventorySp$name)

traits<-apply(traits[,c("L_thickness","L_toughness","L_chloro","SLA","WD","Bark_thick")],2,function(col){ 
  tapply(col, traits$name, function(gp){mean(gp,na.rm=T)})})

NAcount<-apply(traits,2,function(col){length(which(is.na(col)))})

NAcount/nrow(traits)*100

NAcount_ind<-apply(traits[,c("L_thickness","L_toughness","L_chloro","SLA","WD","Bark_thick")],2,function(col){ 
  length(which(is.na(col)))})

NAcount_ind/nrow(traits)*100


splist<-lapply(LivingStand_all,function(yr){
  return(unlist(lapply(yr,function(plo){
    freq<-as.character(plo[which(plo[,"name"]%in%traits[,"name"]),"name"])
    freq<-sum(tapply(freq,freq,length))
    return(freq/sum(tapply(as.character(plo[,"name"]),as.character(plo[,"name"]),length)))
  })))
})
splist<-lapply(splist,function(yr){return(as.numeric(yr[which(names(yr)%in%1:12)]))})
splist<-splist[which(unlist(lapply(splist,function(yr){length(yr)!=0})))]

Pby<-c("1996","1998","2000","2002","2004","2006","2008","2010","2012","2014","2016","2017")
splist<-unlist(splist[which(!names(splist)%in%Pby)])
min(splist)
mean(splist)

