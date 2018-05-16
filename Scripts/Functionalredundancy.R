library("MASS");library("plotrix")
source("Scripts/Vernacular_handle.R")           # Script for vernacular names replacement
source("Scripts/TraitsMiceFilling.R")            # Traits Gap-filling

load("DB/Alpha_Plots")                  # Alpha vector for vernacular replacement
load("DB/Paracou_R_Subdivided_ok")      # Paracou dataset, lists by inventoried year, sublists by plot, etc...
dates<-sort(names(LivingStand_all)) # dates considered
dates<-dates[which(!dates%in%c("1996","1998","2000","2002","2004","2006","2008","2010","2012","2014","2016","2017"))]

Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")       # Brigde functional traits dataset
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="") # Paracou15 seed and Hmax dataset
TraitsName<-c("L_thickness","L_chloro","L_toughness","SLA","WD","Bark_thick","Hmax") # Considered traits list

InventorySp<-do.call(rbind,lapply(LivingStand_all,function(yr){
  ret<-do.call(rbind,yr)[,c("Famille","Genre","name")]
  ret<-ret[which(!grepl("Indet.",ret[,"name"])),]
  return(ret[which(!duplicated(ret)),])}))
InventorySp<-InventorySp[which(!duplicated(InventorySp)),]

# Treatments definition
T0<-c(1,6,11); T1<-c(2,7,9); T2<-c(3,5,10); T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("Control","T1","T2","T3")

ColorsTr<-c("darkolivegreen2","deepskyblue2","darkorange1","red2")

smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}

# For the species with one individual: reproduce a gaussian kernel,
# centered in the individual location and with sd like the mean sd of all species
#MeanSd<-mean(unlist(lapply(unique(Bigacp[,"name"]),function(Sp){acp<-Bigacp[which(Bigacp[,"name"]==Sp),4:5]
  #if(nrow(acp)!=1){return(mean(apply(acp,2,sd)))}})))
MeanSd<-0.7

## Overall community redundancy

Nrep<-50
RedundancyTraj<-lapply(1:Nrep,function(rep){
  
  traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
  
  ACP<-dudi.pca(scale(traits_filled[,TraitsName]),scannf=FALSE,nf=2)
  Bigacp<-merge(ACP$li,traits_filled[,c("Family","Genus","name" )],by="row.names")[,c("Family","Genus","name","Axis1","Axis2")]

  xgen<-c(min(Bigacp[,"Axis1"]),max(Bigacp[,"Axis1"]))
  ygen<-c(min(Bigacp[,"Axis2"]),max(Bigacp[,"Axis2"]))
  
  SpTdp<-lapply(unique(Bigacp[,"name"]),function(Sp){
    acp<-Bigacp[which(Bigacp[,"name"]==Sp),4:5]
    if(nrow(acp)!=1){
      method1<-method2<-"ste"
      if(tryCatch(width.SJ(acp[,"Axis1"]), error=function(e) "error")=="error")
      {method1<-"dpi"}
      if(tryCatch(width.SJ(acp[,"Axis2"]), error=function(e) "error")=="error")
      {method2<-"dpi"}
      f<-kde2d(acp[,"Axis1"], acp[,"Axis2"], n = 100,lims=c(xgen,ygen),
               h = c(width.SJ(acp[,"Axis1"],method=method1),
                     width.SJ(acp[,"Axis2"],method=method2)) )
      z<-f$z
      colnames(z)<-1:100;rownames(z)<-1:100
    }
    if(nrow(acp)==1){
      z<-outer(seq(xgen[1],xgen[2],length=100),seq(ygen[1],ygen[2],length=100), 
               function(x,y) dnorm(x,acp[,"Axis1"],MeanSd)*dnorm(y,acp[,"Axis2"],MeanSd))
      colnames(z)<-1:100;rownames(z)<-1:100
    }
    return(z)
  })
  names(SpTdp)<-unique(Bigacp[,"name"])
  
  Matrep<-lapply(1:12,function(p){
    
    Redun<-lapply(dates,function(y){
      
      temp<-LivingStand_all[[which(names(LivingStand_all)==y)]]
      if(!any(names(temp)==p)){return(NA)}
      if(any(names(temp)==p)){
        temp<-temp[[which(names(temp)==p)]]
        temp<-temp[!duplicated(temp),]
        temp<-Replacement(temp,Alpha=alphas_plot[[which(names(alphas_plot)==p)]])
        temp<-as.ProbaVector(tapply(temp,temp,length))}
      
      FunSpc<-SpTdp[intersect(names(temp),names(SpTdp))]
      temp<-temp[intersect(names(temp),names(SpTdp))]
      
      FunSpc<-lapply(names(FunSpc),function(sp){
        return(FunSpc[[sp]]*temp[which(names(temp)==sp)])})
      FunSpc<-array(unlist(FunSpc),dim=c(nrow(FunSpc[[1]]),ncol(FunSpc[[1]]),length(FunSpc)),
                    dimnames=list(rownames(FunSpc[[1]]),colnames(FunSpc[[1]]),intersect(names(temp),names(SpTdp))))
      FunSpc<-apply(FunSpc,c(1,2),function(pix){return(sum(pix))})#-max(pix)
      return(sum(FunSpc)-1)})
    Redun<-unlist(Redun)
    names(Redun)<-dates
    return(Redun)
  })
  Matrep<-do.call(rbind,Matrep)
  rownames(Matrep)<-1:12
  return(Matrep)
})

save(RedundancyTraj,file="DB/Redundancy")


## Redundancy restricted to the initial community functional space

Nrep<-2

# Pb species, "Pouteria_sp.42CAY-ATDN", "Inga_nobilis"
Sptdp<-lapply(1:Nrep,function(rep){
  traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
  
  ACP<-dudi.pca(scale(traits_filled[,TraitsName]),scannf=FALSE,nf=2)
  Bigacp<-merge(ACP$li,traits_filled[,c("Family","Genus","name" )],by="row.names")[,c("Family","Genus","name","Axis1","Axis2")]
  
  xgen<-c(min(Bigacp[,"Axis1"]),max(Bigacp[,"Axis1"]))
  ygen<-c(min(Bigacp[,"Axis2"]),max(Bigacp[,"Axis2"]))
  
  sptdp<-lapply(unique(Bigacp[,"name"]),function(Sp){
    #sptdp<-lapply(c("Pouteria_sp.42CAY-ATDN", "Inga_nobilis"),function(Sp){
    
    acp<-Bigacp[which(Bigacp[,"name"]==Sp),4:5]
    if(nrow(acp)!=1){
      method1<-method2<-"ste"
      if(tryCatch(width.SJ(acp[,"Axis1"]), error=function(e) "error")=="error")
      {method1<-"dpi"}
      if(tryCatch(width.SJ(acp[,"Axis2"]), error=function(e) "error")=="error")
      {method2<-"dpi"}
      f<-kde2d(acp[,"Axis1"], acp[,"Axis2"], n = 100,lims=c(xgen,ygen),
               h = c(width.SJ(acp[,"Axis1"],method=method1),
                     width.SJ(acp[,"Axis2"],method=method2)) )
      z<-f$z
      colnames(z)<-1:100;rownames(z)<-1:100
    }
    if(nrow(acp)==1){
      z<-outer(seq(xgen[1],xgen[2],length=100),seq(ygen[1],ygen[2],length=100), 
               function(x,y) dnorm(x,acp[,"Axis1"],MeanSd)*dnorm(y,acp[,"Axis2"],MeanSd))
      colnames(z)<-1:100;rownames(z)<-1:100
    }
    if(is.null(z)){print(Sp)}
    return(z)
  })
  names(sptdp)<-unique(Bigacp[,"name"])
  #names(sptdp)<-c("Pouteria_sp.42CAY-ATDN", "Inga_nobilis")
  return(sptdp)
})

SpTdp<-lapply(names(Sptdp[[1]]),function(Sp){
  ret<-lapply(Sptdp,function(rep){return(rep[[Sp]])})
  if(any(lapply(ret,is.null))){print(Sp)}
  ret<-array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),length(ret)),
             dimnames=list(rownames(ret[[1]]),colnames(ret[[1]]),1:length(ret)))
  return(apply(ret,c(1,2),median))})
names(SpTdp)<-names(Sptdp[[1]])

RefSpaces<-lapply(1:12,function(p){   
  temp<-LivingStand_all[[which(names(LivingStand_all)==1984)]]
  temp<-temp[[which(names(temp)==p)]]
  temp<-temp[!duplicated(temp),]
  temp<-Replacement(temp,Alpha=alphas_plot[[which(names(alphas_plot)==p)]])
  temp<-as.ProbaVector(tapply(temp,temp,length))
  
  Funref<-SpTdp[intersect(names(temp),names(SpTdp))]
  
  Funref<-lapply(names(Funref),function(sp){
    return(Funref[[sp]]*temp[which(names(temp)==sp)])})
  Funref<-array(unlist(Funref),dim=c(100,100,length(Funref)),
                dimnames=list(1:100,1:100,intersect(names(temp),names(SpTdp))))
  Funref<-apply(Funref,c(1,2),function(pix){return(sum(pix))})#-max(pix)
  SpaceRef<-contourLines(1:100,1:100,Funref,levels=c(max(Funref),0.01))
  #image(Funref);rect(SpaceRef[2,1]*0.01,SpaceRef[1,1]*0.01,SpaceRef[2,2]*0.01,SpaceRef[1,2]*0.01,border=1)
  #contour(1:100,1:100,Funref,levels=0.01)
  #Sys.sleep(1)
  SpaceRefx<-unlist(lapply(SpaceRef,function(lev){
    return(c(min(lev$x),max(lev$x)))
  }))
  SpaceRefy<-unlist(lapply(SpaceRef,function(lev){
    return(c(min(lev$y),max(lev$y)))
  }))
  SpaceRef<-rbind(c(min(SpaceRefy),max(SpaceRefy)),c(min(SpaceRefx),max(SpaceRefx)))
  rownames(SpaceRef)<-c("y","x");colnames(SpaceRef)<-c("min","max")
  SpaceRef[,"min"]<-floor(SpaceRef[,"min"]);SpaceRef[,"max"]<-ceiling(SpaceRef[,"max"])
  
  
  return(SpaceRef)})
names(RefSpaces)<-1:12

RedundancyTraj_restricted<-lapply(1:Nrep,function(rep){
Matrep<-lapply(1:12,function(p){
  
  Redun<-lapply(dates,function(y){
    
temp<-LivingStand_all[[which(names(LivingStand_all)==y)]]
    if(!any(names(temp)==p)){return(NA)}
     if(any(names(temp)==p)){
        temp<-temp[[which(names(temp)==p)]]
        temp<-temp[!duplicated(temp),]
        temp<-Replacement(temp,Alpha=alphas_plot[[which(names(alphas_plot)==p)]])
        temp<-as.ProbaVector(tapply(temp,temp,length))}

FunSpc<-SpTdp[intersect(names(temp),names(SpTdp))]
temp<-temp[intersect(names(temp),names(SpTdp))]

FunSpc<-lapply(names(FunSpc),function(sp){
  return(FunSpc[[sp]]*temp[which(names(temp)==sp)])})
FunSpc<-array(unlist(FunSpc),dim=c(nrow(FunSpc[[1]]),ncol(FunSpc[[1]]),length(FunSpc)),
              dimnames=list(rownames(FunSpc[[1]]),colnames(FunSpc[[1]]),intersect(names(temp),names(SpTdp))))
FunSpc<-apply(FunSpc,c(1,2),function(pix){return(sum(pix))})#-max(pix)
FunSpc<-FunSpc[RefSpaces[[p]][2,1]:RefSpaces[[p]][2,2],RefSpaces[[p]][1,1]:RefSpaces[[p]][1,2]]
return(sum(FunSpc)-1)})
Redun<-unlist(Redun)
names(Redun)<-dates
return(Redun)
})
Matrep<-do.call(rbind,Matrep)
rownames(Matrep)<-1:12
return(Matrep)
})

save(RedundancyTraj_restricted,file="DB/Redundancy_restricted")





#####################################################
traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)

ACP<-dudi.pca(scale(traits_filled[,TraitsName]),scannf=FALSE,nf=2)
ACP_Indiv<-merge(ACP$li,traits_filled[,c("Family","Genus","name" )],by="row.names")[,c("Family","Genus","name","Axis1","Axis2")]
ACP_Indiv[which(ACP_Indiv[,"name"]=="Lecythis_corrugata subsp. corrugata"),"name"]<-"Lecythis_corrugata"
ACP_traits<-ACP$co
Eigen<-as.vector(ACP$eig*100/sum(ACP$eig))
DataACP<-list(ACP_Indiv,ACP_traits,Eigen)
names(DataACP)<-c("Indiv","Traits","Eigen")

save(DataACP,file="DB/FunctionalPCA")
#Bigacp<-aggregate(Bigacp[,c("Axis1","Axis2")],list(Bigacp$name),median) # Median values for species level data

save(Bigacp,file="DB/FunctionalPCA")
BigacpIndiv<-merge(ACP$li,traits_filled[,c("Family","Genus","name" )],by="row.names")[,c("Family","Genus","name","Axis1","Axis2")]


plot(DataAcp_traits[,"Comp1"],DataAcp_traits[,"Comp2"],type="n",xlab="",ylab="",frame.plot=F,lwd=0)
text(DataAcp_traits[,"Comp1"],DataAcp_traits[,"Comp2"],labels=rownames(DataAcp_traits))
abline(h=0,v=0)
mtext("Axis 1",1,at=0,line=2)
mtext(paste(round(vp[1]),"% of variance",sep=""),1,at=0,line=3,cex=0.8)
mtext("Axis 2",2,at=0,las=1,line=2.5)
mtext(paste(round(vp[2]),"% of variance",sep=""),2,at=-1.5,las=1,line=2.5,cex=0.8)
arrows(0,0,x1=DataAcp_traits[,"Comp1"], y1=DataAcp_traits[,"Comp2"], col="grey", length=0.1)
title(main="Traits in the main PCA plan",adj=0,cex.main=0.9)

