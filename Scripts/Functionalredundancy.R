library("MASS")
source("Scripts/Vernacular_handle.R")           # Script for vernacular names replacement
source("Scripts/TraitsMiceFilling.R")                                   # Traits Gap-filling

load("DB/Alpha_Plots")                  # Alpha vector for vernacular replacement
load("DB/Paracou_R_Subdivided_ok")      # Paracou dataset, lists by inventoried year, sublists by plot, etc...
dates<-sort(names(LivingStand_all))[-1] # dates considered

Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")       # Brigde functional traits dataset
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="") # Paracou15 seed and Hmax dataset
TraitsName<-c("L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Bark_thick","Hmax") # Considered traits list

InventorySp<-do.call(rbind,lapply(LivingStand_all,function(yr){
  ret<-do.call(rbind,yr)[,c("Famille","Genre","name")]
  ret<-ret[which(!grepl("Indet.",ret[,"name"])),]
  return(ret[which(!duplicated(ret)),])}))
InventorySp<-InventorySp[which(!duplicated(InventorySp)),]

traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)

### Whole community mapping, individuals: NMDS ou ACP?

BigNmds<-as.matrix(daisy(scale(traits_filled[,TraitsName]),metric="gower"))
BigNmds<-scores(metaMDS(BigNmds,distance="bray"),display="sites")
BigNmds<-merge(BigNmds,traits_filled[,c("Family","Genus","name" )],by="row.names")[,c("Family","Genus","name","NMDS1","NMDS2")]
save(BigNmds,file="DB/BigNmds")

ACP<-dudi.pca(scale(traits_filled[,TraitsName]))
Bigacp<-merge(ACP$li,traits_filled[,c("Family","Genus","name" )],by="row.names")[,c("Family","Genus","name","Axis1","Axis2")]

plot(Bigacp[,"Axis1"], Bigacp[,"Axis2"],type="n")
text(Bigacp[,"Axis1"], Bigacp[,"Axis2"],labels=Bigacp[,"name"])

### Species TPD

xgen<-c(min(Bigacp[,"Axis1"]),max(Bigacp[,"Axis1"]))
ygen<-c(min(Bigacp[,"Axis2"]),max(Bigacp[,"Axis2"]))
  
# For the species with one individual: reproduce a gaussian kernel,
# centered in the individual location and with sd like the mean sd of all species
MeanSd<-mean(unlist(lapply(unique(Bigacp[,"name"]),function(Sp){
  acp<-Bigacp[which(Bigacp[,"name"]==Sp),4:5]
  if(nrow(acp)!=1){
    return(mean(apply(acp,2,sd)))
}})))

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
  return(f$z)
  }
  if(nrow(acp)==1){
    z<-outer(seq(xgen[1],xgen[2],length=100),seq(ygen[1],ygen[2],length=100), 
              function(x,y) dnorm(x,acp[,"Axis1"],MeanSd)*dnorm(y,acp[,"Axis2"],MeanSd))
    return(z)
  }
})
names(SpTdp)<-unique(Bigacp[,"name"])

### Community TPD

# Treatments definition
T0<-c(1,6,11); T1<-c(2,7,9); T2<-c(3,5,10); T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("Control","T1","T2","T3")

Pby<-c("1996","1998","2000","2002","2004","2006","2008","2010","2012","2014","2016","2017")

ColorsTr<-c("darkolivegreen2","deepskyblue2","darkorange1","red2")

smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}

RefSpaces<-lapply(1:12,function(p){   
  temp<-LivingStand_all[[which(names(LivingStand_all)==1989)]]
  temp<-temp[[which(names(temp)==p)]]
  temp<-temp[!duplicated(temp),]
  temp<-Replacement(temp,Alpha=alphas_plot[[which(names(alphas_plot)==p)]])
  temp<-as.ProbaVector(tapply(temp,temp,length))
  
  Funref<-SpTdp[intersect(names(temp),names(SpTdp))]
  
  Funref<-lapply(names(Funref),function(sp){
    return(Funref[[sp]]*temp[which(names(temp)==sp)])})
  Funref<-array(unlist(Funref),dim=c(100,100,length(Funref)),
                dimnames=list(1:100,1:100,intersect(names(temp),names(SpTdp))))
  Funref<-apply(FunSpc,c(1,2),function(pix){return(sum(pix))})#-max(pix)
  SpaceRef<-contourLines(1:100,1:100,Funref,levels=c(max(Funref),0.05))
  SpaceRefx<-unlist(lapply(SpaceRef,function(lev){
    return(c(min(lev$x),max(lev$x)))
  }))
  SpaceRefy<-unlist(lapply(SpaceRef,function(lev){
    return(c(min(lev$y),max(lev$y)))
  }))
  SpaceRef<-rbind(c(min(SpaceRefy),max(SpaceRefy)),c(min(SpaceRefx),max(SpaceRefx)))
  rownames(SpaceRef)<-c("y","x");colnames(SpaceRef)<-c("min","max")
  return(SpaceRef)})
names(RefSpaces)<-1:12
  
Nrep<-3
RedundancyTraj<-lapply(1:Nrep,function(rep){
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
FunSpc<-array(unlist(FunSpc),dim=c(100,100,length(FunSpc)),
              dimnames=list(1:100,1:100,intersect(names(temp),names(SpTdp))))
FunSpc<-apply(FunSpc,c(1,2),function(pix){return(sum(pix))})#-max(pix)
return(sum(FunSpc)-1)})
Redun<-unlist(Redun)
names(Redun)<-dates
return(Redun)
})
Matrep<-do.call(rbind,Matrep)
Matrep<-Matrep[,which(!colnames(Matrep)%in%Pby)]
rownames(Matrep)<-1:12
return(Matrep)
})

save(RedundancyTraj,file="DB/Redundancy2")

load("DB/Redundancy")


RedundancyPlot<-function(Red){
  
  
  Red<-RedundancyTraj
  abs<-as.numeric(colnames(Red[[1]]))-1984
  Red<-lapply(Red,function(rep){return(smooth(rep,2))})
  Red<-lapply(Red,function(rep){return(apply(rep,2,function(col){return(col-rep[,5])}))})
  Red<-array(unlist(Red),dim=c(12,ncol(Red[[1]]),length(Red)),
             dimnames=list(1:12,colnames(Red[[1]]),1:length(Red)))
  Red<-lapply(c(0.025,0.5,0.975),function(quant){return(apply(Red,c(1,2),function(rep){return(quantile(rep,probs=quant))}))})
  Red<-array(unlist(Red),dim=c(12,ncol(Red[[1]]),3),
             dimnames=list(1:12,colnames(Red[[1]]),c(0.025,0.5,0.975)))
  Red<-Red[,5:ncol(Red),]
  
  plot(abs[abs>=5],Red[1,,1],type='n',ylim=c(min(Red),max(Red)),xlab="",ylab="")
  invisible(lapply(1:4,function(tr){
    toplot<-Red[which(rownames(Red)%in%treatments[[tr]]),,"0.5"]
    apply(toplot,1,function(li){
      lines(abs[abs>=5],li,col=ColorsTr[tr],lwd=2)
    })
  }))

  
  
}

RedundancyPlot(Red)















attach(geyser)
plot(duration, waiting, xlim = c(0.5,6), ylim = c(40,100))
f1 <- kde2d(duration, waiting, n = 50, lims = c(0.5, 6, 40, 100))
image(f1, zlim = c(0, 0.05))
f2 <- kde2d(duration, waiting, n = 50, lims = c(0.5, 6, 40, 100),
            h = c(width.SJ(duration), width.SJ(waiting)) )
image(f2, zlim = c(0, 0.05))
persp(f2, phi = 30, theta = 20, d = 5)

plot(duration[-272], duration[-1], xlim = c(0.5, 6),
     ylim = c(1, 6),xlab = "previous duration", ylab = "duration")
f1 <- kde2d(duration[-272], duration[-1],
            h = rep(1.5, 2), n = 50, lims = c(0.5, 6, 0.5, 6))
contour(f1, xlab = "previous duration",
        ylab = "duration", levels  =  c(0.05, 0.1, 0.2, 0.4) )
f1 <- kde2d(duration[-272], duration[-1],
            h = rep(0.6, 2), n = 50, lims = c(0.5, 6, 0.5, 6))
contour(f1, xlab = "previous duration",
        ylab = "duration", levels  =  c(0.05, 0.1, 0.2, 0.4) )
f1 <- kde2d(duration[-272], duration[-1],
            h = rep(0.4, 2), n = 50, lims = c(0.5, 6, 0.5, 6))
contour(f1, xlab = "previous duration",
        ylab = "duration", levels  =  c(0.05, 0.1, 0.2, 0.4) )


### other kernel

RadSym <- function(u)
  exp(-rowSums(u^2)/2) / (2*pi)^(ncol(u)/2)

Scott <- function(data)
  t(chol(cov(data))) * nrow(data) ^ (-1/(ncol(data)+4))

mvkde <- function(x, data, bandwidth=Scott, kernel=RadSym) {
  if(is.function(bandwidth))
    bandwidth <- bandwidth(data)
  u <- t(solve(bandwidth, t(data) - x))
  mean(kernel(u))
}

smvkde <- function(x, ...) apply(x, 1, mvkde, ...)
