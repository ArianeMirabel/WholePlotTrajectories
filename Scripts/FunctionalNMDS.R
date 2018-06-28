source("Scripts/Vernacular_handle.R")           # Script for vernacular names replacement
source("Scripts/TraitsMiceFilling.R")                                   # Traits Gap-filling

load("DB/Alpha_Plots")                  # Alpha vector for vernacular replacement
load("DB/Paracou_R_Subdivided_ok")      # Paracou dataset, lists by inventoried year, sublists by plot, etc...
dates<-sort(names(LivingStand_all))[-1] # dates considered
dates<-dates[which(!dates%in%c("1996","1998","2000","2002","2004","2006","2008","2010","2012","2014","2016"))]


Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")       # Brigde functional traits dataset
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="") # Paracou15 seed and Hmax dataset
TraitsName<-c("Hmax","L_thickness","L_chloro",
              "L_toughness","L_DryMass","SLA","WD","Bark_thick") # Considered traits list

# Treatments definition
T0<-c(1,6,11); T1<-c(2,7,9); T2<-c(3,5,10); T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("Control","T1","T2","T3")

# list of all species recorded: to have global identical NMDS classification
AllInv<-unique(unlist(lapply(LivingStand_all,function(yr){
  return(unique(unlist(lapply(yr,function(pl){return(unique(pl[,"name"]))}))))})))
AllInv<-as.data.frame(unique(AllInv),row.names=as.character(AllInv))

InventorySp<-do.call(rbind,lapply(LivingStand_all,function(yr){
  ret<-do.call(rbind,yr)[,c("Famille","Genre","name")]
  ret<-ret[which(!grepl("Indet.",ret[,"name"])),]
  return(ret[which(!duplicated(ret)),])}))
InventorySp<-InventorySp[which(!duplicated(InventorySp)),]

Nrep<-50 # Iterations for the vernacular names replacement and traits gap filling

MatrepFun<-lapply(1:Nrep,function(rep){
  
  Mat<-lapply(1:12,function(p){
    
    # For on plot (p),the species proba vector of successive inventories
    Trajnmds<-lapply(dates,function(y){
      
      temp<-LivingStand_all[[which(names(LivingStand_all)==y)]]
      if(!any(names(temp)==p)){return(NA)}
      if(any(names(temp)==p)){
        temp<-temp[[which(names(temp)==p)]]
        temp<-temp[!duplicated(temp),]
        temp<-Replacement(temp,Alpha=alphas_plot[[which(names(alphas_plot)==p)]]) # Verncular names replacement
        temp[is.na(temp)]<-0  
        return(as.ProbaVector(tapply(temp,temp,length)))}})
    names(Trajnmds)<-dates
    
    # Remove the tricky years (different inventory protocols): we discard years where more than 90 trees disapeared
    Ok<-names(Trajnmds)[which(unlist(lapply(Trajnmds,function(yr){return(any(!is.na(yr)))})))]
    Trajnmds<-Trajnmds[which(names(Trajnmds)%in%Ok)]
    Ntrees<-unlist(lapply(Trajnmds,length))
    Ok<-names(which(unlist(lapply(2:length(Trajnmds),function(step){return(Ntrees[step]-Ntrees[step-1])}))<=-90))
    Trajnmds<-Trajnmds[which(!names(Trajnmds)%in%Ok)]
    namesOk<-names(Trajnmds)
    
    # Homogenize the inventories, according to the general species list
    Trajnmds<-do.call(cbind,lapply(Trajnmds,function(yr){
      ret<-merge(AllInv,yr,by="row.names",all=TRUE)[,c(1,3)]
      ret[is.na(ret)]<-0;rownames(ret)<-ret[,"Row.names"]
      return(ret[2])}))
    colnames(Trajnmds)<-paste(p,"_",namesOk,sep="")
    return(Trajnmds)})
    
# Final matrix: one column per plot, per sampled years
Mat<-do.call(cbind,Mat)
    
# Traits gap filling
##### Warning: very long command!
traits_filled<-Traits_filling(Traits1,Traits2)
#####
traits_filled<-aggregate(traits_filled[,TraitsName],list(traits_filled$name),median) # Median values for species level data
rownames(traits_filled)<-traits_filled[,1];traits_filled<-traits_filled[,TraitsName]
    
tra<-traits_filled[which(rownames(traits_filled)%in%rownames(Mat)),]

# Attribute traits to none recorded species: sample traits of another same genus / family species
Corrinv<-do.call(rbind,lapply(LivingStand_all,function(yr){
  Ret<-do.call(rbind,lapply(yr,function(pl){
    ret<-pl[,c("Famille","Genre","name")];return(ret[!duplicated(ret),])}))
    return(Ret[!duplicated(Ret),])}))
Corrinv<-Corrinv[!duplicated(Corrinv),]
Corrinv<-as.data.frame(Corrinv,row.names=as.character(Corrinv[,"name"]))
missings<-setdiff(rownames(Mat),rownames(traits_filled))
TraitComp<-lapply(missings,function(sp){
  traitDB<-Corrinv[which(Corrinv[,"Genre"]==Corrinv[sp,"Genre"]),]
  if(!any(!rownames(traitDB)%in%missings)){
    traitDB<-Corrinv[which(Corrinv[,"Famille"]==Corrinv[sp,"Famille"]),]}
  if(any(!rownames(traitDB)%in%missings)){
    ret<-traits_filled[sample(setdiff(rownames(traitDB),missings),1),]
    rownames(ret)<-sp
  return(ret)
  }
  })
TraitComp<-do.call(rbind,TraitComp[which(!unlist(lapply(TraitComp,is.null)))])
TraitComp<-TraitComp[which(!apply(TraitComp,1,anyNA)),]

tra<-rbind(tra,TraitComp);tra<-tra[order(rownames(tra)),]
tra<-tra[which(!is.na(tra[,"Hmax"])),]

# Homogenize the datasets, and measure the CWM of every inventory
Mat<-Mat[which(rownames(Mat)%in%rownames(tra)),];Mat<-Mat[order(rownames(Mat)),]
Mat2<-t(apply(Mat,2,function(Inv){Inv<-colSums(apply(scale(tra),2,function(col){col<-col*Inv}))}))
# Compute distance matrix according to inventories CWM
Mat2<-as.matrix(daisy(Mat2,metric="gower"))
    
# Return the coordinates
    return(scores(metaMDS(Mat2,distance="bray"),display="sites"))
})

save(MatrepFun,file = "DB/FunctionalComposition_forGraph")


###########################################################
## traits density distance

library(sfsmisc)

traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
traits_filled<-aggregate(traits_filled[,TraitsName],list(traits_filled$name),median)
rownames(traits_filled)<-traits_filled[,1];traits_filled<-traits_filled[,TraitsName]

Nrep<-3

MatrepFun<-lapply(1:Nrep,function(rep){
  
Mat<-do.call(rbind,lapply(1:12,function(p){
    
    #traits_filled<-aggregate(traits_filled[,TraitsName],list(traits_filled$name),median) # Median values for species level data
    #rownames(traits_filled)<-traits_filled[,1];traits_filled<-traits_filled[,TraitsName]
        
densities<-lapply(TraitsName,function(Trt){
    
    ref<-LivingStand_all[[which(names(LivingStand_all)=="1985")]]
    ref<-ref[[which(names(ref)==p)]]
    ref<-Replacement(ref,Alpha=alphas_plot[[which(names(alphas_plot)==p)]])
    traRef<-traits_filled[which(rownames(traits_filled)%in%ref),Trt]
    names(traRef)<-rownames(traits_filled)[which(rownames(traits_filled)%in%ref)]
    traRef<-traRef[order(names(traRef))]
    bwRef<-bw.SJ(traRef)
    
    ref<-ref[which(ref%in%names(traRef))]
    ref<-as.ProbaVector(tapply(ref,ref,length))
    ref<-ref[order(names(ref))]
    
    dRef<-density(traRef,weights=ref,bw=bwRef)
    
     Trajnmds<-lapply(dates,function(y){
      
       temp<-LivingStand_all[[which(names(LivingStand_all)==y)]]
       if(!any(names(temp)==p)){return(NA)}
       if(any(names(temp)==p)){
         temp<-temp[[which(names(temp)==p)]]
         temp<-temp[!duplicated(temp),]
         temp<-Replacement(temp,Alpha=alphas_plot[[which(names(alphas_plot)==p)]])
    
         tra<-traits_filled[which(rownames(traits_filled)%in%temp),Trt]
         names(tra)<-rownames(traits_filled)[which(rownames(traits_filled)%in%temp)]
         tra<-tra[order(names(tra))]
         temp<-temp[which(temp%in%names(tra))]
         temp<-as.ProbaVector(tapply(temp,temp,length))
         temp<-temp[order(names(temp))]
         
         dYear<-density(tra,weights=temp,bw=bwRef)
         
         adRef <- approx(dRef$x, dRef$y, xout=dYear$x)
         adRef$y[is.na(adRef$y)] <- 0
         
         adYear <- approx(dYear$x, dYear$y, xout=dRef$x)
         adYear$y[is.na(adYear$y)] <- 0
         x <- c(dRef$x, dYear$x)
         yRef <- c(dRef$y, adRef$y)
         yYear <- c(adYear$y, dYear$y)
         sortedXindex <- sort(x, index.return=TRUE)$ix
         yRef <- yRef[sortedXindex]
         yYear <- yYear[sortedXindex]
         x <- x[sortedXindex]
         
         w <- pmin(yRef, yYear)
         total <- integrate.xy(x, yRef) + integrate.xy(x, yYear)
         intersection <- integrate.xy(x, w)
         
         return(1 -  2*intersection/total)}
      })         
      
    names(Trajnmds)<-dates
    Ok<-names(Trajnmds)[which(unlist(lapply(Trajnmds,function(yr){return(any(!is.na(yr)))})))]
    Trajnmds<-unlist(Trajnmds[which(names(Trajnmds)%in%Ok)])
    #Ntrees<-unlist(lapply(Trajnmds,length))
    #Ok<-names(which(unlist(lapply(2:length(Trajnmds),function(step){return(Ntrees[step]-Ntrees[step-1])}))<=-90))
    #Trajnmds<-Trajnmds[which(!names(Trajnmds)%in%Ok)]
    #namesOk<-names(Trajnmds)
    return(Trajnmds)
 })
  
if(!unique(unlist(lapply(densities,length)))>1){
  print(p)
  return(densities)
}
if(unique(unlist(lapply(densities,length)))>1){
  densities<-do.call(rbind,densities)
return(colSums(densities))
}

}))

rownames(Mat)<-1:12

return(Mat)
})


Matrep<-Mat


Matrep<-array(unlist(Matrep),dim=c(nrow(Matrep[[1]]),ncol(Matrep[[1]]),length(Matrep)),
                  dimnames=list(1:nrow(Matrep[[1]]),colnames(Matrep[[1]]),1:length(Matrep)))

ret<-lapply(c(0.025,0.5,0.975),function(quant){
    return(apply(Matrep,c(1,2),function(col){return(quantile(col,probs=quant))}))})
  names(ret)<-c(0.025,0.5,0.975)
  
  plot(colnames(ret[[2]]),ret[[2]][1,],type="n",xlab="",ylab="",
       ylim=c(min(unlist(ret)),max(unlist(ret))),cex.axis=0.7)
  invisible(lapply(1:length(treatments),function(tr){
    Toplot<-lapply(treatments[[tr]],function(plo){return(do.call(rbind,lapply(ret,function(quant){return(quant[plo,])})))})
    invisible(lapply(Toplot,function(plo){
      lines(colnames(plo),plo["0.5",],col=ColorsTr[[tr]],lwd=2)
      polygon(c(colnames(plo),rev(colnames(plo))),c(plo["0.025",],rev(plo["0.975",])),
              col=rgb(0,0,0,alpha=0.1),border=NA)
    }))
  }))
  mtext("Distance from 1989 inventory",side=2,padj=0,line=2,cex=0.8)
  mtext("(c)",side=3,adj=0,line=0.5)

save(FunEuclid,file="DB/FunDistance_ForGraphs")
windows()
par(mfrow=c(1,2))
EuclidDist(FunEuclid)
EuclidDist(TaxoEuclid)

