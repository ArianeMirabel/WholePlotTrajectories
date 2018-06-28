library(vegan)
library(reshape2)
library(ade4)
library(cluster)
source("Scripts/Vernacular_handle.R")
load("DB/Alpha_Plots")
load("DB/Paracou_R_Subdivided_ok")

T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3);names(treatments)<-c("Control","T1","T2","T3")

ColorsTr<-c("darkolivegreen2","deepskyblue2","darkorange1","red2")

dates<-sort(names(LivingStand_all))[-1]
dates<-dates[which(!dates%in%c("1996","1998","2000","2002","2004","2006","2008","2010","2012","2014","2016"))]

AllInv<-do.call(rbind,lapply(LivingStand_all,function(yr){
  Ret<-do.call(rbind,lapply(yr,function(pl){
    ret<-pl[,c("Genre","name")]
    return(ret[which(!duplicated(ret)),])}))
  return(Ret[which(!duplicated(Ret)),])}))
AllInv<-AllInv[which(!duplicated(AllInv)),]
AllInv<-as.data.frame(unique(AllInv),row.names=as.character(AllInv[,"name"]))

Nrep<-50

Matrep<-lapply(1:Nrep,function(rep){
  Mat<-lapply(1:12,function(p){
    Trajnmds<-lapply(dates,function(y){
      temp<-LivingStand_all[[which(names(LivingStand_all)==y)]]
      if(!any(names(temp)==p)){return(NA)}
      if(any(names(temp)==p)){
        temp<-temp[[which(names(temp)==p)]]
        temp<-temp[!duplicated(temp),]
        temp<-as.data.frame(Replacement(temp,Alpha=alphas_plot[[which(names(alphas_plot)==p)]]))
        colnames(temp)<-"name"
        temp<-as.character(merge(AllInv,temp,by.x="row.names",by.y="name",all.y=TRUE)[,"Genre"])
        return(as.ProbaVector(tapply(temp,temp,length)))}})
    names(Trajnmds)<-dates
    Ok<-names(Trajnmds)[which(unlist(lapply(Trajnmds,function(yr){return(any(!is.na(yr)))})))]
    Trajnmds<-Trajnmds[which(names(Trajnmds)%in%Ok)]
    Ntrees<-unlist(lapply(Trajnmds,length))
    Ok<-names(which(unlist(lapply(2:length(Trajnmds),function(step){return(Ntrees[step]-Ntrees[step-1])}))<=-90))
    Trajnmds<-Trajnmds[which(!names(Trajnmds)%in%Ok)]
    namesOk<-names(Trajnmds)
    ref<-as.data.frame(Trajnmds[["1985"]])
    Trajnmds<-do.call(cbind,lapply(Trajnmds,function(yr){
      All<-as.data.frame(AllInv[, "Genre"]);colnames(All)<-"Genre"
      ret<-merge(All,yr,by.x="Genre",by.y="row.names",all.x=TRUE)
      ret<-merge(ret,ref,by.x="Genre",by.y="row.names",all.x=TRUE)
      ret<-ret[which(!duplicated(ret)),2:3]
      ret[is.na(ret)]<-0
      return(dist(t(ret)))}))
    colnames(Trajnmds)<-namesOk
    return(Trajnmds)})
    
Mat<-do.call(rbind,Mat)
Mat<-t(Mat[which(rowSums(Mat)>0),])

ret<-metaMDS(Mat,distance="bray")
k<-0
while(!ret$converged & k<=10){ret<-metaMDS(Mat,distance="bray");k<-k+1;print(k)}

    return(scores(ret,display="sites"))
})

save(MatrepTaxo,file="DB/TaxoComposition_ForGraphs")

#####################################################
### Taxonomic Euclidean distance

Nrep<-10

TaxoEuclid<-lapply(1:Nrep,function(rep){
  Mat<-lapply(1:12,function(p){
    Trajnmds<-lapply(dates,function(y){
      temp<-LivingStand_all[[which(names(LivingStand_all)==y)]]
      if(!any(names(temp)==p)){return(NA)}
      if(any(names(temp)==p)){
        temp<-temp[[which(names(temp)==p)]]
        temp<-temp[!duplicated(temp),]
        temp<-as.data.frame(Replacement(temp,Alpha=alphas_plot[[which(names(alphas_plot)==p)]]))
        colnames(temp)<-"name"
        temp<-as.character(merge(AllInv,temp,by.x="row.names",by.y="name",all.y=TRUE)[,"Genre"])
        return(as.ProbaVector(tapply(temp,temp,length)))}})
    names(Trajnmds)<-dates
    Ok<-names(Trajnmds)[which(unlist(lapply(Trajnmds,function(yr){return(any(!is.na(yr)))})))]
    Trajnmds<-Trajnmds[which(names(Trajnmds)%in%Ok)]
    Ntrees<-unlist(lapply(Trajnmds,length))
    Ok<-names(which(unlist(lapply(2:length(Trajnmds),function(step){return(Ntrees[step]-Ntrees[step-1])}))<=-90))
    Trajnmds<-Trajnmds[which(!names(Trajnmds)%in%Ok)]
    namesOk<-names(Trajnmds)
    ref<-as.data.frame(Trajnmds[["1985"]])
    Trajnmds<-do.call(cbind,lapply(Trajnmds,function(yr){
      All<-as.data.frame(AllInv[, "Genre"]);colnames(All)<-"Genre"
      ret<-merge(All,yr,by.x="Genre",by.y="row.names",all.x=TRUE)
      ret<-merge(ret,ref,by.x="Genre",by.y="row.names",all.x=TRUE)
      ret<-ret[which(!duplicated(ret)),2:3]
      ret[is.na(ret)]<-0
      return(dist(t(ret)))}))
    colnames(Trajnmds)<-namesOk
    return(Trajnmds)})
  
  return(do.call(rbind,Mat))})

TaxoEuclid<-array(unlist(Matrep),dim=c(nrow(Matrep[[1]]),ncol(Matrep[[1]]),length(Matrep)),
           dimnames=list(1:nrow(Matrep[[1]]),colnames(Matrep[[1]]),1:length(Matrep)))

save(TaxoEuclid,file="DB/TaxoDistance_ForGraphs")


