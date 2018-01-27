#### Loadings
library(vegan)
library(reshape2)
library(ade4)
library(shape)
library(cluster)
library(mvtnorm)
library(entropart)
library(lattice)
library(plotrix)
library(gtools)
library(RcppParallel)
library(parallel)

# Obtain the vernacular/botanical observed association frequency
alpha_construct<-function(Plot){
  mat_eff<-as.data.frame.matrix(with(Plot, xtabs(~name+nomPilote)))
  mat_freq <- apply(mat_eff,2,function(v){v/sum(v)})
  return(mat_freq)
}

# Built the plot-specific association matrix
AlphaPlots<-function(Inv){
  unlisted<-do.call(rbind,lapply(Inv,function(x){do.call(rbind,x)}))
  ListedByPlot<-lapply(unique(unlisted[,"n_parcelle"]),function(y){
    p<-unlisted[which(unlisted[,"n_parcelle"]==y),-which(colnames(unlisted)=="campagne")]
    p<-p[which(!duplicated(p)),]
    ap<-alpha_construct(p); return(ap)
  })
  names(ListedByPlot)<-unique(unlisted[,"n_parcelle"])
  return(ListedByPlot)
}

# Listing all vernacular/genus and family associations 
Correspondences<-function(Plot){
  Corr<-Plot[order(Plot[,"Famille"],Plot[,"Genre"],Plot[,"Espece"]),]
  Corr<-Corr[which(!duplicated(Corr)),]
  return(Corr)
}

# Multinomial/dirichlet process, return a botanical name trialed in proba vector
Dirichlet_draw<-function(V){
  Vdir<-rdirichlet(1,V)
  names(Vdir)<-names(V)
  res <- rmultinom(1,1,Vdir)
  res<-rownames(res)[which(res>0)]
  return(res)
}

####### The function to randomly draw a botanical name replacing
#the undertermined species. Use alpha matrix and botanic table
tirage <- function(alpha,bota,name,eps=0.05){ 
  
  if(name["nomPilote"]=="-"){
    trial<-sample(rownames(alpha),1)
    while(grepl("Indet.",trial)){trial<-sample(rownames(alpha),1)}
  }
  
  alpha_adjust<-alpha[,name["nomPilote"]]
  
  if(name["nomPilote"]!="-" & name["Genre"]!="Indet."){
    mismatch<-bota[intersect(which(bota[,"nomPilote"]==name["nomPilote"]),which(!bota[,"Genre"]%in%name["Genre"])),"name"]
    alpha_adjust[as.character(mismatch)]<-0
    match<-names(alpha_adjust[which(alpha_adjust!=0)])
    
    alpha_adjust[which(alpha_adjust==0)]<-eps/length(which(alpha_adjust==0))
    
    trial<-Dirichlet_draw(alpha_adjust)
    
    if(any(grep("Indet.",match,invert=T)!=0)){
      while(grepl("Indet.",trial)){trial<-Dirichlet_draw(alpha_adjust)
      }
    }
  }
  
  # Restrict the sampling pool to account for all taxonomic information
  
  if(name["nomPilote"]!="-" & name["Genre"]=="Indet." & name["Famille"]!="Indet."){
    mismatch<-bota[intersect(which(bota[,"nomPilote"]==name["nomPilote"]),which(!bota[,"Famille"]%in%name["Famille"])),"name"]
    alpha_adjust[as.character(mismatch)]<-0
    match<-names(alpha_adjust[which(alpha_adjust!=0)])
    alpha_adjust[which(alpha_adjust==0)]<-eps/length(which(alpha_adjust==0))
    
    trial<-Dirichlet_draw(alpha_adjust)
    if(any(grep("Indet.",match,invert=T))!=0){ 
      while(grepl("Indet.",trial)){trial<-Dirichlet_draw(alpha_adjust)
      } 
    }
  }
  if(name["nomPilote"]!="-" & name["Genre"]=="Indet." & name["Famille"]=="Indet."){
    alpha_adjust[which(alpha_adjust==0)]<-eps/length(which(alpha_adjust==0))
    
    trial<-Dirichlet_draw(alpha_adjust)
    while(grepl("Indet.",trial)){trial<-Dirichlet_draw(alpha_adjust)
    }
  }
  return(trial)
}


#Whole replacement process
Replacement<-function(Plot,Alpha){
  Determ<-Plot[which(!Plot["Espece"]=="Indet."),]
  Indet<-Plot[which(Plot["Espece"]=="Indet."),]
  
  Vern<-apply(Indet,2,as.character)
  Bota<-Correspondences(Plot)
  if(nrow(Indet)==0){
    Simu<-as.character(Determ[,"name"])
  }
  if(nrow(Indet)==1){
    Simu<-c(as.character(Determ[,"name"]),tirage(alpha=Alpha,bota=Bota,name=Vern))
  }
  if(nrow(Indet)>1){
    Simu<-unlist(lapply(1:nrow(Vern),function(i){tirage(alpha=Alpha,bota=Bota,name=Vern[i,])}))
    Simu<-c(as.character(Determ[,"name"]),Simu)
  }
  return(Simu)
}
