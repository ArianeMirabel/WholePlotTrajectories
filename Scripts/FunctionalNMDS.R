source("Scripts/Vernacular_handle.R")           # Script for vernacular names replacement
source("Scripts/TraitsMiceFilling.R")                                   # Traits Gap-filling

load("DB/Alpha_Plots")                  # Alpha vector for vernacular replacement
load("DB/Paracou_R_Subdivided_ok")      # Paracou dataset, lists by inventoried year, sublists by plot, etc...
dates<-sort(names(LivingStand_all))[-1] # dates considered

Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")       # Brigde functional traits dataset
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="") # Paracou15 seed and Hmax dataset
TraitsName<-c("Hmax","S_mass","L_thickness","L_chloro",
              "L_toughness","L_DryMass","SLA","WD","Bark_thick") # Considered traits list

# Treatments definition
T0<-c(1,6,11); T1<-c(2,7,9); T2<-c(3,5,10); T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("Control","T1","T2","T3")

# list of all species recorded: to have global identical NMDS classification
AllInv<-unique(unlist(lapply(LivingStand_all,function(yr){
  return(unique(unlist(lapply(yr,function(pl){return(unique(pl[,"name"]))}))))})))
AllInv<-as.data.frame(unique(AllInv),row.names=as.character(AllInv))

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
