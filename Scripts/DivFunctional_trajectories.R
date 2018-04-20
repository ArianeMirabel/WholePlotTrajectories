library(RColorBrewer)
library(mice)
library(ade4)
library(cluster)
source("Scripts/Vernacular_handle.R")
#source("NMDS/TraitsDB.R")
source("Scripts/TraitsMiceFilling.R")
load("DB/Paracou_R_Subdivided_ok")
load("DB/Alpha_Plots")

###################################   Data

Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="")
TraitsName<-c("Hmax","L_thickness","L_chloro","L_toughness","L_DryMass","SLA","WD","Bark_thick")

InventorySp<-do.call(rbind,lapply(LivingStand_all,function(yr){
  ret<-do.call(rbind,yr)[,c("Famille","Genre","name")]
  ret<-ret[which(!grepl("Indet.",ret[,"name"])),]
  return(ret[which(!duplicated(ret)),])}))
InventorySp<-InventorySp[which(!duplicated(InventorySp)),]

#Plots of the different treatments
T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3);names(treatments)<-c("T0","T1","T2","T3")

Pby<-c("1996","1998","2000","2002","2004","2006","2008","2010","2012","2014","2016","2017")
dates<-sort(names(LivingStand_all)[which(!names(LivingStand_all)%in%Pby)])

# Moving average function
smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}


CompleteFun<-lapply(1:4,function(t){
  
  treat<-treatments[[t]]
  Nrep<-50
  
   Mplot<-mclapply(1:Nrep,function(r){
      
      VernDiv_living<-do.call(cbind,lapply(dates,function(y){
        x1<-LivingStand_all[[which(names(LivingStand_all)==y)]] 
        
        M<-x1[which(names(x1)%in%treat)] 
        
        if(length(M)!=0){ #If there are surveys for these year and treatment
          
          Mean<-lapply(1:length(M),function(x3){
            
            Mc<-Replacement(M[[x3]],Alpha=alphas_plot[[names(M)[x3]]])
            
            traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
            traits_filled<-aggregate(traits_filled[,TraitsName],list(traits_filled$name),median)
            #If seed mass
            #traits_filled<-aggregate(traits_filled[,TraitsName],list(rownames(traits_filled)),median)
            rownames(traits_filled)<-traits_filled[,1];traits_filled<-traits_filled[,TraitsName]
            tra<-traits_filled[which(rownames(traits_filled)%in%Mc),];tra<-tra[order(rownames(tra)),]
            Mc<-Mc[which(Mc%in%rownames(tra))]
            
            dissim<-as.matrix(daisy(tra,metric="gower"))
            dissim <- 1 - dissim/max(dissim)
            
            return(expq(Hqz(as.AbdVector(tapply(Mc,Mc,length)), q=2, Z=dissim,Correction="None"),q=2))
          })
      
          names(Mean)<-names(M)
          Mean<-do.call(rbind,Mean) #bind the plots vector of diversity per year
          Mean<-merge(Mean,data.frame(treat),by.x="row.names",by.y="treat",all=T)
          #To have the same number of plots each year, fill missing plots with NA
          
          return(Mean[order(Mean[,"Row.names"]),2]) 
        }
        
        if(length(M)==0){return(rep(NA,length(treat)))}
      })) #returns the binded columns by year
      
      rownames(VernDiv_living)<-treat
      colnames(VernDiv_living)<-dates
      VernDiv_living<-VernDiv_living[,which(!apply(VernDiv_living,2,anyNA))]
      
      return(smooth(VernDiv_living,2))})
  
  return(array(unlist(Mplot),dim=c(length(treat),length(dates),Nrep),
               dimnames=list(treat,dates,1:Nrep)))
})

save(CompleteFun,file="DB/FunctionalTraj_ForGraphs")

