load("DB/LivingStand_carre")
source("Scripts/Vernacular_handle.R")

# Treatments definition
T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3);names(treatments)<-c("T0","T1","T2","T3")

# Years with identified problems of inventories
Pby<-c("1996","1998","2000","2002","2004","2006","2008","2010","2012","2014","2016","2017")
dates<-sort(names(LivingStand_all)[which(!names(LivingStand_all)%in%Pby)])

# Alpha matrix for vernacular names association
alphas_plot<-AlphaPlots(LivingStand_all)

# Reference for botany (Genus/species association)
load("DB/BotanyGenus")

# Moving average function
smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}

#trajectories drawn by treatment
CompleteTaxo<-lapply(1:4,function(t){
  
  treat<-treatments[[t]]
  Nrep<-5   # Number of iteration: replacement of vernacular names
  print(t)
  
  Trajectories<-lapply(c(0,2),function(Q){
    
          Repet<-mclapply(1:Nrep,function(r){
            
    Plot_trajectory<-do.call(cbind,lapply(dates,function(y){
      
      # Select current year and current plots
      year<-LivingStand_all[[which(names(LivingStand_all)==y)]] 
      yearPlot<-year[which(names(year)%in%treat)] 
      
      if(length(yearPlot)!=0){
        
        Div<-lapply(1:length(yearPlot),function(pl){
          
        # Replace vernacular names and get corresponding genus
        Div_carre<-unlist(lapply(unique(yearPlot[[pl]][,"n_carre"]),function(car){
          carr<-yearPlot[[pl]][which(yearPlot[[pl]][,"n_carre"]==car),]
          carr<-as.data.frame(Replacement(carr,Alpha=alphas_plot[[names(yearPlot)[pl]]]))
           colnames(carr)<-"names"
           carr<-as.character(merge(carr,RefBota,by.x="names",by.y="row.names",all.x=T)[,"Genre"])
         
        return(expq(bcTsallis(as.AbdVector(tapply(carr,carr,length)),q=Q,Correction = "None"),q=Q))
        }))
        names(Div_carre)<-unique(yearPlot[[pl]][,"n_carre"])
        return(Div_carre)
     })
        
        names(Div)<-names(yearPlot)
        return(unlist(Div))
    }
      
      if(length(yearPlot)==0){return(rep(NA,length(treat)*4))}
    }))
    
    Plot_trajectory<-Plot_trajectory[,which(!apply(Plot_trajectory,2,anyNA))]
    Plot_trajectory<-smooth(Plot_trajectory,2)
    
    return(smooth(Plot_trajectory,2)) # moving average, path=2
})
          
  Repet<-array(unlist(Repet),dim=c(length(treat)*4,length(dates),Nrep),
               dimnames=list(rownames(Repet[[1]]),dates,1:Nrep))
   
  #Calculate the 95% interval and median accross repettions
  Rnames<-rownames(Repet)
  Repet<-lapply(c(0.025,0.5,0.975),function(p){
    apply(Repet,c(1,2),function(x){return(quantile(x,probs=p,na.rm=T))})})
  Repet<-array(unlist(Repet),dim=c(nrow(Repet[[1]]),ncol(Repet[[1]]),3),
                dimnames=list(Rnames,dates,c(0.025,0.5,0.975)))
  return(Repet)
  })
  
  Trajectories<-array(unlist(Trajectories),dim=c(nrow(Trajectories[[1]]),ncol(Trajectories[[1]]),3,2),
      dimnames=list(rownames(Trajectories[[1]]),dates,c(0.025,0.5,0.975),c("Richness","Simpson")))
})

save(CompleteTaxo,file = "DB/ReplacementTraj_ForGraphs_carre")

