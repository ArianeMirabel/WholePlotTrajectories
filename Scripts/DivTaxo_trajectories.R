load("DB/Paracou_R_Subdivided_ok")
source("Vernacular_handle.R")

# Treatments definition
T0<-c(1,6,11);T1<-c(2,7,9);T2<-c(3,5,10);T3<-c(4,8,12)
treatments<-list(T0,T1,T2,T3);names(treatments)<-c("T0","T1","T2","T3")

# Years with identified problems of inventories
Pby<-c("1996","1998","2000","2002","2004","2006","2008","2010","2012","2014")
dates<-sort(names(LivingStand_all)[which(!names(LivingStand_all)%in%Pby)])

# Alpha matrix for vernacular names association
alphas_plot<-AlphaPlots(LivingStand_all)

# Reference for botany (Genus/species association)
load("DB/BotanyGenus")

#trajectories drawn by treatment
CompleteTaxo<-lapply(1:4,function(t){
  
  treat<-treatments[[t]]
  Nrep<-2   # Number of iteration: replacement of vernacular names
  print(t)
  
  Trajectories<-lapply(0:2,function(Q){
    
          Repet<-mclapply(1:Nrep,function(r){
            
    Plot_trajectory<-do.call(cbind,lapply(dates,function(y){
      
      # Select current year and current plots
      year<-LivingStand_all[[which(names(LivingStand_all)==y)]] 
      yearPlot<-year[which(names(year)%in%treat)] 
      
      if(length(yearPlot)!=0){
        
        Div<-lapply(1:length(yearPlot),function(pl){
          
        # Replace vernacular names and get corresponding genus
        Div_plot<-as.data.frame(Replacement(yearPlot[[pl]],Alpha=alphas_plot[[names(yearPlot)[pl]]]))
        colnames(Div_plot)<-"names"
        Div_plot<-as.character(merge(Div_plot,RefBota,by.x="names",by.y="row.names",all.x=T)[,"Genre"])
         
        return(expq(bcTsallis(as.AbdVector(tapply(Div_plot,Div_plot,length)),q=Q,Correction = "Best"),q=Q))
       })
        
        names(Div)<-names(yearPlot)
        Div<-do.call(rbind,Div)
        Div<-merge(Div,data.frame(treat),by.x="row.names",by.y="treat",all=T)
        
        return(Div[order(Div[,"Row.names"]),2]) 
      }
      
      if(length(yearPlot)==0){return(rep(NA,length(treat)))}
    }))
    
    rownames(Plot_trajectory)<-treat
    colnames(Plot_trajectory)<-dates
    
    return(Plot_trajectory)})
          
  Repet<-array(unlist(Repet),dim=c(length(treat),length(dates),Nrep),dimnames=list(treat,dates,1:Nrep))
   
  #Calculate the 95% interval and median accross repettions
  Repet<-lapply(c(0.025,0.5,0.975),function(p){
    apply(Repet,c(1,2),function(x){return(quantile(x,probs=p,na.rm=T))})})
  Repet<-array(unlist(Repet),dim=c(nrow(Repet[[1]]),ncol(Repet[[1]]),3),
                dimnames=list(treat,dates,c(0.025,0.5,0.975)))
  return(Repet)
  })
  
  Trajectories<-array(unlist(Trajectories),dim=c(nrow(Trajectories[[1]]),ncol(Trajectories[[1]]),3,3),
      dimnames=list(treat,dates,c(0.025,0.5,0.975),c("Richness","Shannon","Simpson")))
})

save(CompleteTaxo,file = "DB/ReplacementTraj_ForGraphs")
