load("DB/Paracou_R_Subdivided_ok")
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
  Nrep<-50   # Number of iteration: replacement of vernacular names
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
         
        return(expq(bcTsallis(as.AbdVector(tapply(Div_plot,Div_plot,length)),q=Q,Correction = "None"),q=Q))
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
    Plot_trajectory<-Plot_trajectory[,which(!apply(Plot_trajectory,2,anyNA))]
    
    return(smooth(Plot_trajectory,2)) # moving average, path=2
})
          
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




EuclidDist<-function(Distances){
  colnames(Distances)<-as.numeric(colnames(Distances))-1986
  ret<-lapply(c(0.025,0.5,0.975),function(quant){
    return(apply(Distances,c(1,2),function(col){
      return(quantile(col,probs=quant))}))})
  names(ret)<-c(0.025,0.5,0.975)
  
  plot(colnames(ret[[2]]),ret[[2]][1,],type="n",xlab="",ylab="",xaxt="n",
       ylim=c(min(unlist(ret)),max(unlist(ret))),cex.axis=0.7)
  axis(1,at=as.numeric(colnames(ret[[1]])),labels=TRUE)  
  invisible(lapply(1:length(treatments),function(tr){
    Toplot<-lapply(treatments[[tr]],function(plo){return(do.call(rbind,lapply(ret,function(quant){return(quant[plo,])})))})
    invisible(lapply(Toplot,function(plo){
      lines(colnames(plo),plo["0.5",],col=ColorsTr[[tr]],lwd=2)
      polygon(c(colnames(plo),rev(colnames(plo))),c(plo["0.025",],rev(plo["0.975",])),
              col=rgb(0,0,0,alpha=0.1),border=NA)
    }))
  }))
}

load("DB/FunctionalComposition_forGraphs");load("DB/TaxoComposition_forGraphs")
load("DB/FunDistance_ForGraphs");load("DB/TaxoDistance_ForGraphs")

windows()
par(mfrow=c(2,2),mar=c(2,2.5,2.5,1),oma=c(1,2,1,4),no.readonly = T)
TaxoCompo(MatrepTaxo)
mtext("Taxonomic composition",side=3,line=1.8,cex=1.05)
FunCompo(MatrepFun)
mtext("Functional composition",side=3,line=1.8,cex=1.05)
EuclidDist(TaxoEuclid)
mtext("Euclidean distance\nfrom 1984 inventory",side=2,padj=0,line=2,cex=0.8)
mtext("(c)",side=3,adj=0,line=0.5)
EuclidDist(FunEuclid)
mtext("(d)",side=3,adj=0,line=0.5)
mtext("Years since disturbance",side=1,line=2.2,adj=1,cex=0.8)
legend("right",inset=-0.4,xpd=NA,legend=c("Control","Low","Interm.","High"),col=ColorsTr,lwd=2.5,bty="n",title="Treatment", cex=0.85)

