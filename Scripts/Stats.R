load("DB/TaxoComposition_ForGraphs")

Data_TaxoComp<-MatrepTaxo

MatrepT<-lapply(Data_TaxoComp,function(Rep){
  Rep<-as.data.frame(Rep)
  Rep$plot<-substr(rownames(Rep),start=1,stop=regexpr("_",rownames(Rep))-1)
  Rep$year<-substr(rownames(Rep),start=regexpr("_",rownames(Rep))+1,
                   stop=nchar(rownames(Rep)))
  return(Rep)
})

DistT<-as.data.frame(unlist(lapply(1:12,function(pl){
  ret<-lapply(MatrepT,function(rep){return(rep[which(rep[,"plot"]==pl),])})
  ret<-lapply(ret,function(rep){return(rep[which(rep[,"year"]%in%
                              ret[[which(unlist(lapply(ret,nrow))==min(unlist(lapply(ret,nrow))))[1]]][,"year"]),])})
  ret<-do.call(rbind,lapply(ret,function(rep){
    ret2<-apply(rep[,c("NMDS1","NMDS2")],1,function(li){
      return(sqrt(sum((rep[1,c("NMDS1","NMDS2")]-li)^2)))})
    names(ret2)<-rep[,"year"]
    return(ret2)}))
  ret<-apply(ret,2,median)
  return(names(ret)[which(ret==max(ret))])})))
treats<-cbind(c(1,6,11,2,7,9,3,5,10,4,8,12),rep(0:3,each=3))
DistT<-cbind(DistT,rownames(DistT),treats[order(treats[,1]),2])
colnames(DistT)<-c("Max","Plot","treat")
cor(DistT[,"Max"],DistT[,"treat"],method="spearman")

TimeMax<-unlist(lapply(1:12,function(pl){
  ret<-lapply(MatrepT,function(rep){return(rep[which(rep[,"plot"]==pl),])})
  ret<-lapply(ret,function(rep){return(rep[which(rep[,"year"]%in%
                                                   ret[[which(unlist(lapply(ret,nrow))==min(unlist(lapply(ret,nrow))))[1]]][,"year"]),])})
  ret<-do.call(rbind,lapply(ret,function(rep){
    ret2<-apply(rep[,c("NMDS1","NMDS2")],1,function(li){
      return(sqrt(sum((rep[1,c("NMDS1","NMDS2")]-li)^2)))})
    names(ret2)<-rep[,"year"]
    return(ret2)}))
  ret<-apply(ret,2,median)
  return(names(ret)[which(ret==max(ret))])}))
mean(as.numeric(TimeMax))-1984


##############################################################

load("DB/ReplacementTraj_ForGraphs")

CompTaxo<-CompleteTaxo

treats<-cbind(c(1,6,11,2,7,9,3,5,10,4,8,12),rep(0:3,each=3))
treats<-treats[order(treats[,1]),]
colnames(treats)<-c("plot","treat"); rownames(treats)<-treats[,"plot"]

for(q in 1:3){     
  Toplot<-lapply(CompTaxo,function(tr){return(tr[,,,q])})
  Toplot<-lapply(Toplot,function(toplot){return(toplot[,which(colnames(toplot)>=1989),])})
  Toplot<-lapply(Toplot,function(tr){
    ret<-lapply(1:dim(tr)[3],function(rep){return(apply(tr[,,rep],2,function(col){col<-col-tr[,1,rep]}))})#
    ret<-array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),length(ret)),
               dimnames=list(rownames(ret[[1]]),as.numeric(colnames(ret[[1]]))-1984,1:length(ret)))
    return(apply(apply(ret,c(1,2),median),1,max))
    })
  Ret<-as.data.frame(unlist(Toplot))
  Ret$plot<-rownames(Ret)
  Ret<-Ret[order(as.numeric(Ret$plot)),]
  Ret<-merge(Ret,treats,by="plot")
  assign(c("Richness","Shannon","Simpson")[q],Ret)
}

colnames(Richness)<-c("plot","Max","treat")
cor(Richness[,"Max"],Richness[,"treat"],method="spearman")
colnames(Shannon)<-c("plot","Max","treat")
cor(Shannon[which(!Shannon$plot%in%c(8,12)),"Max"],Shannon[which(!Shannon$plot%in%c(8,12)),"treat"],method="spearman")
colnames(Simpson)<-c("plot","Max","treat")
cor(Simpson[which(!Simpson$plot%in%c(8,12)),"Max"],Simpson[which(!Simpson$plot%in%c(8,12)),"treat"],method="spearman")


######
# maximums of diversity indices

load("DB/ReplacementTraj_ForGraphs")
CompTaxo<-CompleteTaxo

q<-3
Toplot<-lapply(CompTaxo,function(tr){return(tr[,,,q])})
Toplot<-lapply(Toplot,function(toplot){return(toplot[,which(colnames(toplot)>=1989),])})
Toplot<-lapply(Toplot,function(tr){
  ret<-lapply(1:dim(tr)[3],function(rep){return(apply(tr[,,rep],2,function(col){col<-col-tr[,1,rep]}))})#
  ret<-array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),length(ret)),
             dimnames=list(rownames(ret[[1]]),as.numeric(colnames(ret[[1]]))-1984,1:length(ret)))
  ret<-lapply(c(0.025,0.5,0.975),function(quant){return(apply(ret,c(1,2),function(x){return(quantile(x,probs=quant))}))})
  return(array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),3),
               dimnames=list(rownames(ret[[1]]),colnames(ret[[1]]),c(0.025,0.5,0.975))))})


MaxRich<-apply(do.call(rbind,lapply(Toplot,function(tr){return(tr[,,"0.5"])})),1,max)
MaxRich[order(MaxRich,decreasing = T)][1:2]


######################r
### The plot 7 funcitonal diversity

tab<-do.call(rbind,lapply(1:12,function(pl){
  P<-lapply(LivingStand_all,function(yr){if(any(names(yr)==pl)){return(yr[[which(names(yr)==pl)]])}})
P<-P[which(!names(P)%in%Pby)]
#lapply(P,nrow)

Rao<-do.call(rbind,lapply(P,function(yr){
  
  yr<-Replacement(yr,Alpha=alphas_plot[[which(names(alphas_plot)==pl)]])
      
      #traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
      #Traits_filled<-aggregate(Traits_filled[,TraitsName],list(Traits_filled$name),median)
      #rownames(Traits_filled)<-Traits_filled[,1];Traits_filled<-Traits_filled[,TraitsName]
      tra<-Traits_filled[which(rownames(Traits_filled)%in%yr),];tra<-tra[order(rownames(tra)),]
      
      return(apply(tra,2,sd))
      
      yr<-yr[which(yr%in%rownames(tra))]
      #return(length(which(!yr%in%rownames(tra))))
      
      #dissim<-as.matrix(daisy(tra,metric="gower"))
      #dissim <- 1 - dissim/max(dissim)
      
      #return(expq(Hqz(as.AbdVector(tapply(yr,yr,length)), q=2, Z=dissim,Correction="None"),q=2))
    }))
Rao<-apply(Rao,2,mean)
}))
rownames(tab)<-1:12

PCA<-dudi.pca(tab,scale=T,scan=F)

windows()
plot(PCA$li[,"Axis1"],PCA$li[,"Axis2"],type="n")
text(PCA$li[,"Axis1"],PCA$li[,"Axis2"],labels=rownames(PCA$li))


P<-lapply(LivingStand_all,function(yr){if(any(names(yr)==pl)){return(yr[[which(names(yr)==pl)]])}})
P<-P[which(!names(P)%in%Pby)]
#lapply(P,nrow)

pl<-7
P<-lapply(LivingStand_all,function(yr){if(any(names(yr)==pl)){return(yr[[which(names(yr)==pl)]])}})
P<-P[which(!names(P)%in%Pby)]
Rao<-do.call(rbind,lapply(P,function(yr){
  
  yr<-Replacement(yr,Alpha=alphas_plot[[which(names(alphas_plot)==pl)]])
  
  #traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
  #Traits_filled<-aggregate(Traits_filled[,TraitsName],list(Traits_filled$name),median)
  #rownames(Traits_filled)<-Traits_filled[,1];Traits_filled<-Traits_filled[,TraitsName]
  tra<-Traits_filled[which(rownames(Traits_filled)%in%yr),];tra<-tra[order(rownames(tra)),]
  
  return(apply(tra,2,sd))
  
  #yr<-yr[which(yr%in%rownames(tra))]
  #return(length(which(!yr%in%rownames(tra))))
  
  #dissim<-as.matrix(daisy(tra,metric="gower"))
  #dissim <- 1 - dissim/max(dissim)
  
  #return(expq(Hqz(as.AbdVector(tapply(yr,yr,length)), q=2, Z=dissim,Correction="None"),q=2))
}))

