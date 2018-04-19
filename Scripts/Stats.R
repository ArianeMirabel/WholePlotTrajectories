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
cor(Shannon[which(!Shannon$plot%in%c(8:12)),"Max"],Shannon[which(!Shannon$plot%in%c(8:12)),"treat"],method="spearman")
colnames(Simpson)<-c("plot","Max","treat")
cor(Simpson[which(!Simpson$plot%in%c(8:12)),"Max"],Simpson[which(!Simpson$plot%in%c(8:12)),"treat"],method="spearman")


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





