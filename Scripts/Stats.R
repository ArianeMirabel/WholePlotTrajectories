load("DB/ReplacementTraj_ForGraphs")
load("DB/FunctionalTraj_ForGraphs")
load("DB/FunctionalRichnessTraj_ForGraphs")

SimpsonTaxo<-lapply(CompleteTaxo,function(tr){return(tr[,,,"Simpson"])})
save(SimpsonTaxo,file="DB/SimpsonDraw")

Sim<-lapply(CompTaxo,function(tr){return(tr[,c("1995","2005","2015"),])})
Sim<-lapply(Sim, function(tr){
  Ret<-lapply(c(0.025,0.5,0.975),function(quant){
  return(apply(tr,c(1,2),function(x){return(quantile(x,probs=quant))}))})
  Ret<-array(unlist(Ret),dim=c(nrow(Ret[[1]]),ncol(Ret[[1]]),3),
             dimnames=list(rownames(Ret[[1]]),colnames(Ret[[1]]),c(0.025,0.5,0.975)))})
save(Sim,file="DB/SimpsonIDH")

RichTaxo<-lapply(CompleteTaxo,function(tr){return(tr[,,,"Richness"])})
Rich<-lapply(RichTaxo,function(tr){return(tr[,c("1995","2005","2015"),])})
Rich<-lapply(Rich, function(tr){
  Ret<-lapply(c(0.025,0.5,0.975),function(quant){
    return(apply(tr,c(1,2),function(x){return(quantile(x,probs=quant))}))})
  Ret<-array(unlist(Ret),dim=c(nrow(Ret[[1]]),ncol(Ret[[1]]),3),
             dimnames=list(rownames(Ret[[1]]),colnames(Ret[[1]]),c(0.025,0.5,0.975)))})
save(Rich,file="DB/RichnessIDH")

CompFun<-CompleteRichFun
CompFun[[2]]<-CompFun[[2]][which(rownames(CompFun[[2]])!=7),,]
FunRich<-lapply(CompFun, function(tr){
  Ret<-lapply(c(0.025,0.5,0.975),function(quant){
    return(apply(tr[,c("1995","2005","2015"),],c(1,2),function(x){return(quantile(x,probs=quant))}))})
  Ret<-array(unlist(Ret),dim=c(nrow(Ret[[1]]),ncol(Ret[[1]]),3),
             dimnames=list(rownames(Ret[[1]]),colnames(Ret[[1]]),c(0.025,0.5,0.975)))})
save(FunRich,file="DB/FunRichIDH")

windows()
par(mfrow=c(2,2),mar=c(3,3,2,1),oma=c(2,1,2,3),no.readonly = T)
plotDiv(SimpsonTaxo)
mtext("Simpson diversity",side=3,adj=0,line=1)
plotDiv(CompleteFun,remove=TRUE)
mtext("Rao diversity",side=3,adj=0,line=1)
mtext("Years since disturbance",side=1,line=2.2,adj=1)
legend("right",inset=c(-0.28,0),xpd=NA,legend=c("T0","T1","T2","T3"),col=ColorsTr,lwd=2.5,bty="n",title="Treatment")

plotIDH(Sim)
plotIDH(FunRich)
mtext("initial %AGB lost",side=1,line=2.2,adj=1)
legend("right",inset=c(-0.28,0),xpd=NA,legend=c("10","20","30"),col=colyear,lwd=2.5,bty="n",title="Years")


load("DB/SimpsonIDH");load("DB/RaoIDH");load("DB/LostAGB");load("DB/RichnessIDH")
time<-c("1995","2005","2015")
colyear<-c("darkgoldenrod1","darkorange2","darkred")
Data<-FunRich
AgbLoss<-AGBloss

plotIDH<-function(Data,AgbLoss){
  Ylim<-c(min(unlist(lapply(Data, function(tr){return(tr[,,"0.5"])}))),
          max(unlist(lapply(Data, function(tr){return(tr[,,"0.5"])}))))
plot(AgbLoss[,"AGB"],AgbLoss[,"AGB"],type="n",xlab="",ylab="",
       ylim=Ylim)
leg<-lapply(1:3,function(Ti){
    toplot<-unlist(lapply(Data, function(tr){return(tr[,time[Ti],"0.5"])}))
    abs<-AgbLoss[which(AgbLoss[,"plot"]%in%names(toplot)),"AGB"]
    points(abs,toplot,col=colyear[Ti],pch=20)
    Lm<-lm(toplot~abs)
    Lm2<-lm(toplot~abs+I(abs^2))
    #summary(Lm)
    #anova(Lm)
    #par(mfrow=c(2,2));plot(Lm)
    if(abs(round(summary(Lm)$adj.r.squared,2))<abs(round(summary(Lm2)$adj.r.squared,2))){
    abs_pred<-seq(min(abs),max(abs),length.out=100)
    lines(sort(abs_pred),predict(Lm2,newdata=data.frame(abs=abs_pred)),col=colyear[Ti],lwd=2.5)
    return(round(summary(Lm2)$adj.r.squared,2)) }
    if(abs(round(summary(Lm)$adj.r.squared,2))>abs(round(summary(Lm2)$adj.r.squared,2))){
    abline(a=Lm$coefficients[1],b=Lm$coefficients[2],col=colyear[Ti],lwd=2.5)
    return(round(summary(Lm)$adj.r.squared,2)) }
    
    #legend(round(summary(Lm)$adj.r.squared,2),x=min(AgbLoss),y=seq(Ylim[2],Ylim[1],length.out=10)[Ti],
    #       bty="n",lty=1,col=colyear[Ti],lwd=2.5,cex=0.8)
  })
  legend("topleft",legend=unlist(leg),bty="n",lty=1,col=colyear,lwd=2.5,cex=0.8,title=expression(paste('R'^2,'adjusted')))
}
  
plotIDH(FunRich,AGBloss)

##########################################
## Rho Spearman, Taxo
load("DB/FunDistance_ForGraphs");load("DB/TaxoDistance_ForGraphs");load("DB/LostAGB")

Dist<-TaxoEuclid
ret<-apply(apply(Dist,c(1,2),median),1,function(li){return(max(abs(li)))})
ret<-data.frame(cbind(ret,names(ret)))
colnames(ret)<-c("Dist","plot")
ret<-merge(ret,AGBloss,by="plot")
ret<-apply(ret,2,as.numeric)
cor(ret[,"Dist"],ret[,"AGB"],method="spearman")

#(old)
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

##########################################
## Rho Spearman, redundancy

smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}
treatments<-list(c(1,6,11),c(2,7,9),c(3,5,10),c(4,8,12))
names(treatments)<-c("Control","T1","T2","T3")
ColorsTr<-c("darkolivegreen2","deepskyblue2","darkorange1","red2")
colyear<-c("darkgoldenrod1","darkorange2","darkred")
time<-c("1995","2005","2015")

load("DB/Redundancy_restricted")
RedundancyTraj_restricted->Red

Maxred<-unlist(lapply(1:12,function(pl){
  ret<-do.call(rbind,lapply(Red,function(iter){return(iter[pl,])}))
  return(max(apply(ret,2,median)))
}))

treats<-cbind(c(1,6,11,2,7,9,3,5,10,4,8,12),rep(0:3,each=3))
Maxred<-cbind(Maxred,1:12)
colnames(Maxred)<-c("Max","plot")
Maxred<-merge(Maxred,AGBloss,by="plot")
cor(Maxred[,"Max"],Maxred[,"AGB"],method="spearman")


###############
## CWM distance

load("DB/CWM")

Max<-do.call(rbind,lapply(CWM,function(pl){ return(apply(pl[,,"0.5"],2,max))}))
treats<-cbind(c(1,6,11,2,7,9,3,5,10,4,8,12),rep(0:3,each=3))
Max<-apply(cbind(Max,rownames(Max)),2,as.numeric)#,treats[order(treats[,1]),2])
colnames(Max)<-c(colnames(CWM[[1]]),"plot")
Max<-merge(Max,AGBloss,by="plot")
apply(Max[,colnames(CWM[[1]])],2,function(x){return(cor(x,Max[,"AGB"],method="spearman"))})


load("DB/ReplacementTraj_ForGraphs")

CompTaxo<-CompleteTaxo

q<-2    
  Toplot<-lapply(CompTaxo,function(tr){return(tr[,,,q])})
  Toplot<-lapply(Toplot,function(toplot){return(toplot[,which(colnames(toplot)>=1989),])})
  Toplot<-lapply(Toplot,function(tr){
    ret<-lapply(1:dim(tr)[3],function(rep){return(apply(tr[,,rep],2,function(col){col<-col-tr[,1,rep]}))})#
    ret<-array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),length(ret)),
               dimnames=list(rownames(ret[[1]]),as.numeric(colnames(ret[[1]]))-1984,1:length(ret)))
    ret<-apply(ret,c(1,2),median)
    return(apply(ret,1,function(li){return(as.numeric(names(li)[which(li==max(li))[1]]))}))})
  
lapply(Toplot,mean)

##############################################################

load("DB/ReplacementTraj_ForGraphs")
load("DB/LostAGB")

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
               dimnames=list(rownames(ret[[1]]),as.numeric(colnames(ret[[1]]))-1987,1:length(ret)))
    return(apply(apply(ret,c(1,2),median),1,function(li){return(max(abs(li)))}))
    })
  Ret<-as.data.frame(unlist(Toplot))
  Ret$plot<-rownames(Ret)
  #Ret<-Ret[order(as.numeric(Ret$plot)),]
  Ret<-merge(Ret,AGBloss,by="plot")
  assign(c("Richness","Shannon","Simpson")[q],Ret)
}

colnames(Richness)<-c("plot","Max","AGB")
cor(Richness[,"Max"],Richness[,"AGB"],method="spearman")
colnames(Shannon)<-c("plot","Max","treat")
cor(Shannon[which(!Shannon$plot%in%c(8,12)),"Max"],Shannon[which(!Shannon$plot%in%c(8,12)),"treat"],method="spearman")
colnames(Simpson)<-c("plot","Max","AGB")
cor(Simpson[,"Max"],Simpson[,"AGB"],method="spearman")

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
### The plot 7 functional diversity

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

###### Functional spearman
load("DB/FunctionalRichnessTraj_ForGraphs")

CompFun<-CompleteFun

#if(remove){CompFun[[2]]<-CompFun[[2]][which(rownames(CompFun[[2]])!=7),,]}

CompFun<-do.call(c,lapply(CompFun,function(tr){
  ret<-tr[,which(colnames(tr)>=1989),]
  ret<-apply(ret,c(1,2),median)
  return(apply(ret,1,max))}))

CompFun<-cbind(CompFun,names(CompFun),rep(0:3,each=3))
colnames(CompFun)<-c("Max","Plot","treat")
cor(as.numeric(CompFun[,"Max"]),as.numeric(CompFun[,"treat"]),method="spearman")


### traits dataset coverage:

splist<-lapply(LivingStand_all,function(yr){
  return(unlist(lapply(yr,function(plo){
    freq<-as.character(plo[which(plo[,"name"]%in%traits[,"name"]),"name"])
    freq<-sum(tapply(freq,freq,length))
    return(freq/sum(tapply(as.character(plo[,"name"]),as.character(plo[,"name"]),length)))
  })))
})
splist<-lapply(splist,function(yr){return(as.numeric(yr[which(names(yr)%in%1:12)]))})
splist<-splist[which(unlist(lapply(splist,function(yr){length(yr)!=0})))]

Pby<-c("1996","1998","2000","2002","2004","2006","2008","2010","2012","2014","2016","2017")
splist<-unlist(splist[which(!names(splist)%in%Pby)])
min(splist)
mean(splist)
