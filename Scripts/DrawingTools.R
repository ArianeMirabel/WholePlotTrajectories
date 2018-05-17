### Manuscript drawing functions
Library(c("shape","plotrix"))

treatments<-list(c(1,6,11),c(2,7,9),c(3,5,10),c(4,8,12))
names(treatments)<-c("Control","T1","T2","T3")
ColorsTr<-c("darkolivegreen2","deepskyblue2","darkorange1","red2")
colyear<-c("darkgoldenrod1","darkorange2","darkred")
time<-c("1995","2005","2015")

smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}

TaxoCompo<-function(Data_TaxoComp){
par(mar=c(3,3,2,1))
matT<-as.data.frame(Data_TaxoComp[[1]])
matT$plot<-substr(rownames(matT),start=1,stop=regexpr("_",rownames(matT))-1)
matT$year<-substr(rownames(matT),start=regexpr("_",rownames(matT))+1,
                  stop=nchar(rownames(matT)))
colnames(matT)<-c("X","Y","plot","year")
plot(matT$X,matT$Y,type="n",xlab="",ylab="",cex.axis=0.7)
invisible(lapply(1:length(treatments),function(tr){
  plots<-treatments[[tr]]
  invisible(lapply(plots,function(pl){
    toplot<-matT[which(matT$plot==pl),]
    points(toplot[1,"X"],toplot[1,"Y"],col=ColorsTr[[tr]],pch=16,cex=1.7)
    points(toplot[nrow(toplot),"X"],toplot[nrow(toplot),"Y"],col=ColorsTr[[tr]],pch=1,cex=1.7,lwd=2)
    
    lines(toplot$X,toplot$Y,col=ColorsTr[[tr]],lty=2)
    invisible(lapply(round(seq(2,nrow(toplot),length.out=5)),function(step){
      Arrows(x0=toplot[step,"X"],y0=toplot[step,"Y"],
             x1=toplot[step-1,"X"],y1=toplot[step-1,"Y"],
             col=ColorsTr[tr],code=1,lwd=2,arr.length=0.2,arr.type="simple")}))
  }))
}))

mtext("NMDS 1",side=1,line=2,cex=0.8)
mtext("NMDS 2",side=2,padj=0,line=2,cex=0.8)
mtext("(a) Taxonomic composition",side=3,adj=0,line=0.5)
}

FunCompo<-function(Data_FunComp){mat<-as.data.frame(Data_FunComp[[3]])
mat$plot<-substr(rownames(mat),start=1,stop=regexpr("_",rownames(mat))-1)
mat$year<-substr(rownames(mat),start=regexpr("_",rownames(mat))+1,
                 stop=nchar(rownames(mat)))
colnames(mat)<-c("X","Y","plot","year")
plot(mat$X,mat$Y,type="n",xlab="",ylab="",cex.axis=0.7)
invisible(lapply(1:length(treatments),function(tr){
  plots<-treatments[[tr]]
  invisible(lapply(plots,function(pl){
    toplot<-mat[which(mat$plot==pl),]
    points(toplot[1,"X"],toplot[1,"Y"],col=ColorsTr[[tr]],pch=16,cex=1.7)
    points(toplot[nrow(toplot),"X"],toplot[nrow(toplot),"Y"],col=ColorsTr[[tr]],pch=1,cex=1.7,lwd=2)
    
    lines(toplot$X,toplot$Y,col=ColorsTr[[tr]],lty=2)
    invisible(lapply(round(seq(2,nrow(toplot),length.out=5)),function(step){
      Arrows(x0=toplot[step,"X"],y0=toplot[step,"Y"],
             x1=toplot[step-1,"X"],y1=toplot[step-1,"Y"],
             col=ColorsTr[tr],code=1,lwd=2,arr.length=0.2,arr.type="simple")}))
  }))
}))
mtext("(b) Functional composition",side=3,adj=0,line=0.5)
mtext("NMDS 1",side=1,line=2,cex=0.8)
mtext("NMDS 2",side=2,padj=0,line=2,cex=0.8)}

TaxoDist<-function(Data_TaxoComp){
MatrepT<-lapply(Data_TaxoComp,function(Rep){
  Rep<-as.data.frame(Rep)
  Rep$plot<-substr(rownames(Rep),start=1,stop=regexpr("_",rownames(Rep))-1)
  Rep$year<-substr(rownames(Rep),start=regexpr("_",rownames(Rep))+1,
                   stop=nchar(rownames(Rep)))
  return(Rep)
})

DistT<-lapply(1:12,function(pl){
  ret<-lapply(MatrepT,function(rep){return(rep[which(rep[,"plot"]==pl),])})
  ret<-lapply(ret,function(rep){return(rep[which(rep[,"year"]%in%
                                                   ret[[which(unlist(lapply(ret,nrow))==min(unlist(lapply(ret,nrow))))[1]]][,"year"]),])})
  ret<-do.call(rbind,lapply(ret,function(rep){
    ret2<-apply(rep[,c("NMDS1","NMDS2")],1,function(li){
      return(sqrt(sum((rep[1,c("NMDS1","NMDS2")]-li)^2)))})
    names(ret2)<-rep[,"year"]
    return(ret2)}))
  ret<-do.call(rbind,lapply(c(0.025,0.5,0.975),function(quant){
    return(apply(ret,2,function(col){return(quantile(col,probs=quant))}))}))
  rownames(ret)<-c(0.025,0.5,0.975)
  return(ret)})
names(DistT)<-1:12
DistT<-lapply(DistT,function(pl){
  colnames(pl)<-as.numeric(colnames(pl))-1984;return(pl)})

plot(colnames(DistT[[1]]),DistT[[1]][1,],type="n",xlab="",ylab="",
     ylim=c(min(unlist(DistT)),max(unlist(DistT))),cex.axis=0.7)
invisible(lapply(1:length(treatments),function(tr){
  toplot<-DistT[treatments[[tr]]]
  invisible(lapply(toplot,function(pl){
    lines(colnames(pl),pl["0.5",],col=ColorsTr[[tr]],lwd=2)
    polygon(c(colnames(pl),rev(colnames(pl))),c(pl["0.025",],rev(pl["0.975",])),
            col=rgb(0,0,0,alpha=0.1),border=NA)
  }))}))
mtext("Distance from 1989 inventory",side=2,padj=0,line=2,cex=0.8)
mtext("(c)",side=3,adj=0,line=0.5)
}

FunDist<-function(Data_FunComp){
  MatrepF<-lapply(Data_FunComp,function(Rep){
    Rep<-as.data.frame(Rep)
    Rep$plot<-substr(rownames(Rep),start=1,stop=regexpr("_",rownames(Rep))-1)
    Rep$year<-substr(rownames(Rep),start=regexpr("_",rownames(Rep))+1,
                     stop=nchar(rownames(Rep)))
    return(Rep)
  })
  
  Dist<-lapply(1:12,function(pl){
    ret<-lapply(MatrepF,function(rep){return(rep[which(rep[,"plot"]==pl),])})
    ret<-lapply(ret,function(rep){return(rep[which(rep[,"year"]%in%ret[[1]][,"year"]),])})
    ret<-do.call(rbind,lapply(ret,function(rep){
      ret2<-apply(rep[,c("NMDS1","NMDS2")],1,function(li){
        return(sqrt(sum((rep[1,c("NMDS1","NMDS2")]-li)^2)))})
      names(ret2)<-rep[,"year"]
      return(ret2)}))
    ret<-do.call(rbind,lapply(c(0.025,0.5,0.975),function(quant){
      return(apply(ret,2,function(col){return(quantile(col,probs=quant))}))}))
    rownames(ret)<-c(0.025,0.5,0.975)
    return(ret)})
  names(Dist)<-1:12
  Dist<-lapply(Dist,function(pl){
    colnames(pl)<-as.numeric(colnames(pl))-1984;return(pl)})
  
  plot(colnames(Dist[[1]]),Dist[[1]][1,],type="n",xlab="",ylab="",
       ylim=c(min(unlist(Dist)),max(unlist(Dist))),cex.axis=0.7)
  invisible(lapply(1:length(treatments),function(tr){
    toplot<-Dist[treatments[[tr]]]
    invisible(lapply(toplot,function(pl){
      lines(colnames(pl),pl["0.5",],col=ColorsTr[[tr]],lwd=2)
      polygon(c(colnames(pl),rev(colnames(pl))),c(pl["0.025",],rev(pl["0.975",])),
              col=rgb(0,0,0,alpha=0.1),border=NA)
    }))}))
  mtext("(d)",side=3,adj=0,line=0.5)
  mtext("Years since disturbance",side=1,adj=1,line=2,cex=0.8)
}

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
    abline(a=Lm$coefficients[1],b=Lm$coefficients[2],col=colyear[Ti],lwd=2.5)
    return(round(summary(Lm)$adj.r.squared,2))
   #legend(round(summary(Lm)$adj.r.squared,2),x=min(AgbLoss),y=seq(Ylim[2],Ylim[1],length.out=10)[Ti],
    #       bty="n",lty=1,col=colyear[Ti],lwd=2.5,cex=0.8)
  })
  legend("topleft",legend=unlist(leg),bty="n",lty=1,col=colyear,lwd=2.5,cex=0.8,title=expression(paste('R'^2,'adjusted')))
}

plotDiv<-function(Data,remove=FALSE){
  
  if(remove){Data[[2]]<-Data[[2]][which(rownames(Data[[2]])!=7),,]}
  
  Toplot<-lapply(Data,function(toplot){return(toplot[,which(colnames(toplot)>=1989),])})
  Toplot<-lapply(Toplot,function(tr){
    ret<-lapply(1:dim(tr)[3],function(rep){return(apply(tr[,,rep],2,function(col){col<-col}))})#-tr[,1,rep]
    ret<-array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),length(ret)),
               dimnames=list(rownames(ret[[1]]),as.numeric(colnames(ret[[1]]))-1984,1:length(ret)))
    ret<-lapply(c(0.025,0.5,0.975),function(quant){return(apply(ret,c(1,2),function(x){return(quantile(x,probs=quant))}))})
    return(array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),3),
                 dimnames=list(rownames(ret[[1]]),colnames(ret[[1]]),c(0.025,0.5,0.975))))})
  
  plot(colnames(Toplot[[1]]),Toplot[[1]][1,,1],type="n",xaxt="n",
       xlab="",ylab="",ylim=c(min(unlist(Toplot),na.rm=T),max(unlist(Toplot),na.rm=T)))
  axis(1,at=as.character(seq(5,33,5)),labels=TRUE)  
  
  invisible(lapply(1:4,function(t){  
    toplot<-Toplot[[t]]
    
    invisible(lapply(1:nrow(toplot),function(i){
      absc<-colnames(Toplot[[t]])
      lines(absc,Toplot[[t]][i,,"0.5"], col = ColorsTr[[t]],lty = 1,lwd=2)
      polygon(c(absc,rev(absc)),c(Toplot[[t]][i,,"0.025"],rev(Toplot[[t]][i,,"0.975"])),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    }))
  }))
}

CWMdraw<-function(Cwm){
  par(mfrow=c(2,4),mar=c(2,2,3,1),oma=c(2,1,2,1),no.readonly = T)
invisible(lapply(colnames(Cwm[[1]]),function(trait){
  Toplot<-lapply(Cwm,function(pl){return(t(pl[,trait,]))})
  plot(colnames(Toplot[[1]]),Toplot[[1]]["0.5",], ylim=c(min(unlist(Toplot)),max(unlist(Toplot))),type="n",xlab="years",ylab="")
  #mtext(trait,3,cex=0.8,adj=0,line=0.5)
  
  invisible(lapply(1:4,function(tr){
    toplot<-Toplot[which(names(Toplot)%in%treatments[[tr]])]
    invisible(lapply(toplot,function(plo){
      lines(colnames(plo),plo["0.5",],col=ColorsTr[tr],lwd=1.5)
      polygon(c(colnames(plo),rev(colnames(plo))),c(plo["0.025",],rev(plo["0.975",])),
              col=rgb(0,0,0,alpha=0.05),border=NA)
      
    }))
  }))
}))
mtext("Community Weighted Means",line=0.5,adj=0,outer=TRUE,cex=1.1)
mtext("Years since disturbance",side=1,line=1.2,adj=0.65,cex=0.9,outer=TRUE)
}

legendCWM<-function(){
  mtext("Leaf thickness\n",at=0.13,line=-2.5,outer=TRUE,cex=0.9)
  mtext(expression(paste(mu, "m",sep = "")),at=0.08,line=-3,outer=TRUE,cex=0.9)
  mtext("Leaf cholophyll content\n",at=0.4,line=-2.5,outer=TRUE,cex=0.9)
  mtext(expression(paste("g.",mm^-2,sep = "")),at=0.34,line=-3,outer=TRUE,cex=0.9)
  mtext("Leaf toughness\n",at=0.64,line=-2.5,outer=TRUE,cex=0.9)
  mtext("N",at=0.56,line=-3,outer=TRUE,cex=0.9)
  mtext("SLA\n",at=0.88,line=-2.5,outer=TRUE,cex=0.9)
  mtext(expression(paste(mm^2,".",mg^-1,sep = "")),at=0.84,line=-3,outer=TRUE,cex=0.9)
  
  mtext("WD\n",at=0.13,line=-17.8,outer=TRUE,cex=0.9)
  mtext(expression(paste("g.",cm^-3,sep = "")),at=0.08,line=-18,outer=TRUE,cex=0.9)
  mtext("Bark thickness\n",at=0.4,line=-17.8,outer=TRUE,cex=0.9)
  mtext("mm",at=0.32,line=-18,outer=TRUE,cex=0.9)
  mtext("Hmax\n",at=0.64,line=-17.8,outer=TRUE,cex=0.9)
  mtext("m",at=0.56,line=-18,outer=TRUE,cex=0.9)
}

SeedMassProp<-function(SeedMass){
  par(mfrow=c(1,5),mar=c(1,1,2,1),oma=c(2,1,4,1),no.readonly = T)
  invisible(lapply(1:5,function(clas){
    Smass_class<-lapply(SeedMass,function(pl){return(t(pl[clas,,]))})
    plot(colnames(Smass_class[[1]]),Smass_class[[1]]["0.5",], ylim=c(min(unlist(Smass_class)),max(unlist(Smass_class))),
         type="n",xlab="years",ylab="")
    mtext(paste("class",clas),3,cex=1,line=1)
    invisible(lapply(1:4,function(tr){
      invisible(lapply(treatments[[tr]],function(plo){
        lines(colnames(Smass_class[[plo]]),Smass_class[[plo]]["0.5",],col=ColorsTr[tr],lwd=2)
        polygon(c(colnames(Smass_class[[plo]]),rev(colnames(Smass_class[[plo]]))),
                c(Smass_class[[plo]]["0.05",],rev(Smass_class[[plo]]["0.975",])),
                col=rgb(0,0,0,alpha=0.1),border=NA)}))
    }))
  }))
mtext("Proportion of seed mass classes",line=2,adj=0,outer=TRUE,cex=0.9)
mtext("Years since disturbance",side=1,line=1.2,adj=1,cex=0.9,outer=TRUE)
}

RedundancyPlot<-function(Red){
  Red<-lapply(Red,function(rep){
    ret<-smooth(rep,2)
   colnames(ret)<-colnames(rep)
   ret<-apply(ret[,which(as.numeric(colnames(ret))>="1989")],2,function(col){return(col)})#-ret[,"1989"]
   colnames(ret)<-as.numeric(colnames(ret))-1984
   return(ret)})
  
  Red<-array(unlist(Red),dim=c(12,ncol(Red[[1]]),length(Red)),
             dimnames=list(1:12,colnames(Red[[1]]),1:length(Red)))
  Red<-lapply(c(0.025,0.5,0.975),function(quant){return(apply(Red,c(1,2),function(rep){
    return(quantile(rep,probs=quant))}))})
  Red<-array(unlist(Red),dim=c(12,ncol(Red[[1]]),3),
             dimnames=list(1:12,colnames(Red[[1]]),c(0.025,0.5,0.975)))
  
  plot(colnames(Red),Red[1,,1],type='n',ylim=c(min(Red),max(Red)),xlab="",ylab="")
  
  invisible(lapply(1:4,function(tr){
    toplot<-Red[which(rownames(Red)%in%treatments[[tr]]),,]
    apply(toplot,1,function(li){
      lines(colnames(Red),li[,"0.5"],col=ColorsTr[tr],lwd=2)
      polygon(c(colnames(Red),rev(colnames(Red))),
              c(li[,"0.025"],rev(li[,"0.975"])),
              col=rgb(0,0,0,alpha=0.1),border=NA)
    })}))
}

plotPCA<-function(DataAcp){ 
  
  DataAcp_indiv<-DataAcp[["Indiv"]]
  DataAcp_traits<-DataAcp[["Traits"]]
  
  par(mfrow=c(2,3),mar=c(5,8,3,3))
  layout(mat=matrix(c(1,1,2,1,1,3),2,3,byrow=T))
  plot(DataAcp_indiv[,"Axis1"],DataAcp_indiv[,"Axis2"],type="n",xlab="",ylab="",frame.plot=F)
  points(DataAcp_indiv[,"Axis1"],DataAcp_indiv[,"Axis2"],
         col=terrain.colors(length(unique(DataAcp_indiv[,"name"])),alpha=0.4),pch=20)
  abline(h=0,v=0)
  mtext("Axis 1",1,at=0,line=2.4)
  mtext(paste(round(DataAcp[["Eigen"]][1]),"% of variance",sep=""),1,at=0,line=3.5,cex=0.8)
  mtext("Axis 2",2,at=0,las=1,line=2.5)
  mtext(paste(round(DataAcp[["Eigen"]][2]),"% of variance",sep=""),2,at=-1,las=1,line=0.5,cex=0.8)
  sp<-aggregate(DataAcp_indiv[,c("Axis1","Axis2")],list(DataAcp_indiv$name),median)
  
  outliers<-sp[which(sp[,"Axis1"]>=quantile(sp[,"Axis1"],0.995) |
                       sp[,"Axis1"]<=quantile(sp[,"Axis1"],0.005) |
                       sp[,"Axis2"]>=quantile(sp[,"Axis2"],0.995) | 
                       sp[,"Axis2"]<=quantile(sp[,"Axis2"],0.005)),]
  thigmophobe.labels(outliers[,"Axis1"],outliers[,"Axis2"],labels=outliers[,"Group.1"])
  title(main="(a) Individuals in the main PCA plan",adj=0)
  
  par(mar=c(2,4,2,2))
  barplot(DataAcp[["Eigen"]],col=c("black","black",rep("white",length(DataAcp[["Eigen"]])-2)))
  title(main="(b) Explained variance (%)",adj=0,cex.main=0.9)
  
  plot(DataAcp_traits[,"Comp1"],DataAcp_traits[,"Comp2"],type="n",xlab="",ylab="",frame.plot=F)
  text(DataAcp_traits[,"Comp1"],DataAcp_traits[,"Comp2"],labels=rownames(DataAcp_traits))
  abline(h=0,v=0)
  arrows(0,0,x1=DataAcp_traits[,"Comp1"], y1=DataAcp_traits[,"Comp2"], col="grey", length=0.1)
  title(main="(c) Traits in the main PCA plan",adj=0,cex.main=0.9)
}

### Previous graphs
TaxoTraj<-function(CompTaxo){
  par(mfrow=c(1,3),mar=c(5,2,4,2),oma=c(1,1.5,1,1),no.readonly=TRUE)
  for(q in 1:3){     
    Toplot<-lapply(CompTaxo,function(tr){return(tr[,,,q])})
    Toplot<-lapply(Toplot,function(toplot){return(toplot[,which(colnames(toplot)>=1989),])})
    Toplot<-lapply(Toplot,function(tr){
      ret<-lapply(1:dim(tr)[3],function(rep){return(apply(tr[,,rep],2,function(col){col<-col-tr[,1,rep]}))})#
      ret<-array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),length(ret)),
                 dimnames=list(rownames(ret[[1]]),as.numeric(colnames(ret[[1]]))-1984,1:length(ret)))
      ret<-lapply(c(0.025,0.5,0.975),function(quant){return(apply(ret,c(1,2),function(x){return(quantile(x,probs=quant))}))})
      return(array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),3),
                   dimnames=list(rownames(ret[[1]]),colnames(ret[[1]]),c(0.025,0.5,0.975))))})
    
    plot(colnames(Toplot[[1]]),Toplot[[1]][1,,1],type="n",xaxt="n",
         xlab="",ylab="",ylim=c(min(unlist(Toplot),na.rm=T),max(unlist(Toplot),na.rm=T)))
    axis(1,at=as.character(seq(5,33,5)),labels=TRUE)  
    mtext(paste("(",c("a","b","c")[q],") ",c("Richness","Shannon","Simpson")[q],sep=""),
          line=1,side=3,adj=0)
    
    invisible(lapply(1:4,function(t){  
      toplot<-Toplot[[t]]
      
      invisible(lapply(1:3,function(i){
        absc<-colnames(Toplot[[t]])
        lines(absc,Toplot[[t]][i,,"0.5"], col = ColorsTr[[t]],lty = 1,lwd=2)
        polygon(c(absc,rev(absc)),c(Toplot[[t]][i,,"0.025"],rev(Toplot[[t]][i,,"0.975"])),
                col=rgb(0,0,0,alpha=0.1),border=NA)
      }))
    }))
  }
mtext("Years since disturbance",side=1,adj=1,cex=0.8,line=-2,outer=TRUE)
mtext("Equivalent diversity",side=2,padj=1,cex=0.8,line=1.5,outer=TRUE)}

FunTraj<-function(CompFun,remove=TRUE){
  if(remove){CompFun[[2]]<-CompFun[[2]][which(rownames(CompFun[[2]])!=7),,]}
  
  CompFun<-lapply(CompFun,function(tr){
    ret<-tr[,which(colnames(tr)>=1989),]
    #ret<-lapply(1:dim(ret)[3],function(rep){return(apply(ret[,,rep],2,function(col){col<-col-ret[,1,rep]}))})
    #ret<-array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),length(ret)),dimnames=list(rownames(ret[[1]]),colnames(ret[[1]]),1:length(ret)))
    ret<-lapply(c(0.025,0.5,0.975),function(quant){return(apply(ret,c(1,2),function(x){return(quantile(x,probs=quant))}))})
    return(array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),3),
                 dimnames=list(rownames(ret[[1]]),as.numeric(colnames(ret[[1]]))-1984,c(0.025,0.5,0.975))))})
  
  plot(as.numeric(colnames(CompFun[[1]])),CompFun[[1]][1,,"0.5"],type="n",xaxt="n",
       xlab="",ylab="",ylim=c(min(unlist(CompFun),na.rm=T),max(unlist(CompFun),na.rm=T)),cex.axis=0.7)
  axis(1,at=as.character(seq(5,33,5)),labels=TRUE,cex.axis=0.7)  
  mtext("Rao diversity",3,adj=0,line=1,cex=1.5) 
  mtext("Years since disturbance",side=1,adj=1,line=2)
  mtext("Equivalent diversity",side=2,padj=1,line=3)
  
  invisible(lapply(1:length(CompFun),function(t){  
    toplot<-CompFun[[t]]
    toplot05<-smooth(toplot[,,"0.5"],1)
    toplot25<-smooth(toplot[,,"0.025"],1)
    toplot75<-smooth(toplot[,,"0.975"],1)
    invisible(lapply(1:nrow(toplot),function(i){
      lines(colnames(toplot),toplot05[i,], col = ColorsTr[[t]],lty = 1,lwd=2)
      polygon(c(colnames(toplot),rev(colnames(toplot))),c(toplot75[i,],rev(toplot25[i,])),
              col=rgb(0,0,0,alpha=0.05),border=NA)
    }))
  }))
}
