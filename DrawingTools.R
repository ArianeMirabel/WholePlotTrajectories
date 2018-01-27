### Manuscript drawing functions
Library(c("shape"))

treatments<-list(c(1,6,11),c(2,7,9),c(3,5,10),c(4,8,12))
names(treatments)<-c("Control","T1","T2","T3")
ColorsTr<-c("darkolivegreen2","deepskyblue2","darkorange1","red2")

TaxoCompo<-function(Data_TaxoComp){
par(mfrow=c(1,2))
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
mtext("Euclidean distance from 1989 inventory",side=2,padj=0,line=3,cex=0.9)
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

TaxoDist<-function(Data_TaxoComp){par(mfrow=c(1,2))
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
mtext("Euclidean distance from 1989 inventory",side=2,padj=0,line=2,cex=0.8)
mtext("(a) Taxonomic composition",side=3,adj=0,line=0.5)
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
  mtext("(b) Functional composition",side=3,adj=0,line=0.5)
  mtext("Years since disturbance",side=1,adj=1,line=2,cex=0.8)
}

TaxoTraj<-function(CompTaxo){par(mfrow=c(1,3),oma=c(1,1,1,1),no.readonly=TRUE)
for(q in 1:3){     
  Toplot<-lapply(CompTaxo,function(tr){return(tr[,,,q])})
  Toplot<-lapply(Toplot,function(toplot){return(toplot[,which(colnames(toplot)>=1989),])})
  Toplot<-lapply(Toplot,function(tr){
    ret<-lapply(1:dim(tr)[3],function(rep){return(apply(tr[,,rep],2,function(col){col<-col-tr[,1,rep]}))})#
    ret<-array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),dim(tr)[3]),
               dimnames=list(rownames(ret[[1]]),colnames(ret[[1]]),c(0.025,0.5,0.975)))
    return(ret)})
  Toplot<-lapply(Toplot,function(tr){return(tr[,which(!is.na(Toplot[[2]][1,,"0.5"])),])})
  
  plot(as.numeric(colnames(Toplot[[1]]))-1984,Toplot[[1]][1,,1],type="n",xaxt="n",
       xlab="",ylab="",ylim=c(min(unlist(Toplot),na.rm=T),max(unlist(Toplot),na.rm=T)))
  axis(1,at=as.character(seq(5,33,5)),labels=TRUE)  
  mtext(paste("(",c("a","b","c")[q],") ",c("Richness","Shannon","Simpson")[q],sep=""),
        line=1,side=3,adj=0)
  
  invisible(lapply(1:4,function(t){  
    toplot<-Toplot[[t]]
    
    invisible(lapply(1:3,function(i){
      absc<-as.numeric(colnames(Toplot[[t]]))-1984
      lines(absc,Toplot[[t]][i,,"0.5"], col = ColorsTr[[t]],lty = 1,lwd=2)
      polygon(c(absc,rev(absc)),c(Toplot[[t]][i,,"0.025"],rev(Toplot[[t]][i,,"0.975"])),
              col=rgb(0,0,0,alpha=0.1),border=NA)
    }))
  }))
}
mtext("Years since disturbance",side=1,adj=1,cex=0.8,line=-2,outer=TRUE)
mtext("Equivalent diversity",side=2,padj=1,cex=0.8,line=-1,outer=TRUE)}

FunTraj<-function(CompFun){plot(as.numeric(colnames(CompFun[[1]])),CompFun[[1]][1,,"0.5"],type="n",xaxt="n",
     xlab="",ylab="",ylim=c(min(unlist(CompFun),na.rm=T),max(unlist(CompFun),na.rm=T)),cex.axis=0.7)
axis(1,at=as.character(seq(5,33,5)),labels=TRUE,cex.axis=0.7)  
mtext("Rao diversity",3,adj=0,line=1) 
mtext("Years since disturbance",side=1,adj=1,cex=0.8,line=2)
mtext("Equivalent diversity",side=2,padj=1,cex=0.8,line=3)

invisible(lapply(1:length(CompFun),function(t){  
  
  toplot<-CompFun[[t]]
  
  invisible(lapply(1:nrow(toplot),function(i){
    lines(colnames(toplot),toplot[i,,"0.5"], col = ColorsTr[[t]],lty = 1,lwd=2)
    polygon(c(colnames(toplot),rev(colnames(toplot))),c(toplot[i,,"0.975"],rev(toplot[i,,"0.025"])),
            col=rgb(0,0,0,alpha=0.1),border=NA)
  }))
}))}