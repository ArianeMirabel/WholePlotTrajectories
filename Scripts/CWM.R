library(vegan)
library(reshape2)
library(ade4)
source("Scripts/Vernacular_handle.R")
load("DB/Alpha_Plots")
load("DB/Paracou_R_Subdivided_ok")

Traits1<-read.csv("DB/BridgeOK.csv",sep=";",na.strings="")
Traits2<-read.csv("DB/DataLifeTraits.csv",sep=";",na.strings="")
source("Scripts/TraitsMiceFilling.R")

T0<-c(1,6,11)
T1<-c(2,7,9)
T2<-c(3,5,10)
T3<-c(4,8,12)

treatments<-list(T0,T1,T2,T3)
names(treatments)<-c("Control","T1","T2","T3")

TraitsName<-c("L_thickness","L_chloro","L_toughness","SLA","WD","Bark_thick","Hmax")

ColorsTr<-c("darkolivegreen2","deepskyblue2","darkorange1","red2")#"chartreuse3",

InventorySp<-do.call(rbind,lapply(LivingStand_all,function(yr){
  ret<-do.call(rbind,yr)[,c("Famille","Genre","name")]
  ret<-ret[which(!grepl("Indet.",ret[,"name"])),]
  return(ret[which(!duplicated(ret)),])}))
InventorySp<-InventorySp[which(!duplicated(InventorySp)),]

dates<-sort(names(LivingStand_all))[-1]

Nrep<-2
CWM<-lapply(1:12,function(p){
  matplot<-lapply(1:Nrep,function(rep){
    Trajnmds<-lapply(dates,function(y){
      temp<-LivingStand_all[[which(names(LivingStand_all)==y)]]
      if(!any(names(temp)==p)){return(NA)}
      if(any(names(temp)==p)){
        temp<-temp[[which(names(temp)==p)]]
        temp<-temp[!duplicated(temp),]
        #temprep<-temp[which(temp[,"name"]!="Indet."),"name"]
        temprep<-Replacement(temp,Alpha=alpha_construct(temp))
        #temprep<-as.data.frame(Replacement(temp,Alpha=alpha_construct(temp)));colnames(temprep)<-"name"
        #temp<-temp[,c("Genre","name")];temp<-as.data.frame(temp[which(!duplicated(temp)),],col.names="name")
        #temprep<-merge(temprep,temp,by="name",all.x=T);temprep<-as.character(temprep[,"Genre"])
        return(as.ProbaVector(tapply(temprep,temprep,length)))}})
    
    names(Trajnmds)<-dates
    Ok<-names(Trajnmds)[which(unlist(lapply(Trajnmds,function(yr){return(any(!is.na(yr)))})))]
    Trajnmds<-Trajnmds[which(names(Trajnmds)%in%Ok)]
    Ntrees<-unlist(lapply(Trajnmds,length))
    Ok<-names(which(unlist(lapply(2:length(Trajnmds),function(step){return(Ntrees[step]-Ntrees[step-1])}))<=-90))
    Trajnmds<-Trajnmds[which(!names(Trajnmds)%in%Ok)]
    
    init<-Reduce(c,lapply(Trajnmds,function(yr){return(names(yr))}))
    init<-as.data.frame(unique(init),row.names=unique(init))
    
    namesOk<-names(Trajnmds)
    
    Trajnmds<-do.call(cbind,lapply(Trajnmds,function(yr){
      ret<-merge(init,yr,by="row.names",all=TRUE)
      ret[is.na(ret)]<-0;rownames(ret)<-ret[,"Row.names"]
      return(ret[3])}))
    colnames(Trajnmds)<-namesOk
    
    #traits_filled<-Traits_filling(Traits1,Traits2,InventorySp)
    traits_filled<-aggregate(Traits_filled[,TraitsName],list(traits_filled$name),median)
    rownames(traits_filled)<-traits_filled[,1];traits_filled<-traits_filled[,TraitsName]
    
    tra<-traits_filled[which(rownames(traits_filled)%in%rownames(Trajnmds)),];tra<-tra[order(rownames(tra)),]
    Trajnmds<-Trajnmds[which(rownames(Trajnmds)%in%rownames(tra)),];Trajnmds<-Trajnmds[order(rownames(Trajnmds)),]
    return(t(apply(Trajnmds,2,function(Inv){Inv<-colSums(apply(tra,2,function(col){col<-col*Inv}))})))
    })
  if(length(unique(unlist(lapply(matplot,function(rep){return(nrow(rep))}))))!=1){
    print(paste(p,y))
    }
  return(array(unlist(matplot),dim=c(nrow(matplot[[1]]),ncol(matplot[[1]]),Nrep),
      dimnames=list(rownames(matplot[[1]]),colnames(matplot[[1]]),1:Nrep)))
})
names(CWM)<-1:12

CWM<-lapply(CWM,function(pl){
  ret<-lapply(c(0.025,0.5,0.975),function(quant){apply(pl,c(1,2),function(x){return(quantile(x,probs=quant))})})
  return(array(unlist(ret),dim=c(nrow(pl),ncol(pl),3),dimnames=list(rownames(pl),colnames(pl),c(0.025,0.5,0.975))))
})

windows()
CWMdraw<-function(Cwm){
  par(mfrow=c(2,4),mar=c(2,2,3,1),oma=c(2,1,2,1),no.readonly = T)
  invisible(lapply(colnames(Cwm[[1]]),function(trait){
    Toplot<-lapply(Cwm,function(pl){return(t(pl[,trait,]))})
    plot(colnames(Toplot[[1]]),Toplot[[1]]["0.5",], ylim=c(min(unlist(Toplot)),max(unlist(Toplot))),type="n",xlab="years",ylab="")
    mtext(trait,3,cex=0.8,adj=0,line=0.5)
    
    invisible(lapply(1:4,function(tr){
      toplot<-Toplot[which(names(Toplot)%in%treatments[[tr]])]
      invisible(lapply(toplot,function(plo){
        lines(colnames(plo),plo["0.5",],col=ColorsTr[tr],lwd=2)
        polygon(c(colnames(plo),rev(colnames(plo))),c(plo["0.025",],rev(plo["0.975",])),
                col=rgb(0,0,0,alpha=0.1),border=NA)
        
      }))
    }))
  }))
  mtext("Community Weighted Means",line=0.5,adj=0,outer=TRUE,cex=1.1)
  mtext("Years since disturbance",side=1,line=1.2,adj=0.65,cex=0.9,outer=TRUE)
}


save(CWM_genus,file="DB/CWM_genus")

Naeffectives<-unlist(lapply(TraitsName,function(trait){return(length(which(is.na(traits[,trait])))/nrow(traits))}))*100
names(Naeffectives)<-TraitsName

Nrep<-10
Traits2["name"]<-sub(" ","_",Traits2[,"Name"])
seedM<-Traits2[which(!is.na(Traits2[,"Masse"])),c("name","Masse")]
rownames(seedM)<-seedM[,"name"];seedM<-seedM["Masse"]

#Do with mass class percentage
Smass<-lapply(1:12,function(p){
  matplot<-lapply(1:Nrep,function(rep){
    Trajnmds<-lapply(dates,function(y){
      temp<-LivingStand_all[[which(names(LivingStand_all)==y)]]
      if(!any(names(temp)==p)){return(NA)}
      if(any(names(temp)==p)){
        temp<-temp[[which(names(temp)==p)]]
        temp<-temp[!duplicated(temp),]
        temp<-Replacement(temp,Alpha=alpha_construct(temp))
        return(as.ProbaVector(tapply(temp,temp,length)))}})
    
    names(Trajnmds)<-dates
    Ok<-names(Trajnmds)[which(unlist(lapply(Trajnmds,function(yr){return(any(!is.na(yr)))})))]
    Trajnmds<-Trajnmds[which(names(Trajnmds)%in%Ok)]
    Ntrees<-unlist(lapply(Trajnmds,length))
    Ok<-names(which(unlist(lapply(2:length(Trajnmds),function(step){return(Ntrees[step]-Ntrees[step-1])}))<=-90))
    Trajnmds<-Trajnmds[which(!names(Trajnmds)%in%Ok)]
    
    init<-Reduce(c,lapply(Trajnmds,function(yr){return(names(yr))}))
    init<-as.data.frame(unique(init),row.names=unique(init))
    
    namesOk<-names(Trajnmds)
    
    Trajnmds<-do.call(cbind,lapply(Trajnmds,function(yr){
      ret<-merge(init,yr,by="row.names",all=TRUE)
      ret[is.na(ret)]<-0;rownames(ret)<-ret[,"Row.names"]
      return(ret[3])}))
    colnames(Trajnmds)<-namesOk
    
    ret<-merge(Trajnmds,seedM,by="row.names",all=F)
    rownames(ret)<-ret[,"Row.names"];ret<-ret[,-1]
  
    ret<-do.call(rbind,lapply(sort(unique(ret[,"Masse"])),function(clas){
      apply(ret[which(ret[,"Masse"]==clas),-ncol(ret)],2,sum)}))
    rownames(ret)<-1:5
    return(ret)
    })
  matplot<-array(unlist(matplot),dim=c(5,ncol(matplot[[1]]),Nrep),
               dimnames=list(1:5,as.numeric(colnames(matplot[[1]]))-1984,1:Nrep))
  matplot<-lapply(c(0.025,0.5,0.975),function(quant){
    return(apply(matplot,c(1,2),function(rep){return(quantile(rep,probs=quant))}))})
  matplot<-array(unlist(matplot),dim=c(5,ncol(matplot[[1]]),3),
                 dimnames=list(1:5,colnames(matplot[[1]]),c(0.05,0.5,0.975)))
})

Smass<-lapply(Smass,function(pl){return(pl[,which(colnames(pl)%in%colnames(Smass[[2]])),])})

save(Smass,file="DB/SeedMassClasses_proportion")
windows()
par(mfrow=c(2,3))
lapply(1:5,function(clas){
  Smass_class<-lapply(Smass,function(pl){return(t(pl[clas,,]))})
  plot(colnames(Smass_class[[1]]),Smass_class[[1]]["0.5",], ylim=c(min(unlist(Smass_class)),max(unlist(Smass_class))),
     type="n",xlab="years",ylab="")
mtext(paste("Seed Mass, Proportion class",clas),3,cex=1.3,line=1)
invisible(lapply(1:4,function(tr){
  invisible(lapply(treatments[[tr]],function(plo){
    lines(colnames(Smass_class[[plo]]),Smass_class[[plo]]["0.5",],col=ColorsTr[tr],lwd=3)
  polygon(c(colnames(Smass_class[[plo]]),rev(colnames(Smass_class[[plo]]))),
          c(Smass_class[[plo]]["0.05",],rev(Smass_class[[plo]]["0.975",])),
          col=rgb(0,0,0,alpha=0.1),border=NA)}))
}))
})



load("DB/CWM_genus")

CWM<-CWM_genus

lapply(CWM,function(pl){return(list(unique(unlist(lapply(pl,nrow))),unique(unlist(lapply(pl,ncol)))))})
pl<-CWM[[5]]
dates<-rownames(pl[[5]])
unlist(lapply(pl,function(rep){setdiff(dates,rownames(rep))}))
CWM[[5]]<-lapply(CWM[[5]],function(rep){return(rep[which(!rownames(rep)%in%c("2000","2002")),])})
CWM[[7]]<-lapply(CWM[[7]],function(rep){return(rep[which(!rownames(rep)%in%c("2000","2002")),])})


CWM<-lapply(CWM,function(pl){
  pl<-lapply(pl,function(rep){return(rep[which(!rownames(rep)%in%c("1998","2000","2002")),])})
  #ret<-lapply(pl,function(rep){return(t(apply(rep,1,function(li){return(li-rep[6,])})))})
  #ret<-array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),length(ret)))
  ret<-array(unlist(pl),dim=c(nrow(pl[[1]]),ncol(pl[[1]]),length(pl)))
  ret<-lapply(c(0.025,0.5,0.975),function(quant){apply(ret,c(1,2),function(x){return(quantile(x,probs=quant))})})
  return(array(unlist(ret),dim=c(nrow(ret[[1]]),ncol(ret[[1]]),3),
          dimnames=list(as.numeric(rownames(pl[[1]]))-1984,colnames(pl[[1]]),c(0.025,0.5,0.975))))
  })
CWM<-lapply(CWM,function(pl){return(pl[6:nrow(pl),,])})

windows()
par(mfrow=c(2,4))
invisible(lapply(colnames(CWM[[1]]),function(trait){
Toplot<-lapply(CWM,function(pl){return(t(pl[,trait,]))})
plot(colnames(Toplot[[1]]),Toplot[[1]]["0.5",], ylim=c(min(unlist(Toplot)),max(unlist(Toplot))),type="n",xlab="years",ylab="")
legend("bottomright",c("Control","T1","T2","T3"),lty=1,lwd=3,col=ColorsTr,bty="n")
mtext(trait,3,cex=1.3,line=2)

invisible(lapply(1:4,function(tr){
    toplot<-Toplot[which(names(Toplot)%in%treatments[[tr]])]
  invisible(lapply(toplot,function(plo){
    lines(colnames(plo),plo["0.5",],col=ColorsTr[tr],lwd=3)
    polygon(c(colnames(plo),rev(colnames(plo))),c(plo["0.025",],rev(plo["0.975",])),
            col=rgb(0,0,0,alpha=0.1),border=NA)
    
  }))
}))
}))


save(CWM,file="P:/Private/Taff/These/Redaction/2_WholePlotTrajectories/DB/CWM")

windows(width=40,height=20)
par(mfrow=c(2,4))
invisible(lapply(colnames(CWM[[1]]),function(trait){
  
  cwm<-lapply(CWM[1:12],function(pl){
    rownames(pl)<-as.numeric(rownames(pl))-1984
    pl<-pl[which(as.numeric(rownames(pl))>=10),,]
    return(pl)
  })
  absc<-as.numeric(rownames(cwm[[12]]))
  cwm<-lapply(cwm,function(pl){
    return(pl[which(rownames(pl)%in%absc),,])
  })
  
  Ylim<-c(min(unlist(lapply(cwm,function(pl){return(pl[,trait,])}))),
          max(unlist(lapply(cwm,function(pl){return(pl[,trait,])}))))
  cwm<-lapply(cwm,function(pl){
    ret<-lapply(c(0.025,0.5,0.975),function(quant){
      return(apply(pl,c(1,2),function(x){return(quantile(x,probs=quant))}))})
    names(ret)<-c(0.025,0.5,0.975)
    return(ret)})
  
  plot(absc,cwm[[8]][[0.5]][,trait], ylim=Ylim,type="n",xlab="years",ylab="")
  #legend("bottomright",c("Control","T1","T2","T3"),lty=1,lwd=3,col=ColorsTr[c(1,3:5)],bty="n")
  mtext(trait,3,cex=1.3,line=1)
  
  invisible(lapply(1:4,function(tr){
    invisible(lapply(cwm[treatments[[tr]]],function(plo){
      lines(absc,plo[["0.5"]][,trait],col=ColorsTr[tr],lwd=3)
      polygon(c(absc,rev(absc)),c(plo[["0.975"]][,trait],rev(plo[["0.025"]][,trait])),
              col=rgb(0,0,0,alpha=0.1),border=NA)
    }))
  }))
  
}))

#By treatment
windows(width=40,height=20)
par(mfrow=c(2,4))
invisible(lapply(colnames(CWM[[1]]),function(trait){
  
  cwm<-lapply(treatments,function(tr){
    pl<-do.call(cbind,lapply(CWM[tr],function(plo){return(plo[,trait,])}))
    rownames(pl)<-as.numeric(rownames(pl))-1984
    pl<-pl[which(as.numeric(rownames(pl))>=10),]
    return(pl)
  })
  
  cwm<-lapply(cwm,function(tr){return(t(apply(tr,1,function(li){return(li-tr[1,])})))})
  
  Ylim<-c(min(unlist(cwm)),max(unlist(cwm)))
  cwm<-lapply(cwm,function(pl){
    ret<-do.call(rbind,lapply(c(0.025,0.5,0.975),function(quant){
      return(apply(pl,1,function(x){return(quantile(x,probs=quant))}))}))
    rownames(ret)<-c(0.025,0.5,0.975)
    return(ret)})
  
  plot(colnames(cwm[[1]]),cwm[[1]]["0.5",], ylim=Ylim,type="n",xlab="years",ylab="")
  #legend("bottomright",c("Control","T1","T2","T3"),lty=1,lwd=3,col=ColorsTr[c(1,3:5)],bty="n")
  mtext(trait,3,cex=1.3,line=1)
  
  invisible(lapply(1:4,function(tr){
      lines(colnames(cwm[[tr]]),cwm[[tr]]["0.5",],col=ColorsTr[tr],lwd=3)
      polygon(c(colnames(cwm[[tr]]),rev(colnames(cwm[[tr]]))),c(cwm[[tr]]["0.975",],rev(cwm[[tr]]["0.025",])),
              col=rgb(0,0,0,alpha=0.1),border=NA)
  }))
  
}))

