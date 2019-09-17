library("shape")
library("plotrix")

treatments<-list(c(1,6,11),c(2,7,9),c(3,5,10),c(4,8,12))
names(treatments)<-c("Control","Low","Intermediate","High")
ColorsTr<-c("darkolivegreen2","gold","orangered","darkred")
colyear<-c("deepskyblue","cornflowerblue","darkslateblue")
time<-c("1995","2005","2015")

smooth<-function(mat,larg){return(do.call(cbind,lapply(1:ncol(mat),function(step){
  range<-max(1,step-larg):min(ncol(mat),step+larg)
  rowSums(mat[,range])/length(range)})))}


load("DB/RichnessIDH");load("DB/SimpsonIDH");load("DB/RaoIDH");load("DB/LostAGB");load("DB/FunRichIDH")
  
windows()
  
Data<-Rich
AgbLoss<-AGBloss

Ylim<-c(min(unlist(lapply(Data, function(tr){return(tr[,,"0.5"])}))),
          max(unlist(lapply(Data, function(tr){return(tr[,,"0.5"])}))))

plot(AgbLoss[,"AGB"],AgbLoss[,"AGB"],type="n",xlab="",ylab="",
       ylim=Ylim,xaxt='n',cex.axis=1.5)
Axis(side=1, labels=FALSE)
  
Ti<-2

toplot<-unlist(lapply(Data, function(tr){return(tr[,time[Ti],"0.5"])}))

abs<-AgbLoss[which(AgbLoss[,"plot"]%in%names(toplot)),"AGB"]
points(abs,toplot,pch=20)
Lm<-lm(toplot~abs)
Lm2<-lm(toplot~abs+I(abs^2))

abs_pred<-seq(min(abs),max(abs),length.out=100)
lines(sort(abs_pred),predict(Lm2,newdata=data.frame(abs=abs_pred)),lwd=5)





windows()

Data<-FunRich
AgbLoss<-AGBloss

Ylim<-c(min(unlist(lapply(Data, function(tr){return(tr[,,"0.5"])}))),
        max(unlist(lapply(Data, function(tr){return(tr[,,"0.5"])}))))

plot(AgbLoss[,"AGB"],AgbLoss[,"AGB"],type="n",xlab="",ylab="",
     ylim=Ylim,xaxt='n',cex.axis=1.5)
Axis(side=1, labels=FALSE)

Ti<-2

toplot<-unlist(lapply(Data, function(tr){return(tr[,time[Ti],"0.5"])}))

abs<-AgbLoss[which(AgbLoss[,"plot"]%in%names(toplot)),"AGB"]
points(abs,toplot,pch=20)
Lm<-lm(toplot~abs)
Lm2<-lm(toplot~abs+I(abs^2))

abs_pred<-seq(min(abs),max(abs),length.out=100)
lines(sort(abs_pred),predict(Lm2,newdata=data.frame(abs=abs_pred)),lwd=5)
