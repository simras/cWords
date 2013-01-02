# Rscript lib/plot12.R 2120211539BRIC-miR-21.rnk /home/teichmann/2120211539BRIC-miR-21.rnk_8cuFV/
#
# 
# Simon Horskjær Rasmussen
# The Bioinformatics Centre
# University of Copenhagen
################################################3

# Arguments
labtextS = 1
axistextS = 1.2
ccc = 1.4
textS = 2
args = tail(commandArgs(),3)
print(args)
datfile = args[1]
path = args[2]
num_overview =  as.integer(args[3])

  # Plot as SVG or jpeg
  #dev = jpg or svg or pdf
  dev="jpg"
  if(dev == "svg"){
    library(RSvgDevice)
  }
print(paste("Dev", dev))
# functions
running_mean = function(nums){
     return(as.double(sum(nums))/length(nums))
}

# Make filenames
datfileNeg = paste(path,datfile,".neg.plotfile.dat", sep="")
datfilePos = paste(path,datfile,".pos.plotfile.dat", sep="")

# See if the plot data has been generated
bneg = system(paste("ls",datfileNeg))
bpos = system(paste("ls",datfilePos))

# Get plot data for down-regulating motifs
if(bpos != 2){
  datPos = read.table(file=datfilePos ,header=TRUE,row.names=1)
  zscrsP = datPos[,1]
  npos = length(datPos[,1])
  # Number of Genes
  mpos = length(datPos[1,])
  datPos = datPos[1:npos,2:mpos]
  # arrange graphs
  posOrder = sort.list(zscrsP,decreasing=T)
}

# Get plot data for down-regulating motifs
if(bneg != 2){
  datNeg = read.table(file=datfileNeg ,header=TRUE,row.names=1)
  zscrsN = datNeg[,1]
  nneg = length(datNeg[,1])
  mneg = length(datNeg[1,])
  datNeg = rev(datNeg[1:nneg,2:mneg])
  negOrder = sort.list(zscrsN,decreasing=T)
}

# Is the data file empty?
NoDataPos = is.na(datPos[1,1])
NoDataNeg = is.na(datNeg[1,1])

# Colors of the graph
colspos = rainbow(n=npos)
colsneg = rainbow(n=nneg)
colsOver = rainbow(n=num_overview)

# Plot for each word
#win = 30

num_overviewneg = num_overview
num_overviewpos = num_overview

if(num_overviewpos == 0 | num_overviewpos > npos){
  num_overviewpos = npos
}

if(num_overviewneg == 0 | num_overviewneg > nneg){
  num_overviewneg = nneg
}

for (i in 1:npos){ 
  if(bpos != 2 && !NoDataPos){
    if(dev == "svg"){
      devSVG(file=paste(path,"pos.plot_",i,".svg",sep=""),width=10,height=9)
    }
    if(dev == "jpg"){
      jpeg(file=paste(path,"pos.plot_",i,".jpg",sep=""),width=1600,height=1300,pointsize=175,res = 10,quality=100)
    }
    if(dev == "pdf"){
      pdf(file=paste(path,"pos.plot_",i,".pdf",sep=""),width=10,height=9)
    }
#    if(i == 1){
#    	 plot(as.numeric(datPos[posOrder[i],]),col = colspos[1],type="l",main=rownames(datPos[posOrder[i],]),ylab="Z-Score",xlab="Gene index",lwd=13)
#         posMean = c()
#	 for(k in 1:mpos){
#	   if (j < 1){
#	     j = 1
#	   }	   
#	   posMean = c(posMean,mean(as.numeric(datPos[posOrder[i],j:k])))
#	 }
#	 lines(posMean,col = "blue",main=rownames(datPos[posOrder[i],]),ylab="Z-Score",xlab="Gene index",lwd=13)
#    }else{
    if(dev == "svg"){
      plot(as.numeric(datPos[posOrder[i],]),col = colspos[i],type="l",ylab="Z-Score",xlab="Gene index",cex=ccc,cex.axis=axistextS,cex.main=textS,cex.lab=labtextS,main=paste("Rank",i-1,rownames(datPos[posOrder[i],])))
    }
    if(dev == "jpg"){
      plot(as.numeric(datPos[posOrder[i],]),col = colspos[i],type="l",main=paste("Rank",i-1,rownames(datPos[posOrder[i],])),ylab="Z-Score",xlab="Gene index",lwd=13)
    }
    if(dev == "pdf"){
      plot(as.numeric(datPos[posOrder[i],]),col = colspos[i],type="l",main=paste("Rank",i-1,rownames(datPos[posOrder[i],])),ylab="Z-Score",xlab="Gene index",lwd=2)
    }
#    }
    dev.off()
  }
}
for (i in 1:nneg){ 
    
  if(bneg != 2  && !NoDataNeg){
    if(dev == "svg"){
      devSVG(file=paste(path,"neg.plot_",i,".svg",sep=""),width=10,height=9)
    }
    if(dev == "jpg"){
      jpeg(file=paste(path,"neg.plot_",i,".jpg",sep=""),height=1300,width=1600,pointsize=175,res = 10,quality=100)
    }
    if(dev == "pdf"){
      pdf(file=paste(path,"neg.plot_",i,".pdf",sep=""),width=10,height=9)
    }
                                        #    if(i == 1){
                                        #      plot(as.numeric(datNeg[i,]),col = colsneg[1],type="l",main=rownames(datNeg[i,]),ylab="Z-Score",xlab="Gene index",lwd=13)
                                        #         negMean = c()
                                        #	 for(k in 1:mneg){
                                        #	   j = k - win
                                        #	   if (j < 1){
                                        #	     j = 1
                                        #	   }	   
                                        #	   negMean = c(negMean,mean(as.numeric(datNeg[i,j:k])))
#	 }
                                        #	 lines(negMean,col = "blue",main=rownames(datNeg[i,]),ylab="Z-Score",xlab="Gene index",lwd=13)
                                        #    }else{
    if(dev == "svg"){
      plot(as.numeric(datNeg[negOrder[i],]),col = colsneg[i],type="l",cex.axis=axistextS,cex.lab=labtextS,cex.main=textS,ylab="Z-Score",xlab="Gene index",cex=ccc,main=rownames(datNeg[negOrder[i],]))
    }
    if(dev == "pdf"){
      plot(as.numeric(datNeg[negOrder[i],]),col = colsneg[i],type="l",main=rownames(datNeg[negOrder[i],]),ylab="Z-Score",xlab="Gene index",lwd=2)
    }
    if(dev == "jpg"){
      plot(as.numeric(datNeg[negOrder[i],]),col = colsneg[i],type="l",main=rownames(datNeg[negOrder[i],]),ylab="Z-Score",xlab="Gene index",lwd=13)
    }
    dev.off()
  }
}
  
for(j in 1:2){
  if (j == 1){
    dev = "pdf"
  }else{
    dev = "jpg"
  }
  
  # Plot Landscape summary for motifs in up-regulation
  if(bpos != 2 && !NoDataPos){
    ylimitpos = c(min(as.matrix(datPos[posOrder[1:npos],])),max(as.matrix(datPos[posOrder[1:npos],])))
    if(dev == "svg"){
      devSVG(file=paste(path,"pos.plot",".svg",sep=""),height=9,width=10)
    plot(as.numeric(datPos[posOrder[1],]),col = colsOver[1],ylim=ylimitpos ,type="l",main = paste("Top ",num_overviewpos," words involved in up regulation",sep=""),xlab="Gene Index",ylab="Z-Score",cex=ccc,cex.axis=axistextS,cex.main=textS,cex.lab=labtextS)
    legend(mpos - (mpos/14),y=ylimitpos[2] - 0.1,legend=rownames(datPos[posOrder[1:num_overviewpos],]),colpos=colsOver[1:num_overviewpos],lty=c("solid"))
  }
    if(dev == "jpg"){
      jpeg(file=paste(path,"pos.plot",".jpg",sep=""),height=1300,width=1600,pointsize=175,res = 10,quality=100)
      plot(as.numeric(datPos[posOrder[1],]),col = colsOver[1],ylim=ylimitpos ,type="l",main = paste("Top ",num_overviewpos," words involved in up regulation",sep=""),xlab="Gene Index",ylab="Z-Score",lwd=13)
      legend(mpos - (mpos/8),y=ylimitpos[2] - 0.1,legend=rownames(datPos[posOrder[1:num_overviewpos],]),col=colsOver[1:num_overviewpos],lty=c("solid"),lwd=13)
    }
    if(dev == "pdf"){
      pdf(file=paste(path,"pos.plot",".pdf",sep=""),width=10,height=9)
      plot(as.numeric(datPos[posOrder[1],]),col = colsOver[1],ylim=ylimitpos ,type="l",main = paste("Top ",num_overviewpos," words involved in up regulation",sep=""),xlab="Gene Index",ylab="Z-Score",lwd=2)
      legend(mpos - (mpos/8),y=ylimitpos[2] - 0.1,legend=rownames(datPos[posOrder[1:num_overviewpos],]),col=colsOver[1:num_overviewpos],lty=c("solid"),lwd=2)
    }
    for (i in 2:num_overviewpos){ 
      if(dev == "svg"){
        lines(as.numeric(datPos[posOrder[i],]),col = colsOver[i],cex=ccc,cex.axis=labtextS,cex.main=textS,cex.lab=labtextS)
      }
      if(dev == "jpg"){
      lines(as.numeric(datPos[posOrder[i],]),col = colsOver[i],lwd=13)
    }
      if(dev == "pdf"){
        lines(as.numeric(datPos[posOrder[i],]),col = colsOver[i],lwd=2)
      }
    }
    dev.off()
  }
  
  # Plot Landscape summary for motifs in down-regulation
  if(bneg != 2 && !NoDataNeg){
    ylimitneg = c(min(as.matrix(datNeg[1:nneg,])),max(as.matrix(datNeg[1:nneg,])))
    if(dev == "svg"){
      devSVG(file=paste(path,"neg.plot",".svg",sep=""),height=9,width=10)
      plot(as.numeric(datNeg[negOrder[1],]),col = colsOver[1],ylim=ylimitneg ,type="l",main = paste("Top ",num_overviewneg," words involved in down regulation",sep=""),xlab="Gene Index",ylab="Z-Score",cex=ccc,cex.axis=axistextS,cex.main=textS,cex.lab=labtextS)
      legend(mneg-(mneg/14),y=ylimitneg[2] - 0.1,legend=rownames(datNeg[negOrder[1:num_overviewneg],]),col=colsOver[1:num_overviewneg],lty=c("solid"))
    }
    
    if(dev == "jpg"){
      jpeg(file=paste(path,"neg.plot",".jpg",sep=""),height=1300,width=1600,pointsize=175,res = 10,quality=100)
      plot(as.numeric(datNeg[negOrder[1],]),col = colsOver[1],ylim=ylimitneg ,type="l",main = paste("Top ",num_overviewneg," words involved in down regulation",sep=""),xlab="Gene Index",ylab="Z-Score",lwd=13)
      legend(mneg-(mneg/8),y=ylimitneg[2] - 0.1,legend=rownames(datNeg[negOrder[1:num_overviewneg],]),col=colsOver[1:num_overviewneg],lty=c("solid"),lwd=13)
    }
    if(dev == "pdf"){
      pdf(file=paste(path,"neg.plot",".pdf",sep=""),width=10,height=9)
      plot(as.numeric(datNeg[negOrder[1],]),col = colsOver[1],ylim=ylimitneg ,type="l",main = paste("Top ",num_overviewneg," words involved in down regulation",sep=""),xlab="Gene Index",ylab="Z-Score",lwd=2)
      legend(mneg-(mneg/8),y=ylimitneg[2] - 0.1,legend=rownames(datNeg[negOrder[1:num_overviewneg],]),col=colsOver[1:num_overviewneg],lty=c("solid"),lwd=2)
    }
  
    for (i in 2:num_overviewneg){
      if(dev == "svg"){    
        lines(as.numeric(datNeg[negOrder[i],]),col = colsOver[i],cex=ccc,cex.axis=labtextS,cex.main=textS,cex.lab=labtextS)
      }
      if(dev == "jpg"){
        lines(as.numeric(datNeg[negOrder[i],]),col = colsOver[i],lwd=13)
      }
      if(dev == "pdf"){
        lines(as.numeric(datNeg[negOrder[i],]),col = colsOver[i],lwd=2)
      }
    }
    dev.off()
  }
}
