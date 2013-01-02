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
  dev="pdf"
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
  norm = datPos[1,2]
  datPos = 2 * norm - rev(datPos[1:npos,2:mpos])

  # arrange graphs
  order = sort.list(zscrsP,decreasing=T)
  
}

# Get plot data for down-regulating motifs
if(bneg != 2){
  datNeg = read.table(file=datfileNeg ,header=TRUE,row.names=1)
  zscrsN = datNeg[,1]
  nneg = length(datNeg[,1])
  mneg = length(datNeg[1,])
  datNeg = rev(datNeg[1:nneg,2:mneg])
#  negOrder = sort.list(zscrsN,decreasing=T)
}
  zscrs = c(zscrsN,zscrsP)
print(datNeg)
print(datPos)
  dat = rbind(datNeg,datPos)
  order = sort.list(zscrs,decreasing=T)

# Is the data file empty?
NoData = is.na(dat[1,1])

# Colors of the graph
n=npos+nneg
cols = rainbow(n=(n))

# Plot for each word
win = 30

if(num_overview == 0 | num_overview > n){
  num_overview = n
}

# Plot Landscape summary for motifs in down-regulation
if((bpos != 2 || bneg != 2) && !NoData){
  ylimit = c(min(as.matrix(dat[1:n,])),max(as.matrix(dat[1:n,])))
  if(dev == "svg"){
    devSVG(file=paste(path,datfile,"_one.plot",".svg",sep=""),height=9,width=10)
    plot(as.numeric(dat[order[1],]),col = cols[1],ylim=ylimit ,type="l",main = paste("Top ",num_overview,", Word correlations in genes ordered down to up",sep=""),xlab="Gene Index",ylab="Z-Score",cex=ccc,cex.axis=axistextS,cex.main=textS,cex.lab=labtextS)
    legend(mneg-(mneg/14),y=ylimit[2] - 0.1,legend=rownames(dat[order[1:num_overview],]),col=cols[1:num_overview],lty=c("solid"))
  }
  if(dev == "jpg"){
    jpeg(file=paste(path,datfile,"_one.plot",".jpg",sep=""),height=1300,width=1600,pointsize=175,res = 10,quality=100)
    plot(as.numeric(dat[order[1],]),col = cols[1],ylim=ylimit ,type="l",main = paste("Top ",num_overview,", Word correlations in genes ordered down to up",sep=""),xlab="Gene Index",ylab="Z-Score",lwd=13)
    legend(mneg-(mneg/8),y=ylimit[2] - 0.1,legend=rownames(dat[order[1:num_overview],]),col=cols[1:num_overview],lty=c("solid"),lwd=13)
  }
  if(dev == "pdf"){
    pdf(file=paste(path,datfile,"_one.plot",".pdf",sep=""),width=10,height=9)
    plot(as.numeric(dat[order[1],]),col = cols[1],ylim=ylimit ,type="l",main = paste("Top ",num_overview,", Word correlations in genes ordered down to up",sep=""),xlab="Gene Index",ylab="Z-Score",lwd=2)
    legend(mneg-(mneg/8),y=ylimit[2] - 0.1,legend=rownames(dat[order[1:num_overview],]),col=cols[1:num_overview],lty=c("solid"),lwd=2)
  }

  for (i in 2:num_overview){
    if(dev == "svg"){    
      lines(as.numeric(dat[order[i],]),col = cols[i],cex=ccc,cex.axis=labtextS,cex.main=textS,cex.lab=labtextS)
    }
    if(dev == "jpg"){
      lines(as.numeric(dat[order[i],]),col = cols[i],lwd=13)
    }
    if(dev == "pdf"){
      lines(as.numeric(dat[order[i],]),col = cols[i],lwd=2)
    }

  }
  dev.off()
}
