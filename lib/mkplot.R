# ################################################################
# mkplot.R - Generates Words Cluster Plot
# Produces the Words Cluster Plot in cWords
# 2011 - Simon Horskjaer Rasmussen
# Center For Bioinformatics - University of Copenhagen
##################################################################

sites = function(str,lt) {
# Extracts seed sites of length 6, 7 and 8 of an 8mer seed site
# str is the 8mer seed site and lt commaseparated list of lengths
len = as.numeric(strsplit(lt,",")[[1]])
   out = c()
  if(nchar(str) < 8){
    return(c())
  }
  strlist = strsplit(str,"")[[1]]
  l = length(strlist)
  if(length(which(len==6)) > 0){
    out = c(out,paste(strlist[1:(l-2)],collapse=""),paste(strlist[2:(l-1)],collapse=""))
  }
  if(length(which(len==7)) > 0){
    out = c(out,paste(strlist[1:(l-1)],collapse=""),paste(strlist[2:l],collapse=""))
  }
  if(length(which(len==8)) > 0){
    out = c(out,str)
  }
return(out)
}

ledge_plot = function(fileName,set,seed,len,ngenes,path,base_path,species,annotations,dev){
  # Function that generates the plot 
  # Arguments
  labtextS = 1
  axistextS = 1.2
  ccc = 1.4

  # Plot as SVG or jpeg
  #dev = jpg or svg or pdf or both (jpg and pdf)
  #devv = "both"
  #dev1 = "pdf"
  #dev2 = "jpg"
  if(dev == "svg"){
    library(RSvgDevice)
  }
  # print(annotations)
  # Read annotation file put it into a hash
  if(annotations != ""){
    #annotFName = paste(base_path,"resources/annotations_",species,".tsv",sep="")
    annotFName = annotations
#    print(annotFName)
    # different annotation files
    # annotFName = "resources/annot_cons.tsv"
    # annotFName = paste(base_path,"resources/annot_gap2.tsv",sep="")
    annotTable = read.table(annotFName)
    aHash = new.env(hash=TRUE, parent=emptyenv(), size=100L)
    for(i in 1:length(annotTable[,1])){
      assign(as.character(annotTable[i,1]),as.character(annotTable[i,2]),aHash)
    }
  }else{
    aHash = new.env(hash=TRUE, parent=emptyenv(), size=100L)
  }
  # Layout variable
  Xdiv = 27
  # number of word to write in text in plot
  textWords = 100
  # UPGMA threshold
  tt = 0.36
  threshold = tt
  # Number of words to be clustered
  mtfs = 300
  # Make histogram?
  his = FALSE
  
  for(m in 1:length(len)){
    llen = gsub(pattern=",", replacement="",x=len[m])
    tab = try(read.table(paste(path,fileName,llen,".",set,".0",sep=""),fill=TRUE,sep=" "))
    if(class(tab)=="try-error"){
      print("Could not read the rank file in mkplot.R")
      return(0)
    }
    dat = as.matrix(tab[2:(length(tab[,1])-1),c(2,3,6)])
    
    if(length(dat[,1]) < 10){
      print("Too little data to plot")
      return(0)
    }

    if(dev == "svg"){
      devSVG(file=paste(path,fileName,".",llen,set,".svg",sep=""),width=10,height=9)
      plot(dat[,c(3,2)],xlim=c(0,ngenes),col=rgb(0,100,130,50,maxColorValue=255), cex = ccc,cex.lab = labtextS, cex.axis = axistextS, pch=16, xlab="Gene index", ylab="Z-score")
    }
    if(dev == "jpg"){
      jpeg(filename=paste(path,fileName,".",llen,set,".jpg",sep=""),height=1300,width=1600,pointsize=50,res = 30,quality=100)
      plot(dat[,c(3,2)],xlim=c(0,ngenes),col=rgb(0,100,130,100,maxColorValue=255), pch=16,xlab="Gene index",ylab="Z-score")
    }
                                        #print(paste("here......",dat))
    if(dev == "pdf"){
      pdf(file=paste(path,fileName,".",llen,set,".pdf",sep=""),height=9,width=10)
      plot(dat[,c(3,2)],xlim=c(0,ngenes),col=rgb(0,100,130,100,maxColorValue=255), pch=16,xlab="Gene index",ylab="Z-score")
    }
    
                                        # Make title
    if(set == "negative"){
      sset = "Down regulation"
    }else{
      sset = "Up regulation"
    }
    
    title = gsub("^[0-9]+","",fileName,perl=T)
    
    title(main=paste(sset," in ", title," for word length(s) ",len[m],sep=""))
    
    Ymax = as.numeric(dat[1,2])
    Ymin = as.numeric(tail(dat[,2],1))
      
                                        # Calculate a one unit distance to be able to adjust plot according to different word sizes
    unitDist = (Ymax - Ymin)/24
                                        #print(base_path)
                                        # Do clustering and motif alignment
    system(paste("java -cp ",base_path,"lib Upgma ",path," ",mtfs," " ,0.40," ",set," ",fileName,llen,".",set,".0", " 1> /dev/null 2> /dev/null",sep=""))
    flist = list.files(path=path,patt=c(paste("motifs",set,".dat",sep="")))
                                        #print(paste("flist ",flist))
    if (length(flist) > 0){
      motifs = read.table(file = paste(path,"motifs",set,".dat",sep=""),fill = T,stringsAsFactors = F,quote="\"",sep=" ")
      system(paste("mv ",path,"motifs",set,".dat"," ",path,"motifs",set,".txt",sep=""))
                                        # print(motifs)
      numberOfMotifs = length(motifs[,1])
                                        # Write motif legend		    
      text(x=.8*ngenes,y=Ymax,pos=4,labels=paste("Top 10 Motifs of first",mtfs,"words",sep=" "))
      
      if(numberOfMotifs > 10){
        n = 10
      }else{
        n = numberOfMotifs
      }
      cols = rainbow(n=n)
      for(i in 1:n){
        les = c()
        zs = c()
        annotated = c()
        annotatedZ = c()
        mfs = as.character(motifs[i,!is.na(motifs[i,])])
        mfs = mfs[nchar(mfs) > 0 & !is.na(mfs)]
                                        #print(mfs)
        mm = length(mfs)
        if(mm == 0){
          break
        }
        for(j in 1:mm){
          lls = which(as.character(dat[,1]) == mfs[j])
          testmotif = as.character(mfs[j])
          if(exists(x=testmotif,env=aHash)){
                                        # marks annotated words
            llll = which(as.character(dat[,1]) == mfs[j])
            annotated = c(annotated,dat[llll,3])
            annotatedZ = c(annotatedZ,dat[llll,2])
          }
          if(length(lls) > 0){
            les = c(les,dat[lls,3])
            zs = c(zs,dat[lls,2])
          }
        }
        
        rs = sample(-100:100,3)
        
        if(dev == "svg"){
                                        # plot annotation and words as text
          points(x=as.numeric(les),y=as.numeric(zs),xlim=c(0,ngenes),col = cols[i],cex=1, pch=16)
                                        # Legend points
          points(x=.85*ngenes,y=(Ymax-(i)*unitDist),xlim=c(0,ngenes),col = cols[i],cex=1, pch=16)
                                        # Annotated points
          points(x=as.numeric(annotated),y=as.numeric(annotatedZ),xlim=c(0,ngenes),col = cols[i],cex=1.95, pch=2)
        }
        if(dev == "jpg" | dev == "both"){
                                        # plot annotation and words as text
          points(x=as.numeric(les),y=as.numeric(zs),xlim=c(0,ngenes),col = cols[i], pch=16)
          points(x=.85*ngenes,y=(Ymax-(i)*unitDist),xlim=c(0,ngenes),col = cols[i], pch=16)
          points(x=as.numeric(annotated),y=as.numeric(annotatedZ),xlim=c(0,ngenes),col = cols[i], pch=2,lwd=5)
        }
        if(dev == "pdf" | dev == "both"){
                                        # plot annotation and words as text
          points(x=as.numeric(les),y=as.numeric(zs),xlim=c(0,ngenes),col = cols[i], pch=16)
          points(x=.85*ngenes,y=(Ymax-(i)*unitDist),xlim=c(0,ngenes),col = cols[i], pch=16)
          points(x=as.numeric(annotated),y=as.numeric(annotatedZ),xlim=c(0,ngenes),col = cols[i], pch=2,lwd=2)
        }
        text(x=0.85*ngenes,y=(Ymax-(i)*unitDist),pos=4,labels=paste("Motif ",i,"Count",mm))
        text(x=0.85*ngenes,y=(Ymax-(i + 0.5)*unitDist),pos=4,labels=paste("Top word ",mfs[1]))
      }
    }
                                        # Scheme to avoid that words are not plottet on top of exch other (approximation very difficult problem)
    oldX4 = 0
    oldY4 = 0
    oldX3 = 0
    oldY3 = 0
    oldX2 = 0
    oldY2 = 0
    oldX = 0
    oldY = 0
    noPlot = 0
    Xint = ngenes/max(Xdiv-as.numeric(as.numeric(strsplit(len[m],",")[[1]])))
    if(length(dat[,3]) < textWords){
      textWords = length(dat[,3])
    }   
    for(j in 1:textWords){
      if((abs(as.numeric(dat[j,3]) - oldX) < Xint) && abs(as.numeric(dat[j,2]) - oldY) < (2/5.0)*unitDist){
        noPlot = noPlot + 1
      }else{
        if((abs(as.numeric(dat[j,3]) - oldX2) < Xint) && abs(as.numeric(dat[j,2]) - oldY2) < (2.0/5.0)*unitDist){
          noPlot = noPlot + 1
        }else{
          if((abs(as.numeric(dat[j,3]) - oldX3) < Xint) && abs(as.numeric(dat[j,2]) - oldY3) < (2.0/5.0)*unitDist){
            noPlot = noPlot + 1
          }else{
            if((abs(as.numeric(dat[j,3]) - oldX4) < Xint) && abs(as.numeric(dat[j,2]) - oldY4) < (2.0/5.0)*unitDist){
              noPlot = noPlot + 1
            }else{
              text(x=as.numeric(dat[j,3]),y=as.numeric(dat[j,2]),pos=4,labels=dat[j,1])
              oldX4 = oldX3
              oldY4 = oldY3
              oldX3 = oldX2
              oldY3 = oldY2
              oldX2 = oldX
              oldY2 = oldY
              oldX = as.numeric(dat[j,3])
              oldY = as.numeric(dat[j,2])
            }
          }
        }
      }
    }
                                        #Only mark seed site given as input
    for(s in sites(seed,len[m])){
                                        # Marks all annotated seed sites
                                        #for(s in as.character(annotTable[,1]))
      if(dev == "svg"){
        points(x=as.numeric(dat[which(dat[,1]==s),c(3,2)])[1],y=as.numeric(dat[which(dat[,1]==s),c(3,2)])[2],xlim=c(0,ngenes),main = paste(title,len[m],sep=""),col="black", pch=2,cex=2)
      }
      if(dev == "jpg" | dev == "both"){
        points(x=as.numeric(dat[which(dat[,1]==s),c(3,2)])[1],y=as.numeric(dat[which(dat[,1]==s),c(3,2)])[2],xlim=c(0,ngenes),main = paste(title,len[m],sep=""),col="black", pch=2,cex=1,lwd=5)
      }
      if(dev == "pdf" | dev == "both"){
        points(x=as.numeric(dat[which(dat[,1]==s),c(3,2)])[1],y=as.numeric(dat[which(dat[,1]==s),c(3,2)])[2],xlim=c(0,ngenes),main = paste(title,len[m],sep=""),col="black", pch=2,cex=1,lwd=2)
      }
    }
    dev.off()
                                        # plot histograms to show how many words have same leading edge
    if(his){
      jpeg(file=paste(path,"hist",fileName,".",len[m],".jpg",sep=""),height=800,width=1000,quality=90)
      hist(as.numeric(dat[,3]),xlim=c(1,ngenes),breaks=seq(0,ngenes,20),main=paste(title,len[m],sep=""),col=rgb(0,100,130,100,maxColorValue=255), pch=16)
      dev.off()
    }
  } #end for(m in 1:length(len)){
} #end ledge()

                                        # Get arguments
args = tail(commandArgs(),8)
print(args)
file = args[1]
set  = c("positive","negative")
seed = args[2]
len = args[3]
ngenes = args[4]
Ngens = as.numeric(strsplit(ngenes,",")[[1]])
print(paste("Number of Genes",Ngens))
base_path = args[5]
path = args[6]
species = args[7]
annotations = args[8]
dev = c("jpg","pdf")
for(j in 1:length(dev)){
  for (i in 1:2){
    ledge_plot(file, set[i], seed, len,Ngens[i],path,base_path,species,annotations,dev[j])
  }
}
#print("Done")
# example call
# Rscript-2.12.2 lib/mkplot.R Linsley2007-miR-16-KO.rnk "" 6,7,8 1590,1590 "" /data/simras/methods/ResultMatrix_Gapmer_MA08003_new.txt_j4npy/  human T
