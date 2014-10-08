######
## Create NMDS plot
######

setMethod("plot_nmds",signature=c(x="data.frame",y="data.frame"),
definition = function(x,y, show_cluster=FALSE,type="p",main="",
                      col_nmds="black",cex=0.6,pos=1,group,legend_pos="topleft",abiotic,
                      p.max=0.05,col_abiotic="magenta",verbose=FALSE,...){

  newlist<-list()
  
  x$Label<-gsub("overlap_","",x$Label)  #removes "overlap_" from file labels
  x$Label<-gsub(".png","",x$Label)  #removes ".png" from file labels
  
  label.split<-strsplit(as.character(x$Label),"_&_")  #creates a file with proper label names
  
  
  nam<-NULL  #creates an empty vector file
  for (i in 1:length(label.split)) nam<-rbind(nam, cbind(label.split[[i]][1],
                                                         label.split[[i]][2]))  #creates a file ("nam") with all combinations of label names in it 
  
  dats<-data.frame(Label=x$Label, Label.1=nam[,1], Label.2=nam[,2],
                   IntDen=y$IntDen, Area=x$Area, IntDen.Area=y$IntDen/x$Area/100)  #final data file ("dats") is created by merging and conversion of XOR and overlay data
  
  size<-length(unique(dats$Label.1))+1  #number of unique data labels is estimated
  mat<-matrix(nrow=size, ncol=size)  #empty matrix ('mat') for disimilariy matrix is created
  j=1  #counter J created and set to 1
  for (i in unique(dats$Label.1)) {  #dissimilarity according to formula 1 in MS for each combination is calculated and stored in matrix 'mat'
    mat[,j]<-c(rep(NA,size-length(dats$IntDen.Area[dats$Label.1==i])),
               dats$IntDen.Area[dats$Label.1==i])
    j<-j+1
  }
  
  diag(mat)<-0  #diagonals are filled up with '0' values
  colnames(mat)<-rownames(mat)<-union(unique(dats$Label.1),
                                      unique(dats$Label.2))  #label names are written into matrix
  print(colnames(mat))
  mat.dist<-as.dist(mat)  #matrix mat is converted to dissimilarity matrix ('mat.dist') to be used with metaMDS
  
  if(show_cluster){  
    par(mfrow=c(1,2))
    plot(hclust(mat.dist, "average"), hang=-1)  #calculates and plots cluster analysis from the dissimilarity matrix using 'average' distances within cluster grouping
  }
  
  mds.out<-metaMDS(mat.dist, autotransform=FALSE)  #NMDS analysis is performed using dissimilarity matrix 'mat.dist'
  newlist$mds.out<-mds.out
  
  
  if(missing(group)){    
    plot(mds.out, type="n",display="sites",main=main)  ## plots result
    points(mds.out, type=type,col=col_nmds,...)  ## add points
    text(mds.out,cex=cex,pos=pos) 
    
  
  } else{
    color<-group[,1]
    plot(mds.out, type="n",display="sites",main=main)
    points(mds.out, col=color , pch=as.numeric(color))
    text(mds.out,cex=cex,pos=pos) 
    legend(x=legend_pos, legend=unique(color),  pch=as.numeric(unique(color)),
           col=as.numeric(unique(color)))
    ordihull(mds.out, groups=color, lty=2, col="darkgrey")
  }
  
  #######The following section gives a description if you want to include environmental data for MDS analysis###################

  if(!missing(abiotic)){
    ef<-envfit(mds.out, abiotic, permutation=999)  #calculates the relevant environmental parameters for the nMDS #result based on 999 Monte-Carlo permutations
    plot(ef,cex=cex,p.max=p.max, col=col_abiotic)  #plots the most relevant environmental parameters (with a p.max <=0.05) into the MDS
    newlist$ef<-ef
  }
  
  
  if(verbose){
    return(newlist)
  }
 
})