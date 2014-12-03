######
## Calculate overlap and XOR images
######

setMethod("calculate_overlaps_xor",signature=(subsets="character"),

definition = function (subsets,verbose=FALSE) {
  
  if(verbose){
  cat(paste("Label","Area",sep="\t"),file="Results_overlaps.txt",sep="\n") # create table
  cat(paste("Label","IntDen",sep="\t"),file="Results_xor.txt",sep="\n") # create table
  }
  
  newlist<-list() 
  
  inverse = function(x) (x-1)*(-1) # inverse an image
  
  bXOR = function(x, y) {         # perform XOR on two images
    x = round(x@.Data * (2^8-1),digits=3)
    y = round(y@.Data * (2^8-1),digits=3)
    x = x * (2^8-1)
    y = y * (2^8-1)
    xo = bitwXor(x, y)
    Image(xo / (2^8-1), dim = dim(x))
  }
  df_overlap = data.frame(Label=character(),Area=numeric())
  df_xor = data.frame(Label=character(),IntDen=numeric())
  
  for(i in 1:(length(subsets)-1)){
    c=i+1
    x <- readImage(subsets[i])
    xgray<-channel(x,"gray")
    neg1<-inverse(xgray)
    for(j in c:length(subsets)) {
      y <- readImage(subsets[j])
      ygray<-channel(y,"gray")
      neg2<-inverse(ygray)
      addim<-(neg1+neg2)
      addimgr<-channel(addim,"gray")
      bin<-addimgr<.004
      l<-length(bin[bin==FALSE]) # calculate overlap area
      
      overname<-paste("overlap_",sub("^([^.]*).*", "\\1",basename(subsets[i])),
                      "_&_",sub("^([^.]*).*", "\\1",basename(subsets[j])),".png",sep="")
      xo = bXOR(xgray, ygray) # perform XOR on two images
      lxor<-length(xo)
      m<-round(mean(xo)*255,digits=3)
      
      xorname<-paste("xor_",sub("^([^.]*).*", "\\1",basename(subsets[i])),
                     "_&_",sub("^([^.]*).*", "\\1",basename(subsets[j])) ,".png",sep="")
      
      if(verbose){
      cat(paste(overname,l,sep="\t"),file="Results_overlaps.txt",append=TRUE,sep="\n") # write values to file
      cat(paste(xorname,round(lxor*m,digits=0),sep="\t"),file="Results_xor.txt",
          append=TRUE,sep="\n") # write values to file
      }
      
      df_overlap = rbind(df_overlap,data.frame(Label=overname,Area=l))
      df_xor = rbind(df_xor,data.frame(Label=xorname,IntDen=round(lxor*m,digits=0)))	    
      
    }
  }
  newlist$overlap<-df_overlap
  newlist$xor<-df_xor
  return(newlist)
})
