######
## Create histogram images
######

setMethod("fcs_to_img",signature=(files="character"),

definition = function(files,transformation=FALSE,ch1="FS.Log",ch2="FL.4.Log",
                      width=300,height=300,...){
  
  
  dir.create("chic_images")
  setwd("chic_images")
  
  for(f in files) {
    frame <- read.FCS(f,alter.names=TRUE,transformation=transformation,...) # read FCS file
    png(filename=sub(pattern="*.fcs",replacement=".png",basename(f)),width=width,height=height,
         units = "px",...)
    
    mat<-exprs(frame) # get the expression values 
    x<- mat[,ch1]
    y<- mat[,ch2]
    bin <- hexbin(x,y, xbin = 128) # perform hexbin
    df <- data.frame(hcell2xy(bin), count = bin@count)
    
    theme_map <- function (base_size = 12, base_family ="") {
                 theme_bw(base_size = base_size, base_family = base_family) %+replace% 
                 theme(legend.position = "none",
                       panel.background=element_blank(),
                       panel.border=element_blank(),
                       panel.margin=unit(0, "lines"),
                       plot.background=element_blank(),...)} # change theme attibutes
    
    print(ggplot(df, aes(x = x, y = y)) +  geom_point(aes(color = count)) + 
    scale_color_gradient(low = "white", high = "black" , na.value = "black", 
    limit = c(0, 200))+theme_map()+xlab(ch1)+ylab(ch2),...) # plot
   
  
    dev.off()
  } 
  setwd("..")
})