######
## Create subsets of the histogram images
######

setMethod("img_sub",signature=(files="character"),
definition = function(files,transformation=FALSE,ch1="FS.Log",ch2="FL.4.Log",x_start=0,
                      x_end=4095,y_start=0,y_end=4095,xbin=128,maxv=200,width=300,
                      height=300,...){

dir.create("chic_subset") # create new folder
setwd("chic_subset") 

for(f in files) {

frame<-read.FCS(f,alter.names=TRUE,transformation=transformation,...) # read FCS file
png(filename=sub(pattern="*.fcs",replacement=".png",basename(f)),width=width,height=height,
    units = "px",...)
coords <- list(c(x_start,x_end), c(y_start,y_end)) # define subset area
names(coords) <- c(ch1, ch2)

Noise <- rectangleGate(filterId="Noise",  .gate = coords) # filter noise
Noise.subset <- Subset(frame, Noise) # create subset

mat<-exprs(Noise.subset) # get expressions
x<- mat[,ch1]
y<- mat[,ch2]
bin <- hexbin(x,y, xbin = xbin) # perform hexbin

df <- data.frame(hcell2xy(bin), count = bin@count)
theme_map <- function (base_size = 12, base_family ="") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position = "none",
          panel.background=element_blank(),panel.border=element_blank(),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.margin=unit(0, "lines"),plot.background=element_blank(),...)} # change theme

print(ggplot(df, aes(x = x, y = y)) +  geom_point(aes(color = count)) +
        scale_color_gradient(low = "white", high = "black" , na.value = "black",
        limit = c(0, maxv))+theme_map(),...) # save subset images

dev.off()
  }
setwd("..")
})