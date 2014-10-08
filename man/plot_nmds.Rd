\name{plot_nmds}
\encoding{UTF-8}
\alias{plot_nmds}
\alias{plot_nmds,data.frame,data.frame-method}
\docType{methods}
\title{NMDS plot of samples}
\description{NMDS plot of samples based on the calculated XOR and overlap values}   
\usage{\S4method{plot_nmds}{data.frame,data.frame}(x,y,show_cluster=FALSE,type="p",main="",col_nmds="black",cex=0.6,pos=1,
group,legend_pos="topleft",abiotic,p.max=0.05,col_abiotic="magenta",verbose=FALSE,\dots)}
\arguments{
\item{x}{Table with calculated overlap data.}
\item{y}{Table with calculated XOR data.}
\item{show_cluster}{logical (default=\emph{FALSE}). Change to TRUE if cluster dendrogram shall be plotted.}
\item{type}{Type of the plot (default=\emph{"p"}). The \emph{"p"} indicates points without connecting lines within the plot. 
Only used if \bold{group}=FALSE. For more details type "?points" into R command line.}
\item{main}{Character string used as title of the NMDS plot (default="").}
\item{col_nmds}{Color used for the plotted data points if \bold{group}=FALSE (default=\emph{"black"}).}
\item{cex}{numeric (default=0.6). Character expansion factor. Used for the final size of the characters.}
\item{pos}{Position of the text (default=2). Values of 1, 2, 3 and 4, respectively indicate positions below, 
to the left of, above and to the right of the specified coordinates.}
\item{group}{Data frame containing group assignments. The order and the number of these groups has to be identical to the order and the number of the samples printed in R. Use only integer values in the range from 1 to 25. See the manual for more details.}
\item{legend_pos}{Position of the legend (default=\emph{"topleft"}) if \bold{group}=TRUE. 
For more details type "?legend" into R command line.}
\item{abiotic}{Table with abiotic data. Should be a tab-delimited text file using '.' as decimal delimiter. 
Use one row for one sample and one column for one abiotic or experimental parameter. Use the first column for the first parameter 
and the first line as header. The order and the number of the lines has to be identical to the order and the number of the 
samples printed in R.}
\item{p.max}{Decimal number defining the significance level of the abiotic parameters (default=0.05) if \bold{abiotic}=TRUE. Only parameters less/equal this value are plotted.}
\item{col_abiotic}{Color used for the plotted abiotic parameters (default=\emph{"magenta"}).}
\item{verbose}{logical (default=FALSE). Do not print additional information. Change to TRUE to print results of the 
metaMDS method and the p-values of the abiotic parameters.}
\item{...}{Additional parameters used for plotting the data points if \bold{group}=FALSE. 
For more details type "?points" into R command line.}
}
\details{
This method is used for calculating the similarities found in the histogram images of cytometric data. A dissimilarity matrix is generated from the pairwise comparison of histogram images based on the values returned by the method \emph{calculate_overlaps_xor} or saved in the files "Results_overlaps.txt" and "Results_xor.txt". See reference Koch et al. 2013 for more details. Ensuing from this matrix nonmetric multidimensional scaling (NMDS) is performed to show the results. The NMDS plot is calculated using the \emph{metaMDS} method of the package "vegan" (see reference Warnes et al. 2013). In addition, a cluster analysis can be performed to reveal samples with high similarities.
}

\value{The \emph{plot_nmds()} method calculates a NMDS plot of the samples and an additional cluster dendrogram.}

\references{\emph{Christin Koch, Ingo Fetzer, Hauke Harms, and Susann Müller. CHIC - An Automated
Approach for the Detection of Dynamic Variations in Complex Microbial Communities.
Cytometry Part A, 2013.}

\emph{Jari Oksanen, F. Guillaume Blanchet, Roeland Kindt, Pierre Legendre, Peter R. Minchin, R. B.
  O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens and Helene Wagner (2013). vegan:
  Community Ecology Package. R package version 2.0-10. http://CRAN.R-project.org/package=vegan}}

\author{Joachim Schumann \email{joachim.schumann@ufz.de}, Christin Koch \email{christin.koch@ufz.de}, Ingo Fetzer \email{ingo.fetzer@stockholmresilience.su.se}, 
Susann Müller \email{susann.mueller@ufz.de}}
\examples{
require(flowCore)
require(ggplot2)
require(vegan)

## Show the NMDS plot of the sample files
## included to the package

# Get a list of the filenames of the FCS files
files <- list.files(system.file("extdata",package="flowCHIC"),
full=TRUE,pattern="*.fcs")
# Create histogram images
fcs_to_img(files)
# Create subsets
img_sub(files,x_start=200,x_end=3500,y_start=1000,y_end=3000,maxv=160)
# Get a list of the filenames of the subset PNG files
subsets <- list.files(path=paste(getwd(),"chic_subset",sep="/"),full=TRUE,pattern="*.png")
# Calculate
results<-calculate_overlaps_xor(subsets)
# Show NMDS plot ensuing from the returned values
plot_nmds(x=results$overlap,y=results$xor)


## Show the NMDS plot of the dataset "mix""

data(Results_overlaps_mix)
data(Results_xor_mix)
plot_nmds(Results_overlaps_mix,Results_xor_mix)

## Show the NMDS plot of the dataset "incol""

data(Results_overlaps_incol)
data(Results_xor_incol)
plot_nmds(Results_overlaps_incol,Results_xor_incol)

## Show the NMDS plot of the dataset "incol"" with group assignment
## and abiotic data
## Print results of the metaMDS method and the p-values of the abiotic parameters
## Show a cluster dendrogram

data(Results_overlaps_incol)
data(Results_xor_incol)
groups<-data.frame("groups"=c(15,19,19,19,15,22,19,15,22,15,15,22,22,22,22,19,19))
data(abiotic_incol)
plot_nmds(Results_overlaps_incol,Results_xor_incol,show_cluster=TRUE,group=groups,abiotic=abiotic_incol[,-1],verbose=TRUE)

}

\keyword{nmds}
\keyword{file.show(file.path(R.home("doc"),"KEYWORDS"))}
