\name{img_sub}
\encoding{UTF-8}
\alias{img_sub}
\alias{img_sub,character-method}
\docType{methods}
\title{Create image subsets}
\description{Create subsets of FCS files and the resulting histogram images.}
\usage{\S4method{img_sub}{character}(files,transformation=FALSE,ch1="FS.Log",ch2="FL.4.Log",x_start=0,x_end=4095,y_start=0,y_end=4095,xbin=128,maxv=200,width=300,height=300,\dots)}
\arguments{
\item{files}{Character list of .fcs files. All files have to be in one folder. See the manual for more information about creating the list.}
\item{transformation}{Character string to define the type of data transformation (default=\emph{FALSE}). Fore more details type "?read.FCS" into R command line.}
\item{ch1}{Character string indicating the first channel of the histogram (x-axis) (default=\emph{"FS.Log"}). See the manual for more details.}
\item{ch2}{Character string indicating the second channel of the histogram (y-axis) (default=\emph{"FL.4.Log"}). See the manual for more details.}
\item{x_start}{Start of the rectangle gate on the x-axis (default=\emph{0}). See the manual for more details.}
\item{x_end}{End of the rectangle gate on the x-axis (default=\emph{4095}). See the manual for more details.}
\item{y_start}{Start of the rectangle gate on the y-axis (default=\emph{0}). See the manual for more details.}
\item{y_end}{End of the rectangle gate on the y-axis (default=\emph{4095}). See the manual for more details.}
\item{xbin}{Number of bins within the histogram (default=\emph{128}).}
\item{maxv}{Maximal value of the expressions within the histogram that is set to the highest color value (black) (default=\emph{200}). See the manual for more details.}
\item{width}{Width (pixel) of the resulting histogram image (default=\emph{300}).}
\item{height}{Height(pixel) of the resulting histogram image (default=\emph{300}).}
\item{...}{Additional parameters used for reading the FCS files, creating the PNG images and for plotting. For more details type "?read.FCS", "?png" or "??ggplot2" into R command line.}
}
\details{This method creates subsets of FCS files and the resulting histogram images using the \emph{ggplot} method of the package "ggplot2" (see reference Wickham 2009). After creating a list containing the names of the FCS files a new folder called "chic_subset" is created in the working directory  that contains the subset histogram images. Choose the two channels that are used for plotting on the x/y-axis. Define the start and the end of the rectangle gate of both axes. See the manual for more details.}

\value{The \emph{img_sub()} method creates subsets of the histogram images.}

\references{\emph{H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York,2009.}}

\author{Joachim Schumann \email{joachim.schumann@ufz.de}, Christin Koch \email{christin.koch@ufz.de}, Ingo Fetzer \email{info.fetzer@stockholmresilience.su.se}, Susann Müller \email{susann.mueller@ufz.de}}

\examples{
require(flowCore)
require(ggplot2)

## Write the subset histogram images of the FCS files that are included 
## to the package in a new subfolder of the working directory
## called "chic_subset"

# Get a list of the filenames of the FCS files
files <- list.files(system.file("extdata",package="flowCHIC"),
full=TRUE,pattern="*.fcs")
# Create subsets
img_sub(files,x_start=200,x_end=3500,y_start=1000,y_end=3000,maxv=160)
}

\keyword{subset,fcs}
\keyword{file.show(file.path(R.home("doc"),"KEYWORDS"))}