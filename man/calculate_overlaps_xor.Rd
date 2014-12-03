\name{calculate_overlaps_xor}
\encoding{UTF-8}
\alias{calculate_overlaps_xor}
\alias{calculate_overlaps_xor,character-method}
\docType{methods}
\title{Calculate overlap and XOR images}
\description{Calculate overlap and XOR images for each combination of every subset histogram   image.}
\usage{\S4method{calculate_overlaps_xor}{character}(subsets,verbose=FALSE)}
\arguments{
\item{subsets}{List of the subset image files. All files have to be in one folder. See the manual for more information about creating the list.}
\item{verbose}{logical (default=FALSE). Change to TRUE to print the calculated values to two new files called "Results_overlaps.txt" and "Results_xor.txt" to the working directory.}
}
\details{After saving a list containing the filenames of the subset histogram images this method calculates the XOR and overlap images/values for each combination of every image, returns the values and is able to write the values to two new files called "Results_overlaps.txt" and "Results_xor.txt" (see example section). See reference Koch et al. 2013 for more information about the calculation.}

\value{The \emph{calculate_overlaps_xor()} method calculates the overlap and XOR images and returns a list with two data frames containing the calculated data.}

\references{\emph{Christin Koch, Ingo Fetzer, Hauke Harms, and Susann Müller. CHIC - An Automated
Approach for the Detection of Dynamic Variations in Complex Microbial Communities.
Cytometry Part A, 2013.}}

\author{Joachim Schumann \email{joachim.schumann@ufz.de}, Christin Koch \email{christin.koch@ufz.de}, Ingo Fetzer \email{info.fetzer@stockholmresilience.su.se}, Susann Müller \email{susann.mueller@ufz.de}}

\examples{
require(EBImage)

## Calculate the overlap and XOR images
## Save the returned values as a list

# Get a list of the filenames of the FCS files
files <- list.files(system.file("extdata",package="flowCHIC"),
full=TRUE,pattern="*.fcs")
# Create histogram images and save them
fcs_to_img(files)
# Create subsets
img_sub(files,x_start=200,x_end=3500,y_start=1000,y_end=3000,maxv=160)
# Get a list of the filenames of the PNG files
subsets <- list.files(path=paste(getwd(),"chic_subset",sep="/"),full=TRUE,pattern="*.png")
# Calculate and save values as a list
results<-calculate_overlaps_xor(subsets)

## Calculate the overlap and XOR images
## Two new files called "Results_overlaps.txt" and
## "Results_xor.txt" are written to the working directory

# Get a list of the filenames of the FCS files
files <- list.files(system.file("extdata",package="flowCHIC"),
full=TRUE,pattern="*.fcs")
# Create histogram images and save them
fcs_to_img(files)
# Create subsets
img_sub(files,x_start=200,x_end=3500,y_start=1000,y_end=3000,maxv=160)
# Get a list of the filenames of the PNG files
subsets <- list.files(path=paste(getwd(),"chic_subset",sep="/"),full=TRUE,pattern="*.png")
# Calculate
calculate_overlaps_xor(subsets,verbose=TRUE)
}

\keyword{overlaps,xor}
\keyword{file.show(file.path(R.home("doc"),"KEYWORDS"))}