library("devtools")
library(roxygen2)
#https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
setwd(".")
#create("sircle") #add in R scripts into the R folder --> just need to be done once, afterwards for updating not needed again
setwd("./sircle")

document()
setwd("..")

#how to install the package
library("devtools")
library(roxygen2)
install("sircle")
