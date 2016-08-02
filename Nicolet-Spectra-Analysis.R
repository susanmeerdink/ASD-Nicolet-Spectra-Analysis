### Nicolet-Spectra-Analysis ###
## Susan Meerdink
## 8/2/2016
## This code reads in Nicolet spectra that have been averaged using 
## https://github.com/susanmeerdink/ASD-Nicolet-Spectra-Processing.
## The result of this Java program is a .csv file with a two rows for each
## spectra process: averaged spectra and standard deviation.
## This code will analyze the spectra to determine similarities and differences.
##--------------------------------------------------------------------------------
### Reading in Data ###
dataFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\Nicolet_Averaged_Spectra.csv" #Set to spectra file location
metaFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\Nicolet_Averaged_Spectra_Metadata.csv" #Set to metadata file location
data <- read.csv(dataFile) #Read in the averaged and std of spectra
meta <- read.csv(metaFile) #Read in associated metadata of spectra