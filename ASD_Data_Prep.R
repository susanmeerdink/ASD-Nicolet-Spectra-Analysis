### ASD_Data_Prep ###
## Susan Meerdink
## 10/1/2016
## This code reads in ASD spectra that have been averaged using 
## https://github.com/susanmeerdink/ASD-Nicolet-Spectra-Processing.
## The result of this Java program is a .csv file with a two rows for each
## spectra process: averaged spectra and standard deviation.
## This code will format that spectral library for future analysis (pulls out species with 3 or more samples).
###---------------------------------------------------------------------------------------------------------------------------- ###
### Setting Up Environment ###


### Reading in Data ###
directory <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\"
dataFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\ASD_Averaged_Spectra_Huntington.csv" #Set to spectra file location
metaFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\ASD_Averaged_Spectra_Metadata_Huntington.csv" #Set to metadata file location
data <- read.csv(dataFile) #Read in the averaged and std of spectra
allMeta <- read.csv(metaFile) #Read in associated metadata of spectra

### Setting up Variables ###
allSpectra <- data[which(data[,2] %in% "AVG"),] #Get averaged spectra
allSpectra$Code <- NULL #Remove column 2 with AVG values
wavelength <- as.numeric(gsub('X','',names(data)[3:length(names(data))])) #Get wavelengths

### Pulling out Samples with Replicates
##Some samples only have n = 1, ignore those for statistical analysis
replicationCount <- dplyr::count(allMeta,Acronym)
spectra <- vector()
meta <- vector()
for(i in 1:nrow(replicationCount)){
  if(replicationCount$n[i] > 2){ #if there is more than one sample for a species...
    tempSpectra <- allSpectra[which(allMeta$Acronym %in% replicationCount$Acronym[i]),] #Get data for that species
    tempMeta <- allMeta[which(allMeta$Acronym %in% replicationCount$Acronym[i]),] #Get metadata for that species
    spectra <- rbind(spectra, tempSpectra)
    meta <- rbind(meta,tempMeta)
  }
}
fileNameKW <- paste(directory,"ASD_Averaged_Spectra_Huntington_rep3.csv",sep="") #create filename that will hold kw test results
write.table(spectra,file = fileNameKW,sep=",")#write header to file

fileNameKW <- paste(directory,"ASD_Averaged_Spectra_Metadata_Huntington_rep3.csv",sep="")
write.table(meta,file = fileNameKW,sep=",")

#Does require some tweaking of the files - so check them before continuing on to the next step
###-----END----------------------------------------------------------------------------------------------------------------------- ###
