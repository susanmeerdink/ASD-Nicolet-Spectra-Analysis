### Nicolet_Data_Prep ###
## Susan Meerdink
## 10/12/2016
## This code reads in Nicolet spectra that have been averaged using 
## https://github.com/susanmeerdink/ASD-Nicolet-Spectra-Processing.
## The result of this Java program is a .csv file with a two rows for each
## spectra process: averaged spectra and standard deviation.
## This code will format that spectral library for future analysis (pulls out species with 3 or more samples).
###---------------------------------------------------------------------------------------------------------------------------- ###
### Setting Up Environment ###
library(dplyr) #load dplyr for data manipulation 
library(multcompView) # load library (install if necessary)
library(dunn.test) #load dunn.test library for kruskal wallis post test - dunn 
library(pgirmess) #load pgirmess library for kruskal wallis post hoc test
library(ggplot2) #load ggplot2 for graphics/figures

### Reading in Data ###
directory <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\"
dataFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\Nicolet_Averaged_Spectra_Huntington.csv" #Set to spectra file location
metaFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\Nicolet_Averaged_Spectra_Metadata_Huntington.csv" #Set to metadata file location
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
    if(i != 3 & i != 6 & i != 7 & i != 10 & i != 28){ #,7,10,28
      print(i)
      tempSpectra <- allSpectra[which(allMeta$Acronym %in% replicationCount$Acronym[i]),] #Get data for that species
      tempMeta <- allMeta[which(allMeta$Acronym %in% replicationCount$Acronym[i]),] #Get metadata for that species
      spectra <- rbind(spectra, tempSpectra)
      meta <- rbind(meta,tempMeta)
    }
  }
}
noSpeciesList = c('AGAT','ALARF','ALBA','BERE','FICO')

fileNameKW <- paste(directory,"Nicolet_Averaged_Spectra_Huntington_rep3.csv",sep="") #create filename that will hold kw test results
write.table(spectra,file = fileNameKW,sep=",")#write header to file

fileNameKW <- paste(directory,"Nicolet_Averaged_Spectra_Metadata_Huntington_rep3.csv",sep="")
write.table(meta,file = fileNameKW,sep=",")

## Grouping Data by Species ##
acronym <- unique(meta$Acronym) #Get unique listing of acronyms
avgSpectra <- vector() #empty vector that will hold averaged species spectra (all samples into one spectrum)
for(a in acronym){ #Loop through acronyms
  temp <- data.frame(spectra[which(meta$Acronym %in% a),2:ncol(spectra)])
  avgSpectra <- rbind(avgSpectra,colMeans(temp)) #add to avgSpectra variable
}

fileName <- paste(directory,"Nicolet_Averaged_Spectra_Huntington_species.csv",sep="") #create filename that will hold kw test results
write.table(avgSpectra,file = fileName,sep=",")#write header to file
###-----END----------------------------------------------------------------------------------------------------------------------- ###
