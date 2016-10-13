### Nicolet-Spectra-Analysis ###
## Susan Meerdink
## 10/12/2016
## This code reads in Nicolet spectra that have been averaged using 
## https://github.com/susanmeerdink/ASD-Nicolet-Spectra-Processing.
## The result of this Java program is a .csv file with a two rows for each
## spectra process: averaged spectra and standard deviation.
## This code will analyze the spectra to determine similarities and differences.
###---------------------------------------------------------------------------------------------------------------------------- ###
### Setting Up Environment ###
library(dplyr) #load dplyr for data manipulation 
library(multcompView) # load library (install if necessary)
library(dunn.test) #load dunn.test library for kruskal wallis post test - dunn 
library(pgirmess) #load pgirmess library for kruskal wallis post hoc test

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
    tempSpectra <- allSpectra[which(allMeta$Acronym %in% replicationCount$Acronym[i]),] #Get data for that species
    tempMeta <- allMeta[which(allMeta$Acronym %in% replicationCount$Acronym[i]),] #Get metadata for that species
    spectra <- rbind(spectra, tempSpectra)
    meta <- rbind(meta,tempMeta)
  }
}
###---------------------------------------------------------------------------------------------------------------------------- ###
### Kruskal-Wallis ###
##Since data is non-parametric, Kruskal-Wallis is need to compare samples
kwResults <- array(0,dim = c(length(wavelength),2))#Empty array that will hold values below
kwResults[,1] <- wavelength #add wavelength to the first column
dunnResults <- array(3,0)#Empty array that will hold values below

fileNameKW <- paste(directory,"Kruskal_Wallis_Test\\kruskal_wallis_test_results.csv",sep="") #create filename that will hold kw test results
write("Wavelength,KW_pvalue",file = fileNameKW)#write header to file
fileNameKWPH <- paste(directory,"Kruskal_Wallis_Test\\kruskal_wallis_posthoc_test_results.csv",sep="") #create filename that will hold kw post hoc test results
write("Wavelength,Pair,Obs-Diff,Critical-Diff,Flag",file = fileNameKWPH)
fileNameKWDT <- paste(directory,"Kruskal_Wallis_Test\\dunn_posthoc_test_results.csv",sep="") #create filename that will hold kw post hoc test results
write("Wavelength,Pair,P-value",file = fileNameKWDT)

for(x in c(2:length(spectra)) ){ #Loop through wavelengths 
  #All Wavelengths = c(2:length(spectra)) 
  #Only Some Wavelengths = as.numeric(c(2,3,4))
  
  
  ## kruskal-wallis test ##
  kwTest <- kruskal.test(spectra[,x], meta$Acronym) #perform kruskal-wallis test
  kwResults[x-1,2] <- kwTest$p.value #add pvalue to second column
  output <- paste(wavelength[x-1],kwTest$p.value,sep = ",") #Set output for file
  write(output,file = fileNameKW, append=TRUE) #add line to file
  
  ## Post Hoc kruskal-wallis test with pirgimess##
  kwPostTest <- kruskalmc(spectra[,x], meta$Acronym, probs = 0.05) #run post hoc kruskal wallis test 
  kwPostTestDiff <- kwPostTest$dif.com[complete.cases(kwPostTest$dif.com),] #Get rid of NaN values
  pairNames <- row.names(kwPostTestDiff) #Get names of pairs
  
  for(i in c(1:nrow(kwPostTestDiff))){
    if(grepl('TRUE',kwPostTestDiff[i,3]) == TRUE){
      output <- c(wavelength[x-1],pairNames[i],kwPostTestDiff[i,1],kwPostTestDiff[i,2],kwPostTestDiff[i,3])
      write(output,file = fileNameKWPH,append=TRUE,sep = ",") #add line to file 
    }
  }

  ## Post Hoc Kruskal-wallis test with dunn ##
  kwDunnTest <- dunn.test(x = spectra[,x], g = meta$Acronym)
  
  for(i in c(1:length(kwDunnTest$P))){
    if(kwDunnTest$P[i] < 0.05){
      dunnResults <- rbind(dunnResults,c(wavelength[x-1],kwDunnTest$comparisons[i],kwDunnTest$P[i]))
      output <- paste(wavelength[x-1],kwDunnTest$comparisons[i],kwDunnTest$P[i],sep = ",")
      write(output,file = fileNameKWDT,append=TRUE) #add line to file 
    }
  }
}