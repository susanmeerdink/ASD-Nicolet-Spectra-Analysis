### Nicolet-Spectra-Analysis ###
## Susan Meerdink
## 8/2/2016
## This code reads in Nicolet spectra that have been averaged using 
## https://github.com/susanmeerdink/ASD-Nicolet-Spectra-Processing.
## The result of this Java program is a .csv file with a two rows for each
## spectra process: averaged spectra and standard deviation.
## This code will analyze the spectra to determine similarities and differences.
##--------------------------------------------------------------------------------
### Setting Up Environment ###
library(ggplot2) #load ggplot2 library
library(moments) #load moments library

### Reading in Data ###
dataFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\Nicolet_Averaged_Spectra.csv" #Set to spectra file location
metaFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\Nicolet_Averaged_Spectra_Metadata.csv" #Set to metadata file location
data <- read.csv(dataFile) #Read in the averaged and std of spectra
meta <- read.csv(metaFile) #Read in associated metadata of spectra

### Setting up Variables ###
avg <- data[which(data[,2] %in% "AVG"),] #Get averaged spectra
std <- data[which(data[,2] %in% "STD"),] #Get std of spectra
avg$um <- NULL #Remove column 2 with AVG values
std$um <- NULL #Remove column 2 with STD values
names(avg) <- sub("Wavelength","ID",names(avg)) #Rename first column to ID instead of Wavelength
names(std) <- sub("Wavelength","ID",names(std)) #Rename first column to ID instead of Wavelength
wavelength <- as.numeric(gsub('X','',names(data)[3:length(names(data))])) #Get wavelengths

### Plotting Spectra ###
#avgOrder <- avg[order(meta$Acronym),] #Reorder avg spectra based on acronym
acronym <- unique(meta$Acronym) #Get unique listing of acronyms
for(a in acronym){
  plotData <- data.frame(wl = wavelength, t(avg[which(meta$Acronym %in% a),2:ncol(avg)]))
  ggplot(data = plotData, aes(x = wavelength, y = plotData[,2:ncol(plotData)])) +
    geom_line()
}

### Testing for Normality ###
normTest <- array(0,dim = c(length(wavelength),4))#Empty array that will hold values below
for(x in c(2:length(wavelength))){
  ##Shapiro-Wilk Test for normality
  test1 <- shapiro.test(avg[,x])
  if(test1$p.value > 0.05){ #retain the null hypothesis, normally distributed
    normTest[x,1] <- 1 ## 1 = normally distributed
  }else{ #reject the null hypothesis, NOT normally distributed
    normTest[x,1] <- 0 ## 0 = not normally distributed
  }
  normTest[x,2] <- test1$p.value
  
  ##Skewness (asymmetric)
  normTest[x,3] <- skewness(avg[,x])
  
  ##Kurosis (Pointed)
  normTest[x,4] <- kurtosis(avg[,x])
}
##This data is not Normal - using non-parametric tests


