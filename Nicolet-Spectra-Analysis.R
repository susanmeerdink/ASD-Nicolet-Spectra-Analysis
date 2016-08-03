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
library(reshape2) #load reshape2 library

### Reading in Data ###
directory <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\"
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
## WARNING: will produce 46 image files!
acronym <- unique(meta$Acronym) #Get unique listing of acronyms
for(a in acronym){
  plotData <- data.frame(wl = wavelength, t(avg[which(meta$Acronym %in% a),2:ncol(avg)])) #make current data into data.frame for plotting
  colnames(plotData) <- c("wl",as.character(avg[which(meta$Acronym %in% a),1])) #rename columns with ID
  plotDataMelt <-  melt(plotData,id.vars = 'wl') #Melt using wavelength as the ID variable (puts it in format to plot)
  
  ggplot(data = plotDataMelt, aes(x = wl,y = value, group = variable,color=variable)) + #Plotting the spectra
    geom_line()+ #plot as a line
    xlab(expression(paste("Wavelength (",mu,"m)",sep ="")))+ #add xlabel
    ylab("Reflectance (%)")+ #add ylabel
    ggtitle(a) #add title
  
  plotName = paste(directory,"Spectra_Plots\\",a,"_spectra.png",sep="") #create plot name to save the file
  ggsave(plotName,width = 10,height = 5) #save the plot
}

### Testing for Normality ###
## WARNING: will produce 1738 image files!
normTest <- array(0,dim = c(length(wavelength),5))#Empty array that will hold values below
normTest[,1] <- wavelength #add wavelength to the first column
for(x in c(2:length(wavelength))){#Loop through wavelengths
  
  test1 <- shapiro.test(avg[,x]) ##Shapiro-Wilk Test for normality
  normTest[x,3] <- test1$p.value #Add p-value to table
  if(test1$p.value > 0.05){ #retain the null hypothesis, normally distributed
    normTest[x,2] <- 1 ## 1 = normally distributed
  }else{ #reject the null hypothesis, NOT normally distributed
    normTest[x,2] <- 0 ## 0 = not normally distributed
  }
  
  normTest[x,4] <- skewness(avg[,x])##Skewness (asymmetric)
  
  normTest[x,5] <- kurtosis(avg[,x])##Kurosis (Pointed)
  
  plotData <- data.frame(avg[,x]) #add data to data frame for plotting
  colnames(plotData) <- "wl" #change column name for plotting
  ggplot(data = plotData,aes(wl))+ ##Plotting Histogram
    geom_histogram(binwidth = 1) + #plot as histogram
    ggtitle(wavelength[x]) + #add title to histogram
    xlab("Reflectance (%)") #add xlabel
  
  plotName <- paste(directory,"Hist_Plots\\",round(wavelength[x]*1000,digits = 0),"_hist.png",sep="") #create plot name to save the file
  ggsave(plotName,width = 3,height = 3) #save the plot
}
fileName <- paste(directory,"Hist_Plots\\normality_test_results.csv",sep="") 
write(t(normTest),file = fileName,sep = ",")#save normality test data to file
##This data is not Normal - using non-parametric tests

### Kruskal-Wallis ###
##Since data is non-parametric, Kruskal-Wallis is need to compare samples
for(x in c(2:length(wavelength))){ #Loop through wavelengths
  
}