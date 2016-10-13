### Nicolet-Spectra-Analysis ###
## Susan Meerdink
## 8/2/2016
## This code reads in Nicolet spectra that have been averaged using 
## https://github.com/susanmeerdink/ASD-Nicolet-Spectra-Processing.
## The result of this Java program is a .csv file with a two rows for each
## spectra process: averaged spectra and standard deviation.
## This code will analyze the spectra to determine similarities and differences.
###---------------------------------------------------------------------------------------------------------------------------- ###
### Setting Up Environment ###
library(ggplot2) #load ggplot2 library for figures
library(moments) #load moments library
library(reshape2) #load reshape2 library for data manipulation
library(pgirmess) #load pgirmess library for kruskal wallis post hoc test
library(dplyr) #load dplyr for data manipulation 

### Reading in Data ###
directory <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\"
dataFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\Nicolet_Averaged_Spectra_Huntington.csv" #Set to spectra file location
metaFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\Nicolet_Averaged_Spectra_Metadata_Huntington.csv" #Set to metadata file location
data <- read.csv(dataFile) #Read in the averaged and std of spectra
allMeta <- read.csv(metaFile) #Read in associated metadata of spectra

### Setting up Variables ###
allSpectra <- data[which(data[,2] %in% "AVG"),] #Get averaged spectra
std <- data[which(data[,2] %in% "STD"),] #Get std of spectra
allSpectra$um <- NULL #Remove column 2 with AVG values
std$um <- NULL #Remove column 2 with STD values
names(allSpectra) <- sub("Wavelength","ID",names(allSpectra)) #Rename first column to ID instead of Wavelength
names(std) <- sub("Wavelength","ID",names(std)) #Rename first column to ID instead of Wavelength
wavelength <- as.numeric(gsub('X','',names(data)[3:length(names(data))])) #Get wavelengths
###---------------------------------------------------------------------------------------------------------------------------- ###

### Plotting Spectra  by species ###
## WARNING: will produce 41 image files!
acronym <- unique(allMeta$Acronym) #Get unique listing of acronyms
for(a in acronym){
  plotData <- data.frame(wl = wavelength, t(allSpectra[which(allMeta$Acronym %in% a),2:ncol(allSpectra)])) #make current data into data.frame for plotting
  colnames(plotData) <- c("wl",as.character(allSpectra[which(allMeta$Acronym %in% a),1])) #rename columns with ID
  plotDataMelt <-  melt(plotData,id.vars = 'wl') #Melt using wavelength as the ID variable (puts it in format to plot)
  
  ggplot(data = plotDataMelt, aes(x = wl,y = value, group = variable,color=variable)) + #Plotting the spectra
    geom_line()+ #plot as a line
    xlab(expression(paste("Wavelength (",mu,"m)",sep ="")))+ #add xlabel
    coord_cartesian(xlim = c(2.5,13), ylim = c(0,15)) + #Fix axis
    ylab("Reflectance (%)")+ #add ylabel
    ggtitle(a) #add title
  
  plotName = paste(directory,"Spectra_Plots\\",a,"_spectra.png",sep="") #create plot name to save the file
  ggsave(plotName,width = 10,height = 5) #save the plot
}
###---------------------------------------------------------------------------------------------------------------------------- ###

### Plotting Spectra  by Genus ###
## WARNING: will produce many image files!
genus <- unique(allMeta$Genus) #Get unique listing of acronyms
for(a in genus){
  plotData <- data.frame(wl = wavelength, t(allSpectra[which(allMeta$Genus %in% a),2:ncol(allSpectra)])) #make current data into data.frame for plotting
  colnames(plotData) <- c("wl",as.character(allSpectra[which(allMeta$Genus %in% a),1])) #rename columns with ID
  plotDataMelt <-  melt(plotData,id.vars = 'wl') #Melt using wavelength as the ID variable (puts it in format to plot)
  
  ggplot(data = plotDataMelt, aes(x = wl,y = value, group = variable,color=variable)) + #Plotting the spectra
    geom_line()+ #plot as a line
    xlab(expression(paste("Wavelength (",mu,"m)",sep ="")))+ #add xlabel
    coord_cartesian(xlim = c(2.5,13)) + #Fix axis
    ylab("Reflectance (%)")+ #add ylabel
    ggtitle(a) #add title
  
  plotName = paste(directory,"Genus_Plots\\",a,"_spectra.png",sep="") #create plot name to save the file
  ggsave(plotName,width = 10,height = 5) #save the plot
}
###---------------------------------------------------------------------------------------------------------------------------- ###

### Plotting Average & STD of Spectra ###
totalSpectra <- array(0,dim = c(length(wavelength),3))
for(x in c(2:ncol(allSpectra))){ #Loop through wavelengths 
  totalSpectra[x-1,1] <- mean(allSpectra[,(x)])
  totalSpectra[x-1,2] <- min(allSpectra[,(x)]) 
  totalSpectra[x-1,3] <- max(allSpectra[,(x)])
}
plotData <- data.frame(wl = wavelength,totalSpectra)
colnames(plotData) <- c("wl","avg","min","max") #rename columns with ID
plotDataMelt <-  melt(plotData,id.vars = 'wl') #Melt using wavelength as the ID variable (puts it in format to plot)
ggplot(plotDataMelt,aes(x = wl, y = value, group = variable, color= variable)) + #Plotting the spectra  
  geom_line()+ #plot as a line
  geom_ribbon(aes(ymin = value,ymax = value)) +
  xlab(expression(paste("Wavelength (",mu,"m)",sep ="")))+ #add xlabel
  ylab("Reflectance (%)") + #add ylabel
  coord_cartesian(xlim = c(2.5,13), ylim = c(0,25)) +#Fix axis
  scale_y_continuous(expand = c(0,0)) + #remove y axis buffer (starts at 0% now)
  scale_x_continuous(expand = c(0,0)) #remove x axis buffer (starts at 0% now)
plotName = paste(directory,"Spectra_Plots\\","AVG_STD_spectra.png",sep="") #create plot name to save the file
ggsave(plotName,width = 10,height = 5) #save the plot

###---------------------------------------------------------------------------------------------------------------------------- ###

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
spectra$Code <- NULL #remove code column

###---------------------------------------------------------------------------------------------------------------------------- ###

### Testing for Normality ###
## WARNING: will produce 1738 image files!
normTest <- array(0,dim = c(length(wavelength),5))#Empty array that will hold values below
normTest[,1] <- wavelength #add wavelength to the first column
for(x in c(2:ncol(spectra))){#Loop through wavelengths
  
  test1 <- shapiro.test(spectra[,x]) ##Shapiro-Wilk Test for normality
  normTest[x-1,3] <- test1$p.value #Add p-value to table
  if(test1$p.value > 0.05){ #retain the null hypothesis, normally distributed
    normTest[x-1,2] <- 1 ## 1 = normally distributed
  }else{ #reject the null hypothesis, NOT normally distributed
    normTest[x-1,2] <- 0 ## 0 = not normally distributed
  }
  
  normTest[x-1,4] <- skewness(spectra[,x])##Skewness (asymmetric)
  
  normTest[x-1,5] <- kurtosis(spectra[,x])##Kurosis (Pointed)
  
  plotData <- data.frame(spectra[,x]) #add data to data frame for plotting
  colnames(plotData) <- "wl" #change column name for plotting
  ggplot(data = plotData,aes(wl))+ ##Plotting Histogram
    geom_histogram(binwidth = 1) + #plot as histogram
    ggtitle(wavelength[x-1]) + #add title to histogram
    xlab("Reflectance (%)") #add xlabel
  
  plotName <- paste(directory,"Normality_Test\\",round(wavelength[x-1]*1000,digits = 0),"_hist.png",sep="") #create plot name to save the file
  ggsave(plotName,width = 3,height = 3) #save the plot
}
fileName <- paste(directory,"Normality_Test\\normality_test_results.csv",sep="") 
write(t(normTest),file = fileName,sep = ",")#save normality test data to file
##This data is not Normal - using non-parametric tests
###---------------------------------------------------------------------------------------------------------------------------- ###

### Kruskal-Wallis ###
##Since data is non-parametric, Kruskal-Wallis is need to compare samples
kwResults <- array(0,dim = c(length(wavelength),2))#Empty array that will hold values below
kwResults[,1] <- wavelength #add wavelength to the first column

fileNameKW <- paste(directory,"Kruskal_Wallis_Test\\kruskal_wallis_test_results.csv",sep="") #create filename that will hold kw test results
write("Wavelength,KW_pvalue,flag",file = fileNameKW)#write header to file
fileNameKWPH <- paste(directory,"Kruskal_Wallis_Test\\kruskal_wallis_posthoc_test_results.csv",sep="") #create filename that will hold kw post hoc test results
write("Wavlength,Pair,P-Value,Obs-Diff,Critical-Diff",file = fileNameKWPH)

for(x in as.numeric(c(2,3,4))){ #Loop through wavelengths c(2:length(wavelength)) as.numeric(c(2,3,4))
  dataKW <- data.frame(acronymIN = meta$Acronym,spectraIN = spectra[,x])
  
  kwTest <- kruskal.test(spectraIN~acronymIN, data = dataKW) #perform kruskal-wallis test
  kwResults[x-1,2] <- kwTest$p.value #add pvalue to second column
  if(kwTest$p.value < 0.05){#If the kruskal-wallis test yields a pvalue less than 0.05, set flag to TRUE
    flag = 'TRUE' #Values differ in at least one    
  } else{#If the kruskal-wallis test yields a pvalue less than 0.05, set flag to FALSE
    flag = 'FALSE'#values do not differ in location distribution
  }
  output <- paste(wavelength[x],kwTest$p.value,flag,sep = ",") #Set output for file
  write(output,file = fileNameKW, append=TRUE) #add line to file  

  if(kwTest$p.value < 0.05){ #If the kruskal-wallis test yields a pvalue less than 0.05, run post hoc test to determine which samples differ
    kwPostTest <- kruskalmc(spectraIN~acronymIN, data = dataKW, probs = 0.05) #run post hoc kruskal wallis test 
    test <- kwPostTest$dif.com$difference # select logical vector
    names(test) <- row.names(kwPostTest$dif.com)# add comparison names
#     for(i in c(1:length(kwPostTest$dif.com[,3]))){
#       if(grepl('TRUE',kwPostTest$dif.com[i,3]) == TRUE){
#         output <- c(wavelength[x-1],kwPostTest$signif.level,kwPostTest$dif.com$obs.dif[i],kwPostTest$dif.com$critical.dif[i])
#         write(output,file = fileNameKWPH,append=TRUE,sep = ",") #add line to file      
#       }  
#     }
  }
      
#  ggplot(data = dataKW,(aes(x = acronym,y=spectra)))+ 
#    geom_boxplot() +
#    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + #rotate x labels
#    xlab("Plant Species Acronym") + #add x label
#    ylab("Reflectance (%)") + #add y label
#    ggtitle(wavelength[x-1]) #add title
    
#  plotName <- paste(directory,"Box_Plots\\",round(wavelength[x-1]*1000,digits = 0),"_box.png",sep="") #create plot name to save the file
#  ggsave(plotName,width = 10,height = 5) #save the plot
}

###---------------------------------------------------------------------------------------------------------------------------- ###
