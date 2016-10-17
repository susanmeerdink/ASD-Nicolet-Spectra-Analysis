### Nicolet-Spectra-Analysis ###
## Susan Meerdink
## 10/16/2016
## This code reads in Nicolet spectra processed by Nicolet_Data_Prep.R 
## This code will analyze the spectra to determine similarities and differences.
###---------------------------------------------------------------------------------------------------------------------------- ###
### Setting Up Environment ###
library(dplyr) #load dplyr for data manipulation 
library(dunn.test) #load dunn.test library for kruskal wallis post test - dunn 
library(pgirmess) #load pgirmess library for kruskal wallis post hoc test
library(ggplot2) #load ggplot2 for graphics/figures
library(moments) #load moments library

### Reading in Data ###
directory <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\"
dataFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\HyTES_Spectra_Huntington_33.csv" #Set to spectra file location
metaFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\HyTES_Spectra_Huntington_Metadata_33.csv" #Set to metadata file location
spectra <- read.csv(dataFile) #Read in the averaged and std of spectra
meta <- read.csv(metaFile) #Read in associated metadata of spectra
wavelength <- as.numeric(gsub('X','',names(spectra)[2:length(names(spectra))])) #Get wavelengths

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
    geom_histogram() + #plot as histogram
    ggtitle(wavelength[x-1]) + #add title to histogram
    xlab("Reflectance (%)") #add xlabel
  
  plotName <- paste(directory,"Normality_Test\\HyTES_",round(wavelength[x-1]*1000,digits = 0),"_hist.png",sep="") #create plot name to save the file
  ggsave(plotName,width = 3,height = 3) #save the plot
}
fileName <- paste(directory,"Normality_Test\\hytes_normality_test_results.csv",sep="") 
write(t(normTest),file = fileName,sep = ",")#save normality test data to file

### Kruskal-Wallis ###
##Since data is non-parametric, Kruskal-Wallis is need to compare samples
kwResults <- array(0,dim = c(length(wavelength),2))#Empty array that will hold values below
kwResults[,1] <- wavelength #add wavelength to the first column
#dunnResults <- array(3,0)#Empty array that will hold values below

fileNameKW <- paste(directory,"Kruskal_Wallis_Test\\hytes_kruskal_wallis_test_results.csv",sep="") #create filename that will hold kw test results
write("Wavelength,KW_pvalue",file = fileNameKW)#write header to file
fileNameKWPH <- paste(directory,"Kruskal_Wallis_Test\\hytes_kruskal_wallis_posthoc_test_results.csv",sep="") #create filename that will hold kw post hoc test results
write("Wavelength,Pair,Obs-Diff,Critical-Diff,Flag",file = fileNameKWPH)
fileNameKWDT <- paste(directory,"Kruskal_Wallis_Test\\hytes_dunn_posthoc_test_results.csv",sep="") #create filename that will hold kw post hoc test results
write("Wavelength,Pair,P-value",file = fileNameKWDT)

for(x in 2:length(spectra)){ #Loop through wavelengths 
  #All Wavelengths = c(2:length(spectra)) 
  #Only Some Wavelengths = as.numeric(c(2,3,4))
  print(x)
  
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
  kwDunnTest <- dunn.test(spectra[,x],meta$Acronym)
  
  for(i in c(1:length(kwDunnTest$P))){
    if(kwDunnTest$P[i] < 0.05){
      #dunnResults <- rbind(dunnResults,c(wavelength[x-1],kwDunnTest$comparisons[i],kwDunnTest$P[i]))
      output <- paste(wavelength[x-1],kwDunnTest$comparisons[i],kwDunnTest$P[i],sep = ",")
      write(output,file = fileNameKWDT,append=TRUE) #add line to file 
    }
  }
}

### Kruskal-Wallis Post Hoc Analysis ###
dunnFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\Kruskal_Wallis_Test\\hytes_dunn_posthoc_test_results.csv" #Set to metadata file location
dunnResults <- read.csv(dunnFile) #Read in the averaged and std of spectra

#Plot of frequency of pairs for wavelength
windows() #open it in a new window
ggplot(data = dunnResults, aes(dunnResults$Wavelength)) + #Plot the data
  geom_histogram(binwidth = .4, fill = "gray",col = "black",boundary = 0) + #Plot it as a Histogram
  theme_set(theme_bw(base_size=20)) + #Set the theme 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3000)) + #removing offset on y axis
  scale_x_continuous(breaks = seq(8,12,.5)) + #changing tick mark frequency in x axis
  labs(x = expression(paste("Wavelength (",mu,"m)")), y = "Frequency") #Add labels to plot
plotName = paste(directory,"Kruskal_Wallis_Test\\hytes_Pair_Wavelength_Frequency.png",sep="") #create plot name to save the file
ggsave(plotName,width = 10,height = 7) #save the plot

