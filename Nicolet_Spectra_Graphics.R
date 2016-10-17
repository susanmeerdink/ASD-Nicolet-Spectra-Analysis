### Nicolet-Spectra-Graphics ###
## Susan Meerdink
## 10/14/2016
## This code reads in Nicolet spectra processed by Nicolet_Data_Prep.R 
## This code will display spectra.
###---------------------------------------------------------------------------------------------------------------------------- ###
### Setting Up Environment ###
library(ggplot2) #load ggplot2 for graphics/figures
library(reshape2) #load for melting function

### Reading in Data ###
directory <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\"
dataFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\Nicolet_Averaged_Spectra_Huntington_rep3.csv" #Set to spectra file location
metaFile <- "C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\Nicolet_Averaged_Spectra_Metadata_Huntington_rep3.csv" #Set to metadata file location
spectra <- read.csv(dataFile) #Read in the averaged and std of spectra
meta <- read.csv(metaFile) #Read in associated metadata of spectra
wavelength <- as.numeric(gsub('X','',names(spectra)[2:length(names(spectra))])) #Get wavelengths

### Plotting Spectra  by species ###
## WARNING: will produce many image files!
acronym <- unique(meta$Acronym) #Get unique listing of acronyms
for(a in acronym){
  plotData <- data.frame(wl = wavelength, t(spectra[which(meta$Acronym %in% a),2:ncol(spectra)])) #make current data into data.frame for plotting
  colnames(plotData) <- c("wl",as.character(spectra[which(meta$Acronym %in% a),1])) #rename columns with ID
  plotDataMelt <-  melt(plotData,id.vars = 'wl') #Melt using wavelength as the ID variable (puts it in format to plot)
  
  ggplot(data = plotDataMelt, aes(x = wl,y = value, group = variable,color=variable)) + #Plotting the spectra
    geom_line()+ #plot as a line
    xlab(expression(paste("Wavelength (",mu,"m)",sep ="")))+ #add xlabel
    coord_cartesian(xlim = c(2.5,13), ylim = c(0,15)) + #Fix axis
    ylab("Reflectance (%)")+ #add ylabel
    ggtitle(a) #add title
  
  plotName = paste(directory,"Spectra_Plots_Nicolet\\",a,"_spectra.png",sep="") #create plot name to save the file
  ggsave(plotName,width = 10,height = 5) #save the plot
}

## Grouping Data by Species ##
avgSpectra <- vector() #empty vector that will hold averaged species spectra (all samples into one spectrum)
for(a in acronym){ #Loop through acronyms
  temp <- data.frame(spectra[which(meta$Acronym %in% a),2:ncol(spectra)])
  avgSpectra <- rbind(avgSpectra,colMeans(temp)) #add to avgSpectra variable
}

## Plotting Overall Species ##
plotData <- data.frame(wl = wavelength, t(avgSpectra)) #make current data into data.frame for plotting
colnames(plotData) <- c("wl",as.character(acronym )) #rename columns with ID
plotDataMelt <-  melt(plotData,id.vars = 'wl') #Melt using wavelength as the ID variable (puts it in format to plot)

windows()
ggplot(data = plotDataMelt, aes(x = wl,y = value, group = variable,color=variable)) + #Plotting the spectra
  geom_line()+ #plot as a line
  theme_set(theme_bw(base_size=20)) + #Set the theme 
  coord_cartesian(xlim = c(2.5,15), ylim = c(0,15)) + #Fix axis
  scale_y_continuous(breaks = seq(0,15,2.5),expand = c(0, 0), limits = c(0, 15)) + #removing offset on y axis
  scale_x_continuous(breaks = seq(2.5,15,1.5),expand = c(0,0)) + #changing tick mark frequency in x axis
  labs(x = expression(paste("Wavelength (",mu,"m)")), y = "Reflectance (%)") #Add labels to plot
plotName = paste(directory,"Spectra_Plots_Nicolet\\all_species_spectra.png",sep="") #create plot name to save the file
ggsave(plotName,width = 10,height = 5) #save the plot

