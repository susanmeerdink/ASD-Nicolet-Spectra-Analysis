# Combining ROIs into one file
# Susan Meerdink
# Created 11/30/2016
# This code reads through a directory that contains formatted ROIs
# This code takes two inputs:
# 1. The directory containing the formatted spectral libraries
# 2. Output file name

# ############INPUTS######################################
directory = 'F:\\Dropbox\\Analysis\\JPL Analysis\\HyTES spectra\\FORMATTED Output from ROIs - Temp GeoRef Image\\'

outputName = 'F:\\Dropbox\\Analysis\\JPL Analysis\\HyTES spectra\\HyTES_Temp_GeoRef.csv'

############ENDINPUTS####################################

import glob, os
os.chdir(directory)
outputFile = open(outputName,'w') #Create the output file for the formatted ROIs
fileCount = 0 #keeps track of number of files processed

for file in glob.glob("*.csv"): #loop through files in directory
    print(file)
    inputFile = open(file,'r') #Open the ROI file
    for line in inputFile: #loop through file
        outputFile.write(line) #save
    inputFile.close()
                
outputFile.close()

print('Files Combined')           

        
    
