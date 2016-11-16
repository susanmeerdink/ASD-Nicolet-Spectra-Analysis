%% Nicolet-HyTES-Comparison-Analysis 
% Susan Meerdink
% 11/11/2016
% This code reads in Nicolet and HyTES spectra.
% This code will analyze the spectra at the leaf and canopy level 
% to determine similarities and differences.
% run on Matlab 2016B
%% Load Data

%Read in Nicolet (leaf level) data
dirNicolet = 'F:\\Dropbox\\Analysis\\JPL Analysis\\Nicolet spectra\\'; %Home directory
dataFileNicolet = strcat(dirNicolet,'Nicolet_Averaged_Huntington_Spectra.csv'); %Set to spectra file location
metaFileNicolet = strcat(dirNicolet,'Nicolet_Averaged_Huntington_Metadata.csv'); %Set to metadata file location
metaNicolet = table2cell(readtable(metaFileNicolet)); %Convert to cell array
waveNicolet = csvread(dataFileNicolet,0,2,[0 2 0 1739]); %pull out nicolet wavelengths
dataNicolet = readtable(dataFileNicolet); %Read in the averaged and std of spectra
spectraNicolet = dataNicolet(strcmp(table2cell(dataNicolet(:,2)),'AVG'),:); %pulling out AVG values NOT STD
spectraNicolet = cell2mat(table2cell(spectraNicolet(:,[3:1740]))); %convert to cell array
spectraNicolet = (100-spectraNicolet)/100; %converting to emissivity

%Read in HyTES (canopy level) data
dirHyTES = 'F:\\Dropbox\\Analysis\\JPL Analysis\\HyTES spectra\\'; %Home directory
dataFileHyTES = strcat(dirHyTES,'HyTES_Spectra_20161101.csv'); %Set to spectra file location
metaFileHyTES = strcat(dirHyTES,'HyTES_Metadata_20161101.csv'); %Set to metadata file location
dataHyTES = readtable(dataFileHyTES); %Read in the averaged and std of spectra
metaHyTES = table2cell(readtable(metaFileHyTES)); %Convert to cell array
waveHyTES = csvread(dataFileHyTES,0,1,[0 1 0 202]); %pull out wavelengths
spectraHyTES = cell2mat(table2cell(dataHyTES(:,[2:203]))); %convert to cell array

species = {'ALAR';'BABE';'BATU';'BRDI';'BRRU';'CACA';'CALE';'CEDE';'CHIN';'CHSP';'FICO';'FITH';'JAMI';'LAIN';'MAGR';'MELI';'PEAF';'PHVI';'POGR';'QUAG';'QUIL';'QURO';'QUSU';'QUVI';'SABA';'TAMU';'TITI'};
dirOut = 'F:\\Dropbox\\Analysis\\JPL Analysis\\Nicolet_HyTES_Comparison\\'; %Output directory

%% Brewer Colors
red = [228 26 28] ./ 255;
blue = [55 126 184] ./ 255;
green = [77 175 74] ./ 255;
purple = [152 78 163] ./ 255;
orange = [255 127 0] ./ 255;

%% Averaging Species
avgSpecNicolet = zeros(27,1738);
avgSpecHyTES = zeros(27,202);
stdSpecNicolet = zeros(27,1738);
stdSpecHyTES = zeros(27,202);

for a = 1:size(species,1) %loop through species
    specTemp = spectraNicolet(strcmp(metaNicolet(:,2),species(a)),:);
    avgSpecNicolet(a,:) = mean(specTemp);
    stdSpecNicolet(a,:) = std(specTemp);
    
    specTemp = spectraHyTES(strcmp(metaHyTES(:,3),species(a)),:);
    avgSpecHyTES(a,:) = mean(specTemp);
    stdSpecHyTES(a,:) = std(specTemp);
end

%% Plots for Comparison

for a = 1:size(species,1) %loop through species
    %Zoomed Out Image
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    
    %Plotting Nicolet
    p = area(waveNicolet, avgSpecNicolet(a,:)+stdSpecNicolet(a,:),'FaceColor',red,'EdgeColor','none','FaceAlpha',0.3);
    area(waveNicolet, avgSpecNicolet(a,:)-stdSpecNicolet(a,:),'FaceColor','w','EdgeColor','none');
    plot(waveNicolet,avgSpecNicolet(a,:),'Color',red,'LineWidth',2)
    
    %Plotting HyTES
    area(waveHyTES, avgSpecHyTES(a,:)+stdSpecHyTES(a,:),'FaceColor',blue,'EdgeColor','none','FaceAlpha',0.3);
    area(waveHyTES, avgSpecHyTES(a,:)-stdSpecHyTES(a,:),'FaceColor','w','EdgeColor','none');
    plot(waveHyTES,avgSpecHyTES(a,:),'Color',blue,'LineWidth',2)
    
    %Figure Features
    set(gca,'FontSize',24,'FontName','Cambria')
    xlabel(('Wavelength (\mum)')) % label x-axis
    ylabel('Emissivity')
    set(gca,'Xlim',[2.5 15],'XTick',[2.5:2.5:15])
    set(gca,'Ylim',[.88 1],'YTick',[.88:.03:1.0])
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
    set(gca,'ygrid','on')
    title(char(species(a)))
    hold off
    
    %Save Image
    set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 7.25 3.5])
    nameFile = strcat(dirOut,'\\Figures\\',char(species(a)),'_Compare_Zoom_Out');
    print(nameFile,'-dpng','-r0')
    
    % Zoomed In
%     figure('units','normalized','outerposition',[0 0 1 1])
%     hold on
%     plot(waveNicolet,avgSpecNicolet(a,:),'Color',red,'LineWidth',2)
%     plot(waveHyTES,avgSpecHyTES(a,:),'Color',blue,'LineWidth',2)
%     set(gca,'FontSize',24,'FontName','Cambria')
%     xlabel(('Wavelength (\mum)')) % label x-axis
%     ylabel('Emissivity')
%     set(gca,'Xlim',[8 11.5],'XTick',[8:.5:11.5])
%     set(gca,'Ylim',[.88 1],'YTick',[.88:.03:1.0])
%     set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
%     set(gca,'ygrid','on')
%     hold off
%     
%     set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 7 3.5])
%     nameFile = strcat(dirOut,'\\Figures\\',char(species(a)),'_Compare_Zoom_In');
%     print(nameFile,'-dpng','-r0')
end

%% END
close all