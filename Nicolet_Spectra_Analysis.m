%% Nicolet-Spectra-Analysis 
% Susan Meerdink
% 10/12/2016
% This code reads in Nicolet spectra that have been averaged using 
% https://github.com/susanmeerdink/ASD-Nicolet-Spectra-Processing.
% The result of this Java program is a .csv file with a two rows for each
% spectra process: averaged spectra and standard deviation.
% This code will analyze the spectra to determine similarities and differences.
%% Import Data
directory = 'F:\\Dropbox\\Analysis\\JPL Analysis\\Nicolet spectra\\'; %Home directory
dataFile = strcat(directory,'Nicolet_Averaged_Huntington_Spectra.csv'); %Set to spectra file location
metaFile = strcat(directory,'Nicolet_Averaged_Huntington_Metadata.csv'); %Set to metadata file location
data = readtable(dataFile); %Read in the averaged and std of spectra
metaTable = readtable(metaFile); %Read in associated metadata of spectra
allMeta = table2cell(metaTable); %Convert to cell array
wavelengths = csvread(dataFile,0,2,[0 2 0 1739]); %pull out nicolet wavelengths

%% Processing Input Data
allSpectra = data(strcmp(table2cell(data(:,2)),'AVG'),:); %pulling out AVG values NOT STD
allSpectra = cell2mat(table2cell(allSpectra(:,[3:1740]))); %convert to cell array
allSpectra = (100-allSpectra)/100; %converting to emissivity

%Finding Species for analysis
%These are harded coded in from HyTES_Spectra_Analysis unique function
species = {'ALAR';'BABE';'BATU';'BRDI';'BRRU';'CACA';'CALE';'CEDE';'CHIN';'CHSP';'FICO';'FITH';'JAMI';'LAIN';'MAGR';'MELI';'PEAF';'PHVI';'POGR';'QUAG';'QUIL';'QURO';'QUSU';'QUVI';'SABA';'TAMU';'TITI'};
repSpectra = []; %holds species spectra that have replication of 3 or more
repMetadata = []; %holds species metadata that have replication of 3 or more
for s = 1:size(species,1) %loop through the species
    index = find(strcmp(allMeta(:,2),species(s)));
    repSpectra = vertcat(repSpectra, allSpectra(index,:));
    repMetadata = vertcat(repMetadata,allMeta(index,:));
end
%% Find average/min/max/std of species
avgSpectra = zeros(27,1738);
stdSpectra = zeros(27,1738);
minSpectra = zeros(27,1738);
maxSpectra = zeros(27,1738);

for a = 1:size(species,1)
    spectra = repSpectra(strcmp(repMetadata(:,2),species(a)),:);
    avgSpectra(a,:) = mean(spectra);
    stdSpectra(a,:) = std(spectra);
    minSpectra(a,:) = min(spectra);
    maxSpectra(a,:) = max(spectra);
    
    % Display average and std for species
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    area(wavelengths, avgSpectra(a,:)+stdSpectra(a,:),'FaceColor',[159/255 182/255 205/255],'EdgeColor','none');
    area(wavelengths, avgSpectra(a,:)-stdSpectra(a,:),'FaceColor','w','EdgeColor','none');
    line(wavelengths, avgSpectra(a,:),'Color','k','LineWidth',2)
    line(wavelengths, minSpectra(a,:),'Color','b','LineWidth',2)
    line(wavelengths, maxSpectra(a,:),'Color','b','LineWidth',2)
    title(species(a)) 
    xlabel('Wavelength');
    ylabel('Emissivity');
    set(gca,'FontSize',24);
    axis([2.5 15 .85 1])
    hold off
end
%% Display species spectra
close all
acronym = unique(allMeta(:,2));
for a = 1:length(acronym)
    figure('units','normalized','outerposition',[0 0 1 1])
    y = allSpectra(strcmp(allMeta(:,2),acronym(a)),:)';
    x = repmat(wavelengths,size(y,2),1)';
    line(x,y);
    text(8.02,0.99,num2str(size(y,2)));
    axis([2.5 15 0.85 1])
    xlabel('Wavelength');
    ylabel('Emissivity');
    set(gca,'FontSize',24);
    title(acronym(a));
end
%% Brewer Colors
red = [228 26 28] ./ 255;
blue = [55 126 184] ./ 255;
green = [77 175 74] ./ 255;
purple = [152 78 163] ./ 255;
orange = [255 127 0] ./ 255;

%% figures for only 5 species - POSTER
close all
figure('units','normalized','outerposition',[0 0 1 1])
hold on
plot(wavelengths,avgSpectra(strcmp(species,'BATU'),:),'Color',red,'LineWidth',2)
plot(wavelengths,avgSpectra(strcmp(species,'JAMI'),:),'Color',blue,'LineWidth',2)
plot(wavelengths,avgSpectra(strcmp(species,'MAGR'),:),'Color',green,'LineWidth',2)
plot(wavelengths,avgSpectra(strcmp(species,'PEAF'),:),'Color',purple,'LineWidth',2)
plot(wavelengths,avgSpectra(strcmp(species,'QURO'),:),'Color',orange,'LineWidth',2)
set(gca,'FontSize',40,'FontName','Cambria')
xlabel(['Wavelength (\mum)']) % label x-axis
ylabel('Emissivity')
set(gca,'Xlim',[2.5 15],'XTick',[2.5:2.5:15])
set(gca,'Ylim',[.85 1],'YTick',[.85:.03:1.0])
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
set(gca,'ygrid','on')
hold off

%% Non Parametric Tests
close all

% Get the pair order
n = size(species,1);
[a,b] = meshgrid(1:n, 1:n);
mask   = triu(ones(n), 1) > 0.5;
pairs  = [a(mask) b(mask)];

pValue = []; %Empty array to hold pvalues from mann - whitney u test
results = [];
results(:,1:2) = pairs; %add onto this array with pvalue results

for w = 1:size(wavelengths,2) %loopthrough wavelengths
    for i = 1: size(pairs,1) %loop through the pair list
        g1Index = strmatch(char(species(pairs(i,1))),repMetadata(:,2));
        g2Index = strmatch(char(species(pairs(i,2))),repMetadata(:,2));
        g1 = repSpectra(g1Index,w);
        g2 = repSpectra(g2Index,w);
        p = ranksum(g1,g2);
        results(i,2+w) = p;
    end 
end

%% Coding up Results
waveSummary = results(:,3:1740) < 0.05; %mark the values with a 1 if they are less then 0.05
waveTotal = sum(waveSummary,1); %get the total number of pairs that were significantly different at that wavelength
waveList = [];
for l = 1:size(wavelengths,2)
    add = repmat(wavelengths(l),[waveTotal(l),1]);
    waveList = vertcat(waveList,add);
end

pairTotal = sum(waveSummary,2); %get the total number of times a pair shows up
pairList = [];
for l = 1: size(pairTotal,1)
    n1 = species(pairs(l,1));
    n2 = species(pairs(l,2));
    add1 = repmat(n1,[pairTotal(l),1]);
    add2 = repmat(n2,[pairTotal(l),1]);
    pairList = vertcat(pairList,add1,add2);
end
%% Histogram
figure('units','normalized','outerposition',[0 0 0.85 1])
hold on

%Create a second axes in the same location as the first axes by setting the position of the second axes equal to the position of the first axes. 
ax1 = gca; %current axis
AxesHandle = findobj(gcf,'Type','axes');
ax1_pos = get(AxesHandle,'Position'); % position of first axes
ax2 = axes('Position',ax1_pos,'YAxisLocation','right','Color','none');
histogram(waveList,'FaceColor',blue,'Parent',ax1) %

input = mean(avgSpectra);
line(wavelengths,input,'Color','r','LineWidth',1.5,'Parent',ax2)

set(ax1,'Xlim',[2.5 15],'XTick',[2.5:2.5:15]) %,'Ylim',[0 ],'YTick',[0:0.01: 0.06]
set(ax2,'Ylim',[.85 1],'Xlim',[2.5 15],'XTick',[]); %'YTick',[.98:.2:1]
xlabel(ax1,['Wavelength ( \mum )']) % label x-axis
ylabel(ax1,'Frequency') % label left y-axis
ylabel(ax2,'Emissivity') % label left y-axis
set(ax1, 'FontSize',30)
set(ax2, 'FontSize',30)
hold off

%% Histogram without Emissivity
figure('units','normalized','outerposition',[0 0 0.85 1])
hold on

%Create a second axes in the same location as the first axes by setting the position of the second axes equal to the position of the first axes. 
histogram(waveList,'FaceColor',blue) %

set(gca,'Xlim',[2.5 15],'XTick',[2.5:2.5:15]) %,'Ylim',[0 ],'YTick',[0:0.01: 0.06]
xlabel(['Wavelength ( \mum )']) % label x-axis
ylabel('Frequency') % label left y-axis
set(gca, 'FontSize',30)
hold off

%% Separable Species
[p1, p2, p3] = unique(pairList);
d = hist(p3, length(p1));
[sorted,indexSorted] = sort(d);

figure('units','normalized','outerposition',[0 0 1 1])
bar(sorted)
set(gca,'XTick',1:1:27,'XTickLabel',species(indexSorted));
%% END
close all