%% Nicolet-Spectra-Analysis 
% Susan Meerdink
% 10/12/2016
% This code reads in Nicolet spectra that have been averaged using 
% https://github.com/susanmeerdink/ASD-Nicolet-Spectra-Processing.
% The result of this Java program is a .csv file with a two rows for each
% spectra process: averaged spectra and standard deviation.
% This code will analyze the spectra to determine similarities and differences.
%% Import Data
directory = 'C:\\Users\\Susan\\Documents\\GitHub\\ASD-Nicolet-Spectra-Analysis\\'; %Home directory
dataFile = strcat(directory,'Nicolet_Averaged_Spectra_Huntington_rep3.csv'); %Set to spectra file location
metaFile = strcat(directory,'Nicolet_Averaged_Spectra_Metadata_Huntington_rep3.csv'); %Set to metadata file location
data = readtable(dataFile); %Read in the averaged and std of spectra
metaTable = readtable(metaFile); %Read in associated metadata of spectra
allMeta = table2cell(metaTable); %Convert to cell array
wavelengths = csvread(dataFile,0,1,[0 1 0 1738]); %pull out nicolet wavelengths
allSpectra = table2cell(data); %convert to cell array

dunnFile = strcat(directory,'\\Kruskal_Wallis_Test\\dunn_posthoc_test_results.csv');
dunnWave = readtable(dunnFile); %pull out nicolet wavelengths
allDunn = table2cell(dunnWave(:,1)); %convert to cell array

%% F-Ratio
fratio = zeros(1,length(wavelengths));
pAll = [];
wAll =[];
bAll = [];
for b = 2: length(wavelengths)
    [d,p,stats] = manova1(cell2mat(allSpectra(:,b)),allMeta(:,2));
    fratio(b) = stats.W/stats.B;
    wAll = vertcat(wAll,stats.W);
    bAll = vertcat(bAll,stats.B);
    pAll = vertcat(p, pAll);
end
normfratio = 1./fratio; %now low fratio = less important wavelengths
%% Plot F-Ratio Results
figure('units','normalized','outerposition',[0 0 1 1])
hold on
plot(wavelengths,fratio)
set(gca,'Xlim',[2.5 15],'Ylim',[0 2],'YTick',[0:.1:1],'XTick',[2.5:2:15]) %
set(gca,'FontSize',30)
set(gca,'Fontname','Calibri')
xlabel(['Wavelength ( \mum )']) % label x-axis
ylabel('F-ratio') % label left y-axis
hold off

%% Histogram
figure('units','normalized','outerposition',[0 0 0.85 1])
hold on

%Create a second axes in the same location as the first axes by setting the position of the second axes equal to the position of the first axes. 
ax1 = gca; %current axis
AxesHandle = findobj(gcf,'Type','axes');
ax1_pos = get(AxesHandle,'Position'); % position of first axes
ax2 = axes('Position',ax1_pos,'YAxisLocation','right','Color','none');
histogram(cell2mat(allDunn),'FaceColor','b','Parent',ax1) %,'Parent',ax2

input = (100 - mean(cell2mat(allSpectra(:,[2:length(allSpectra)]))))/100;
line(wavelengths,input,'Color','r','LineWidth',1.5,'Parent',ax2)

set(ax1,'Xlim',[2.5 15],'XTick',[2.5:2:15]) %,'Ylim',[0 ],'YTick',[0:0.01: 0.06]
set(ax2,'Ylim',[.88 1],'YTick',[.88:.02:1],'Xlim',[2.5 15],'XTick',[]); %
xlabel(ax1,['Wavelength ( \mum )']) % label x-axis
ylabel(ax1,'Frequency') % label left y-axis
ylabel(ax2,'Emissivity') % label left y-axis
set(ax1, 'FontSize',30)
set(ax2, 'FontSize',30)
hold off