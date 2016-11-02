%% HyTES-Spectra-Analysis 
% Susan Meerdink
% 10/16/2016
% This code reads in Nicolet spectra that have been averaged using 
% https://github.com/susanmeerdink/ASD-Nicolet-Spectra-Processing.
% The result of this Java program is a .csv file with a two rows for each
% spectra process: averaged spectra and standard deviation.
% This code will analyze the spectra to determine similarities and differences.
%% Import Data
directory = 'F:\\Dropbox\\Analysis\\JPL Analysis\\HyTES ROI spectra\\'; %Home directory
dataFile = strcat(directory,'HyTES_Spectra_20161101.csv'); %Set to spectra file location
metaFile = strcat(directory,'HyTES_Metadata_20161101.csv'); %Set to metadata file location
data = readtable(dataFile); %Read in the averaged and std of spectra
metaTable = readtable(metaFile); %Read in associated metadata of spectra
allMeta = table2cell(metaTable); %Convert to cell array
wavelengths = csvread(dataFile,0,1,[0 1 0 202]); %pull out wavelengths
allSpectra = cell2mat(table2cell(data)); %convert to cell array

%% Find average of species
acronym = unique(allMeta(:,3));
avgSpectra = zeros(27,202);
stdSpectra = zeros(27,202);

for a = 1:length(acronym)
    spectra = cell2mat(allSpectra(strcmp(allMeta(:,3),acronym(a)),:));
    avgSpectra(a,:) = mean(spectra(:,[2:203]));
    stdSpectra(a,:) = std(spectra(:,[2:203]));   
end

%% Display species spectra

for a = 1:length(acronym)
    figure('units','normalized','outerposition',[0 0 1 1])
    y = cell2mat(allSpectra(strcmp(allMeta(:,3),acronym(a)),[2:203]))';
    x = repmat(wavelengths,size(y,2),1)';
    line(x,y);
    text(8.02,0.99,num2str(size(y,2)));
    axis([8 11.6 0.85 1])
    title(acronym(a));
end
%% Kruskal Wallis Test and Dunn Post Hoc Test
pValue = []; %Empty array to hold pvalues from kruskal wallis
pairs = []; %empty array to hold results from dunn test
for w = 1:size(wavelengths,2) %loopthrough wavelengths
    [p,tbl,stats] = kruskalwallis(allSpectra(:,(w+1)),allMeta(:,3),'off');
    pValue = vertcat(pValue,[wavelengths(w),p]);
    
    if p < 0.05 %if that wavelength is statistically significant, do post hoc test
        c = multcompare(stats,'CType','dunn-sidak','Display','off');
        for i = 1:size(c,1) %loop through resulting pairs
            if c(i,6) < 0.05 %if the pair is significantly different add
                input = [wavelengths(w),stats.gnames(c(i,1)),stats.gnames(c(i,2)),c(i,6)];
                pairs = vertcat(pairs, input);    
            end
        end
    end
end
plot(pValue(:,1),pValue(:,2))
%% Histogram
figure('units','normalized','outerposition',[0 0 0.85 1])
hold on

%Create a second axes in the same location as the first axes by setting the position of the second axes equal to the position of the first axes. 
ax1 = gca; %current axis
AxesHandle = findobj(gcf,'Type','axes');
ax1_pos = get(AxesHandle,'Position'); % position of first axes
ax2 = axes('Position',ax1_pos,'YAxisLocation','right','Color','none');
histogram(pairs,'FaceColor',[0/255 0/255 153/255],'Parent',ax1) %

input = allSpectra(13,[2:size(allSpectra,2)]);
line(wavelengths,input,'Color','r','LineWidth',1.5,'Parent',ax2)

set(ax1,'Xlim',[8 11.6],'XTick',[8:.5:12]) %,'Ylim',[0 ],'YTick',[0:0.01: 0.06]
set(ax2,'Ylim',[.8 1],'Xlim',[8 11.5],'XTick',[]); %'YTick',[.98:.2:1]
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
histogram(cell2mat(allDunn),'FaceColor',[0/255 0/255 153/255]) %

set(gca,'Xlim',[8 11.6],'XTick',[8:.5:12]) %,'Ylim',[0 ],'YTick',[0:0.01: 0.06]
xlabel(['Wavelength ( \mum )']) % label x-axis
ylabel('Frequency') % label left y-axis
set(gca, 'FontSize',30)
hold off