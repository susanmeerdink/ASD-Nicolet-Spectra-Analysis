%% HyTES-Spectra-Analysis 
% Susan Meerdink
% 10/16/2016
% This code reads in Nicolet spectra that have been averaged using 
% https://github.com/susanmeerdink/ASD-Nicolet-Spectra-Processing.
% The result of this Java program is a .csv file with a two rows for each
% spectra process: averaged spectra and standard deviation.
% This code will analyze the spectra to determine similarities and differences.
%% Import Data
directory = 'F:\\Dropbox\\Analysis\\JPL Analysis\\HyTES spectra\\'; %Home directory
dataFile = strcat(directory,'HyTES_Spectra_20161101.csv'); %Set to spectra file location
metaFile = strcat(directory,'HyTES_Metadata_20161101.csv'); %Set to metadata file location
data = readtable(dataFile); %Read in the averaged and std of spectra
metaTable = readtable(metaFile); %Read in associated metadata of spectra
allMeta = table2cell(metaTable); %Convert to cell array
wavelengths = csvread(dataFile,0,1,[0 1 0 202]); %pull out wavelengths
allSpectra = cell2mat(table2cell(data)); %convert to cell array

%% Find average/min/max/std of species
acronym = unique(allMeta(:,3));
avgSpectra = zeros(27,202);
stdSpectra = zeros(27,202);
minSpectra = zeros(27,202);
maxSpectra = zeros(27,202);

for a = 1:length(acronym)
    spectra = allSpectra(strcmp(allMeta(:,3),acronym(a)),:);
    avgSpectra(a,:) = mean(spectra(:,[2:203]));
    stdSpectra(a,:) = std(spectra(:,[2:203]));
    minSpectra(a,:) = min(spectra(:,[2:203]));
    maxSpectra(a,:) = max(spectra(:,[2:203]));
    
    % Display average and std for species
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    area(wavelengths, avgSpectra(a,:)+stdSpectra(a,:),'FaceColor',[159/255 182/255 205/255],'EdgeColor','none');
    area(wavelengths, avgSpectra(a,:)-stdSpectra(a,:),'FaceColor','w','EdgeColor','none');
    line(wavelengths, avgSpectra(a,:),'Color','k','LineWidth',2)
    line(wavelengths, minSpectra(a,:),'Color','b','LineWidth',2)
    line(wavelengths, maxSpectra(a,:),'Color','b','LineWidth',2)
    title(acronym(a)) 
    xlabel('Wavelength');
    ylabel('Emissivity');
    set(gca,'FontSize',24);
    axis([8 11.6 .85 1])
    hold off
end

figure('units','normalized','outerposition',[0 0 .5 1])
plot(wavelengths,avgSpectra((1:9),:))
legend(acronym(1:9))
axis([8 11.5 0.85 1])

figure('units','normalized','outerposition',[0 0 .5 1])
plot(wavelengths,avgSpectra((10:18),:))
legend(acronym(10:18))
axis([8 11.5 0.85 1])

figure('units','normalized','outerposition',[0 0 .5 1])
plot(wavelengths,avgSpectra((19:27),:))
legend(acronym(19:27))
axis([8 11.5 0.85 1])
%% Display species spectra
close all
for a = 1:length(acronym)
    figure('units','normalized','outerposition',[0 0 1 1])
    y = allSpectra(strcmp(allMeta(:,3),acronym(a)),[2:203])';
    x = repmat(wavelengths,size(y,2),1)';
    line(x,y);
    text(8.02,0.99,num2str(size(y,2)));
    axis([8 11.6 0.85 1])
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
%close all
figure('units','normalized','outerposition',[0 0 1 1])
hold on
plot(wavelengths,avgSpectra(strcmp(acronym,'BATU'),:),'Color',red,'LineWidth',2)
plot(wavelengths,avgSpectra(strcmp(acronym,'JAMI'),:),'Color',blue,'LineWidth',2)
plot(wavelengths,avgSpectra(strcmp(acronym,'MAGR'),:),'Color',green,'LineWidth',2)
plot(wavelengths,avgSpectra(strcmp(acronym,'PEAF'),:),'Color',purple,'LineWidth',2)
plot(wavelengths,avgSpectra(strcmp(acronym,'QURO'),:),'Color',orange,'LineWidth',2)
set(gca,'FontSize',40,'FontName','Cambria')
xlabel(['Wavelength (\mum)']) % label x-axis
ylabel('Emissivity')
set(gca,'Xlim',[8 11.5],'XTick',[8:.5:11.5])
set(gca,'Ylim',[.85 1],'YTick',[.85:.03:1.0])
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
set(gca,'ygrid','on')
hold off
%% Kruskal Wallis Test and Dunn Post Hoc Test
close all
pValue = []; %Empty array to hold pvalues from kruskal wallis
pairs = []; %empty array to hold results from dunn test
norm = [];
for w = 1:size(wavelengths,2) %loopthrough wavelengths
    [p,tbl,stats] = kruskalwallis(allSpectra(:,(w+1)),allMeta(:,3),'off'); %
    pValue = vertcat(pValue,[wavelengths(w),p]);
    
    %Multiple comparision test
    %based on seigel and castellan 1988 pg 213
    diff = abs(bsxfun(@minus,stats.meanranks,stats.meanranks'));% Absolute pairwise diifferences
    
    %plotting to look at histogram to determine if it normal
    r = reshape(triu(diff),[1,729]);
    r(find(r ==0)) = [];
    figure
    hist(r)
  
    zScore = @(p) sqrt(2) * erfcinv(p*2); %calculate z score from p value
    nSamples = size(allSpectra,1); %number of samples
    nGroups = size(stats.gnames,1); %number of groups
    n1 = (nSamples*(nSamples+1))/12;
    n2 = bsxfun(@plus,(1./stats.n),(1./stats.n)'); 
    n3 = sqrt(n1*n2);
    p = 0.05/(nGroups*(nGroups-1)); %calculate the p value
    z = zScore(p); %get the zscore of the pvalue
    compare = z.*n3; %Critical value to compare the mean rank differences
    result = diff >= compare; %determine where the mean rank is larger than or equal to the critical value
    [row,col] = find(triu(result) == 1); %get the row (group1) and column (group2) of the significantly different pairs
    
    %Save the significantly different groups with the wavelength
    input = [repmat(wavelengths(w),[1,size(row,1)])',row,col]; 
	pairs = vertcat(pairs, input); 
end
figure
plot(pValue(:,1),pValue(:,2))
%% Histogram
figure('units','normalized','outerposition',[0 0 0.85 1])
hold on

%Create a second axes in the same location as the first axes by setting the position of the second axes equal to the position of the first axes. 
ax1 = gca; %current axis
AxesHandle = findobj(gcf,'Type','axes');
ax1_pos = get(AxesHandle,'Position'); % position of first axes
ax2 = axes('Position',ax1_pos,'YAxisLocation','right','Color','none');
histogram(pairs(:,1),'FaceColor',blue,'Parent',ax1) %cell2mat(

input = allSpectra(13,[2:size(allSpectra,2)]);
line(wavelengths,input,'Color','r','LineWidth',1.5,'Parent',ax2)

set(ax1,'Xlim',[8 11.6],'XTick',[8:.5:12]) %,'Ylim',[0 ],'YTick',[0:0.01: 0.06]
set(ax2,'Ylim',[.85 1],'Xlim',[8 11.5],'XTick',[]); %'YTick',[.98:.2:1]
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
histogram(pairs(:,1),'FaceColor',blue) %cell2mat(

set(gca,'Xlim',[8 11.6],'XTick',[8:.5:12]) %,'Ylim',[0 ],'YTick',[0:0.01: 0.06]
xlabel(['Wavelength ( \mum )']) % label x-axis
ylabel('Frequency') % label left y-axis
set(gca, 'FontSize',30)
hold off

%% Separable Species
allPairs = vertcat(pairs(:,2),pairs(:,3));
[p1 p2 p3] = unique(allPairs);
d = hist(p3, length(p1));
[sorted,indexSorted] = sort(d);

figure('units','normalized','outerposition',[0 0 1 1])
bar(sorted)
set(gca,'XTick',1:1:27,'XTickLabel',acronym(indexSorted));
%% END
close all