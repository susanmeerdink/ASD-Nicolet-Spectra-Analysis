%% Load ENVI Image for Visualization
%Susan Meerdink
%11/30/2016

%% Load Data
directory = 'F:\\Imagery\\HyTES\\2016_01_25 Huntington Gardens\\'; %Home directory
dirOut = 'F:\\Dropbox\\Analysis\\JPL Analysis\\HyTES spectra\\Figures\\'; %Output directory for figures
fileName = '2016-01-25_HuntingtonGardens_LST_GeoRef'; %image name
I = enviread(strcat(directory,fileName));

%% Displaying 
figure('units','normalized','outerposition',[0 0 0.5 1])
hold on
lowestValue = min(I.z(I.z(:)>0));
highestValue = 330; %max(I.z(:));
imagesc(I.x,I.y,I.z)

%Adjust color
cmap = jet(256);
colormap(cmap);
caxis(gca,[lowestValue-2/256, highestValue]);
cmap(1,:)=[0,0,0];% Make less than lowest value black:
colormap(cmap)
c = colorbar;
ylabel(c,'Kelvins')

%Other Graph Features
set(gca,'FontSize',24,'FontName','Cambria')
set(c,'FontSize',24,'FontName','Cambria')
set(gca,'XLim',[min(I.x) max(I.x)],'XTick','')
set(gca,'YLim',[min(I.y) max(I.y)],'YTick','')
hold off

%Save Image
%set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 5 6]) %
nameFile = strcat(dirOut,'HyTES_Temperature');
print(nameFile,'-dpng','-r0')