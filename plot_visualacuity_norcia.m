close all;clear;clc;
set(0,'DefaultAxesFontSize',14,...
    'defaultTextFontName','Calibri',...
    'defaultAxesFontName','Calibri');
bbpath = 'D:\Stanford_infant_EEG\';
cd(bbpath)
datafile = 'Norcia-visualacuity-data.csv';
[num,txt,raw] = xlsread(datafile);
close all;
nsample = size(num,2)/2;
mycolors = [1,0,0;
    0,0,1;
    0.8500 0.3250 0.0980;
    0.3010 0.7450 0.9330;
    1,0,1;
    1,1,1;
    0.466,0.674,0.188;
    0.9290 0.6940 0.1250;
    0.494,0.184,0.556];
markerindx = {'square' 'o' 'diamond' 'v' 'o' 'hexagram' 'square' '^' '*'};
markercolors = [1 0 0 ;
    1 1 1;
    0.8500 0.3250 0.0980;
    1 1 1;
    1 0 1;
    1 1 1;
    1 1 1;
    0.9290 0.6940 0.1250;
    0.494,0.184,0.556];

figure('position',[100 100 600 400],'color','w')
hold on;
p = plot([3*30/7 3*30/7],[1 30], 'LineStyle','--','Color',[0 0 1],'LineWidth',0.8);
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
p = plot([6*30/7 6*30/7],[1 30], 'LineStyle','--','Color',[0.2 0.2 0.2],'LineWidth',0.8);
p.Annotation.LegendInformation.IconDisplayStyle = 'off';

for i = [1:5,7,8]
    ii = 1 + 2*(i-1);
    xdata = num(:,ii);
    ydata = num(:,ii+1);
    pl = plot(xdata,ydata,'color',mycolors(i,:));
    pl.LineWidth = 1.2;
    pl.Marker = markerindx{i};
    pl.MarkerFaceColor = markercolors(i,:);
    legendnames{i} = txt{2*i};
end
legendnames(6) = [];
mylegend = legend(legendnames);
mylegend.Location = 'southeast';
mylegend.Box = 'off';
mylegend.FontSize = 11;
set(gca,'FontSize',14,'linewidth',1)
set(gca,'YScale','log')
yticks([1 3 10 30])
ylim([1 30])
xlim([0 75])
ylabel('Grating Acuity (cpd)')
xlabel('Age (weeks)')
box on
figpath = '';
figname = 'figS1-VisualAcuity-plots.tif';
print(gcf,'-dtiff',fullfile(figpath,figname),'-r300')