%% figS1 plot visual acuity measurements from previous studies
% xy

close all;clear;clc;
bbpath = 'D:\Stanford_infant_EEG\';
cd(bbpath)
datafile = 'Norcia-visualacuity-data.csv';
[num,txt,raw] = xlsread(datafile);
close all;
nsample = size(num,2)/2;
mycolors = [1,0,0;
    0,0,1;
    0 0.6588 0.2627;
    0.3010 0.7450 0.9330;
    1,0,1;
    0.5,0.5,0.5;
    0.466,0.674,0.188;
    0.9290 0.6940 0.1250;
    0.494,0.184,0.556];

markerindx = {'square' 'o' 'diamond' 'v' 'o' 'hexagram' 'square' '^' '*'};
markercolors = [1 0 0 ;
    1 1 1;
    0 0.6588 0.2627;
    1 1 1;
    1 0 1;
    1 1 1;
    1 1 1;
    0.9290 0.6940 0.1250;
    0.494,0.184,0.556];

figure('color','w')
hold on;
p = plot([3*30/7 3*30/7],[1 30], 'LineStyle','--','Color',[0.2 0.2 0.2],'LineWidth',0.5);
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
p = plot([6*30/7 6*30/7],[1 30], 'LineStyle','--','Color',[0.2 0.2 0.2],'LineWidth',0.5);
p.Annotation.LegendInformation.IconDisplayStyle = 'off';

for i =1:8 % [1:5,7,8]
    ii = 1 + 2*(i-1);
    xdata = num(:,ii);
    ydata = num(:,ii+1);
    pl = plot(xdata,ydata,'color',mycolors(i,:));
    pl.LineWidth = 0.8;
    pl.Marker = markerindx{i};
    pl.MarkerFaceColor = markercolors(i,:);
    pl.MarkerSize = 4;
    legendnames{i} = txt{2*i};
end
% legendnames(6) = [];
mylegend = legend(legendnames);
mylegend.Location = 'southeast';
mylegend.Box = 'off';

set(gca,'YScale','log')
yticks([1 3 10 30])
ylim([1 30])
xlim([0 75])
ylabel('Grating Acuity (cpd)')
xlabel('Age (weeks)')
box on
set(gca,'linewidth',0.5)
fontsize(gca,8,"points")
        fontsize(gcf,8,"points")
        fontname(gca,'Arial')
        fontname(gcf,'Arial')
        mylegend.FontSize = 6;

XimPrintsize = 11;
    YimPrintsize = 6;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
figname = 'figS1-VisualAcuity-norcia.eps';
print(gcf,'-depsc',fullfile(figpath,figname),'-r600')