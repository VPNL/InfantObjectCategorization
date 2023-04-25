%% CARRIERS
clear;clc;close all;
adpath = '';
cd(adpath)
set(0,'DefaultAxesFontSize',16,...
    'defaultTextFontName','Calibri',...
    'defaultAxesFontName','Calibri');
load('concat_arranged_data_carrier_adult4Hz.mat');
groupname = 'Adult';

LeftOT_bb = [63 68 57 58 59 64 65]; %removed chan68,73
RightOT_bb = [90 91 95 96 100 99 94];
Posterior_bb = [71 76 75 70 69 74 82 83 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
lists = {'Faces','Cars','Corridors','Limbs','Characters'};
xx = 1:490;
x = 0:100:1200;
xpoints = (x*420)/1000;
thresh = 1.5;
trialmode = 'bbtrial';
trialnb = 60;
EEGData = cellfun(@(x) x(:,:,1:trialnb),cdataOut,'uni',false);
newData = cellfun(@(x) nanmean(x,3),EEGData,'uni',false);
%reshape
% minEpochDurationSamples = 490;
% resampledDataCell =  cellfun(@(x) reshapeTrialToEpochs(x, minEpochDurationSamples), newData, ...
%     'uni', false);
% %average across trials
% mnewData = cellfun(@(x) squeeze(nanmean(x,2)),resampledDataCell,'uni',false);
mVdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);

%combine conditions
comData = combineCells(mVdataIn,2,1);
nsubj = size(comData,1);
data = combineCells(comData,1,0);
figure('position',[100 100 400*2 300*1],'color','w')
for roi = 1:2
    if roi == 1
        rois = OTROI; roiname = 'Occipitotemporal'
    elseif roi == 2
        rois = Posterior_bb; roiname = 'Occipital'
    end
    roiData = cellfun(@(x) squeeze(nanmean(x(:,rois,:),2)), data, 'uni', false);
    roiData = roiData{:};
    meand = nanmean(roiData,2);
    se = std(roiData,[],2)/sqrt(nsubj-1);
    %stats test - cluster analysis
    [realH, realP, realT, corrH, critVal, supraTh, randDist]= ttest_permute_sstats(roiData,1000,'mass');

    % figure('position',[100 100 450 350],'color','w')
    % figure('position',[100 100 600 450],'color','w')
    subplot_tight(1,2,roi,[0.25 0.1])
    hold on;
    mycolor = [182, 3, 252]./255;
    plot(meand,'color',mycolor,'linewidth',1.8,'LineStyle','-');
    err = patch([xx fliplr(xx)],[(meand-se)',fliplr([meand+se]')],...
        mycolor,'FaceAlpha',0.2, 'EdgeColor','none');
    err.Annotation.LegendInformation.IconDisplayStyle = 'off';
    b1 = scatter(xx(find(corrH)),0.05*ones(length(find(corrH)),1),2,mycolor,'*');
    b1.Annotation.LegendInformation.IconDisplayStyle = 'off';

addtp = 1;
if addtp == 1
    time1 = 0;
    time2 = 233;
    tp1 = 420*time1/1000;
    tp2 = 420*time2/1000;
    hold on;
    p = plot([tp1 tp2],[-thresh+0.1 -thresh+0.1],'linewidth',1.5,'color','k')
    p.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
xticks(xpoints);
xticklabels(x)
xtickangle(90)
ylim([-thresh thresh])
yticks([-thresh:0.5:thresh])
xlim([0 490])
xlabel('Time (ms)')
ylabel('Amplitude (µV)')
h = plot([0,490],[0,0],'color',[0.5 0.5 0.5],'linewidth',1,'LineStyle','--');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(gca,'fontsize',16,'linewidth',1.3)
title(roiname,'fontsize',18)
%     title(sprintf('%d-month',months),'fontweight','normal','fontsize',18)
if roi == 2
    h = gca;
    h.YAxis.Visible = 'off';
end
end
% sgtitle(roiname,'fontsize',18)
print(gcf,[adpath,'/Timedomain-adult-OT+OCC-Carrier-pc.tiff'],'-dtiff','-r300')

%% ODDBALL
clear;clc;close all;
adpath = '';
cd(adpath)
load('concat_arranged_data_oddball_adult4Hz.mat');
groupname = 'Adult';
set(0,'DefaultAxesFontSize',16,...
    'defaultTextFontName','Calibri',...
    'defaultAxesFontName','Calibri');
LeftOT_bb = [57 58 59 64 65 68 73]; %removed chan68,73
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 76 75 70 69 74 82 83 89];

parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
parietal_limb = [52 53 60 61 67 62 77 78 85 86 92 72 62];
lists = {'Faces' 'Limbs' 'Corridors'  'Characters' 'Cars'};
xx = 1:490;
x = 0:100:1200;
xpoints = (x*420)/1000;

mycolor = [1,0,0;
    0.929,0.694,0.125;
    0.466,0.674,0.188;
    0.1,0.1,0.1;
    0,0.447,0.741];

trialmode = 'bbtrial';
trialnb = 60;
EEGData = cellfun(@(x) x(:,:,1:trialnb),cdataOut,'uni',false);
newData = cellfun(@(x) nanmean(x,3),EEGData,'uni',false);
permidx = [1,4,3,5,2];
newData = newData(:,permidx);
% %reshape
% minEpochDurationSamples = 490;
% resampledDataCell =  cellfun(@(x) reshapeTrialToEpochs(x, minEpochDurationSamples), newData, ...
%     'uni', false);
% %average across trials
% mnewData = cellfun(@(x) squeeze(nanmean(x,2)),resampledDataCell,'uni',false);
% mVdataIn = cellfun(@(x) x*10^6, mnewData, 'uni', false);
mVdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);
nsubj = size(EEGData,1);

%combine subjs
data = combineCells(mVdataIn,1,0);
% figure('position',[100 100 450 350],'color','w')
% figure('position',[100 100 500 350],'color','w')
allroiData = {};

for roi = 1:2
    figure('position',[10 10 300*5 300],'color','w')
    if roi == 1
        roiname = 'LeftOTROI';
        rois = LeftOT_bb;
    elseif roi == 2
        roiname = 'RightOTROI';
        rois = RightOT_bb;
    elseif roi == 3
        rois = Posterior_bb; roiname = 'OCCROI'
    else
        rois = parietal_limb;roiname = 'ParietalROI';
    end
    for con = 1:5
        roiData = cellfun(@(x) squeeze(nanmean(x(:,rois,:),2)), data, 'uni', false);
        allroiData{con}(:,:,roi) = roiData{con};
        conData = [];
        conData = roiData{con};
        meand = nanmean(conData,2);
        se = std(conData,[],2)/sqrt(nsubj-1);
        %stats test - cluster analysis
        [realH, realP, realT, corrH, critVal, supraTh, randDist]= ttest_permute_sstats(conData,5000,'mass');
        subplot_tight(1,5,con,[0.25 0.05])
        thresh = 1.5;
        hold on;
        plot(meand,'color',mycolor(con,:),'linewidth',1.8,'LineStyle','-');
        err = patch([xx fliplr(xx)],[(meand-se)',fliplr([meand+se]')],...
            mycolor(con,:),'FaceAlpha',0.2, 'EdgeColor','none');
        err.Annotation.LegendInformation.IconDisplayStyle = 'off';
        b2 = scatter(xx(find(corrH)),-thresh*ones(length(find(corrH)),1),3,mycolor(con,:),'*');
        b2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        addtp = 1;
        if addtp == 1
            time1 = 0;
            time2 = 233;
            tp1 = 420*time1/1000;
            tp2 = 420*time2/1000;
            hold on;
            p = plot([tp1 tp2],[-thresh-0.3 -thresh-0.3],'linewidth',1.5,'color','k')
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        xticks(xpoints);
        xticklabels(x)
        xtickangle(90)
        ylim([-thresh-0.5 thresh])
        yticks([-thresh-0.5:0.5:thresh])
        xlim([0 490])
        if roi == 2
            xlabel('Time (ms)')
        end
        ylabel('Amplitude (µV)')
        if con > 1
            h = gca;
            h.YAxis.Visible = 'off';
        end
        h = plot([0,490],[0,0],'color',[0.5 0.5 0.5],'linewidth',1,'LineStyle','--');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        set(gca,'fontsize',16,'linewidth',1.3)
    end   
    set(gca,'fontsize',16,'linewidth',1.3)
    print(gcf,fullfile(adpath,sprintf('Timedomain_adult_oddball_%s_withDiff_pc_new.tiff',roiname)),'-dtiff','-r300');
end
