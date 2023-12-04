%% time domain 
% bb
% xy @ fudan

%% visual: OVERLAY ALL 5 GROUPS IN ONE PLOT 
clear;clc;close all;
figpath = '/Users/xiaoqian/Desktop/arranged_dataset';
cd(figpath)
LeftOT_bb = [57 58 59 64 65 63 68]; 
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 75 70 76 69 74 82 83 89]; %top chans: 66 72 84 
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
lists = {'Faces','Cars','Corridors','Limbs','Characters'};
mycolor = [254,204,92;
    161,218,180;
    65,182,196;
    37,52,148;
    0,0,0]./255;
xx = 1:490;
x = 0:100:1200;
xpoints = (x*420)/1000;
thresh = 6;
roi = 2
if roi == 1
    rois = OTROI; roiname = 'OTROI'
elseif roi == 2
    rois = Posterior_bb; roiname = 'OCCROI'
end

% figure('position',[100 100 450 350],'color','w')
figure('position',[100 100 600 450],'color','w')
hold on;
for run = 1:4
    load('New3Groups_arrange_no_concat_carrier.mat');
    if run == 1
        newData = groupData1;
        subjList_short = idx1;
        groupname = 'Group1-3-4';
    elseif run == 2
        newData = groupData2;
        subjList_short = idx2;
        groupname = 'Group2-4-6';
    elseif run == 3
        newData = groupData3;
        subjList_short = idx3;
        groupname = 'Group3-6-8';
    elseif run == 4
        newData = groupData4;
        subjList_short = idx4;
        groupname = 'Group4-12-15';
    elseif run == 5
        adpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset adult/arranged_dataset/4Hz';
        load(fullfile(adpath,'concat_arranged_data_carrier_adult4Hz.mat'));
        groupname = 'Group5-Adult';
        newData = cellfun(@(x) nanmean(x,3),cdataOut,'uni',false);
    end
    %average across trials
    mnewData = cellfun(@(x) nanmean(x,3),newData,'uni',false);
    mVdataIn = cellfun(@(x) x*10^6, mnewData, 'uni', false);

    nsubj = size(mVdataIn,1);
    newdata = [];
    malldata = [];
    for n = 1:nsubj
        data = mVdataIn(n,:); %data in 3D matrix
        for i = 1:length(data) %combine conditions
            newdata(:,:,i) = squeeze(nanmean(data{i},3));
        end
        mdata = squeeze(nanmean(newdata,4));
        malldata{n,1} = mdata;
    end

    data = malldata;
    temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
    data_temp = cat(2,temp{:});
    meand = nanmean(data_temp,2);
    nsubj = size(data_temp,2);
    se = std(data_temp,[],2)/sqrt(nsubj-1);
    
    if run < 6
        plot(meand,'color',mycolor(run,:),'linewidth',1.8,'LineStyle','-');
    else
        plot(meand,'color',mycolor(run,:),'linewidth',1.8,'LineStyle','--');
    end

    err = patch([xx fliplr(xx)],[(meand-se)',fliplr([meand+se]')],...
        mycolor(run,:),'FaceAlpha',0.2, 'EdgeColor','none');
    err.Annotation.LegendInformation.IconDisplayStyle = 'off';

%     %find peak in each waveform
%     if roi == 1
%         if roi == 1
%             time2 = 233;
%         elseif roi == 2
%             time2 = 120;
%         end
%         tp2 = 420*time2/1000;
%         peakdata = meand(1:tp2);
%         [findmax_value,findmax_ix] = max(peakdata);
%         pmax = plot([findmax_ix findmax_ix],[-thresh thresh],'color',mycolor(run,:),'linestyle','--','linewidth',1.2);
%         pmax.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     end
end
addtp = 1;
if addtp == 1
    time1 = 0;
    time2 = 233;
    tp1 = 420*time1/1000;
    tp2 = 420*time2/1000;
    hold on;
    p = plot([tp1 tp2],[-thresh+0.5 -thresh+0.5],'linewidth',1.5,'color','k')
    p.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
% legend({'3-4 months','4-6 months','6-8 months','12-15 months','adults'},'fontsize',12)
legend({'3-4 months','4-6 months','6-8 months','12-15 months'},'fontsize',12)
% legend({'Occipito-temporal','Posterior medial'},'fontsize',12)
legend boxoff
xticks(xpoints);
xticklabels(x)
xtickangle(90)
ylim([-thresh thresh])
yticks([-thresh:2:thresh])
xlim([0 490])
xlabel('Time (ms)')
ylabel('Amplitude (µV)')
h = plot([0,490],[0,0],'color',[0.5 0.5 0.5],'linewidth',1,'LineStyle','--');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(gca,'fontsize',16,'linewidth',1.3)
%     title(lists{i},'fontweight','normal')
%     title(sprintf('%d-month',months),'fontweight','normal','fontsize',18)
sgtitle(roiname,'fontsize',18)
figname = ['Timedomain-',roiname,'-overlay-4Groups-Carrier-rmtopchan.png'];
saveas(gcf,fullfile(figpath,figname))

% print(gcf,[figpath,'Timedomain-',roiname,'-overlay-4Groups-Carrier-rmtopchan.tiff'],'-dtiff','-r300')


%% CARRIERS: OVERLAY ALL 5 GROUPS IN ONE PLOT  - 5 folds
clear;clc;close all;
figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
cd(figpath)

set(0,'DefaultAxesFontSize',8,...
    'defaultTextFontName','Arial',...
    'defaultAxesFontName','Arial');

LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68,73
RightOT_bb = [90 91 95 96 100 99 94];
Posterior_bb = [71 75 70 76 69 74 82 83 89]; %top chans: 66 72 84 
% Posterior_bb = [66 72 71 76 75 70 69 74 82 83 84 89]; 
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
lists = {'Faces','Cars','Corridors','Limbs','Characters'};

%230,174,53
% mycolor = [230,200,50;
%     161,218,180;
%     65,182,196;
%     37,52,148;
%     0,0,0]./255;
mycolor = [245,210,15;
20,210,95;    
31,144,237;
37,52,148;]./255;

% %plot olds infants first, color order is inversed
% mycolor = [37,52,148;
%     65,182,196;
%     161,218,180;
%     230,174,53]./255;
xx = 1:98;
xxtt = round(xx*1000/420,0);

thresh = 6;
roi = 2
if roi == 1
    rois = OTROI; roiname = 'OTROI'
elseif roi == 2
    rois = Posterior_bb; roiname = 'OCCROI'
end

figure('color','w')
hold on;
for run = 1:4
    load('New3Groups_arrange_no_concat_carrier.mat');
    if run == 1
        newData = groupData1;
        subjList_short = idx1;
        groupname = 'Group1-3-4';
    elseif run == 2
        newData = groupData2;
        subjList_short = idx2;
        groupname = 'Group2-4-6';
    elseif run == 3
        newData = groupData3;
        subjList_short = idx3;
        groupname = 'Group3-6-8';
    elseif run == 4
        newData = groupData4;
        subjList_short = idx4;
        groupname = 'Group4-12-15';
    elseif run == 5
        adpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\adult_4hz\';
        load(fullfile(adpath,'concat_arranged_data_carrier_adult4Hz.mat'));
        groupname = 'Group5-Adult';
        newData = cellfun(@(x) nanmean(x,3),cdataOut,'uni',false);
    end
    %average across trials
    mnewData = cellfun(@(x) nanmean(x,3),newData,'uni',false);
    mVdataIn = cellfun(@(x) x*10^6, mnewData, 'uni', false);

    nsubj = size(mVdataIn,1);
    newdata = [];
    malldata = [];
    for n = 1:nsubj
        data = mVdataIn(n,:); %data in 3D matrix
        for i = 1:length(data) %combine conditions
            newdata(:,:,i) = squeeze(nanmean(data{i},3));
        end
        mdata = squeeze(nanmean(newdata,4));
        malldata{n,1} = mdata;
    end

    data = malldata;
    temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
    data_temp = cat(2,temp{:});
    [timepoints,nsubj] = size(data_temp);
    newdata = reshape(data_temp,[timepoints/5,5,nsubj]);
    new_data_temp = squeeze(nanmean(newdata,2));
    meand = nanmean(new_data_temp,2);
    se = std(new_data_temp,[],2)/sqrt(nsubj-1);

    inData = new_data_temp;
    %use svdnl code for ttest and cluster-based correction
    [realH, realP, realT, corrH, critVal, supraTh, randDist]= ttest_permute_sstats(inData,1000,'mass');
    
    if run < 5
        pp = plot(xxtt,meand,'color',mycolor(run,:),'LineStyle','-','LineWidth',0.8);
        pp.Annotation.LegendInformation.IconDisplayStyle = 'off';
    else
        plot(xxtt,meand,'color',mycolor(run,:),'LineStyle','--');
    end

    err = patch([xxtt fliplr(xxtt)],[(meand-se)',fliplr([meand+se]')],...
        mycolor(run,:),'FaceAlpha',0.2, 'EdgeColor','none');
    err.Annotation.LegendInformation.IconDisplayStyle = 'off';

    b1 = scatter(xxtt(find(corrH)),run*0.2*ones(length(find(corrH)),1),2,mycolor(run,:),'*');
    b1.Annotation.LegendInformation.IconDisplayStyle = 'off';

end

ylim([-thresh thresh])
yticks([-thresh:2:thresh])
xlim([0 233])
xlabel('Time (ms)')
ylabel('Amplitude (µV)')
h = plot([0,233],[0,0],'color',[0.5 0.5 0.5],'linewidth',0.5,'LineStyle','--');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

set(gca,'fontsize',7,'linewidth',0.5,'fontname','arial')
fontsize(gca,7,"points")
fontsize(gcf,7,"points")
fontname(gca,'Arial')
fontname(gcf,'Arial')

%save-skinny
% XimPrintsize = 21/7;
% YimPrintsize = 21/5;
% set(gcf,'PaperUnits','centimeters');
% set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
% figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
% figname = sprintf('Timedomain-%s-overlay-4Groups-Carrier-5folds-skinny.eps',roiname)
% print(gcf,fullfile(figpath,figname),'-depsc','-r600');
 
%save-chunky
XimPrintsize = 21/6;
YimPrintsize = 21/7;
set(gcf,'PaperUnits','centimeters');
set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
figname = sprintf('Timedomain-%s-overlay-4Groups-Carrier-5folds-chunky.eps',roiname)
print(gcf,fullfile(figpath,figname),'-depsc','-r600');

%% ODDBALLS: 4 groups x 5 conditions
clear;clc;close all;
% bbpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset/';
bbpath = '/Users/xiaoqian/Desktop/arranged_dataset';
cd(bbpath)
LeftOT_bb = [50 57 58 59 64 65 63 68]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 101 99 94];
Posterior_bb = [66 72 71 76 75 70 69 74 82 83 84 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

% roiname = 'OTROI'
% rois = OTROI;
% roiname = 'LeftOTROI'
% rois = LeftOT_bb;
roiname = 'RightOTROI'
rois = RightOT_bb;
% roiname = 'OCCROI'
% rois = OCCROI;

Freqlist = [0.8572 4.286];
xlen = 1166.7;
x = 0:100:1200;
xpoints = (x*420)/1000;
sampling = 490;
step = 1/420;
odds = 1/5;
xx = 1:step:sampling/420;
lists = {'Faces','Cars','Corridors','Limbs','Characters'};
mycolor = [1,0,0;...
    0,0.447,0.741;...
    0.466,0.674,0.188;...
    0.929,0.694,0.125;...
    0.1,0.1,0.1];
xx = 1:490;

figure('position',[100 100 1800 450*4],'color','w')
for run = 1:4
    if run == 5
        figure('position',[100 100 1800 450],'color','w')
    end
    if run < 5
        load('New3Groups_arrange_clean_finalbldata.mat');
%         load('New3Groups_arrange_clean_noise_finalbldata.mat');
    end
    if run == 1
        newData = groupData1;
        subjList_short = idx1;
        groupname = 'Group1-3-4';
    elseif run == 2
        newData = groupData2;
        subjList_short = idx2;
        groupname = 'Group2-4-6';
    elseif run == 3
        newData = groupData3;
        subjList_short = idx3;
        groupname = 'Group3-6-8';
    elseif run == 4
        newData = groupData4;
        subjList_short = idx4;
        groupname = 'Group4-12-14';
    elseif run == 5
        adpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset adult/arranged_dataset/4Hz';
        dfile = 'concat_arranged_data_adult4Hz.mat';
        load(fullfile(adpath,dfile));
%         newData = cellfun(@(x) nanmean(x(:,:,[1:46]),3),cdataOut,'uni',false);
%         trialmode = 'bbtrials';
        newData = cellfun(@(x) nanmean(x,3),cdataOut,'uni',false);
        trialmode = 'fulltrials';
        groupname = 'adults';
    end
    mVdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);

   
    for i = 1:5
        if run < 5
            subplot(4,5,i+(run-1)*5);
        else
            subplot(1,5,i)
        end
        hold on;
        data = mVdataIn(:,i);
        temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
        data_temp = cat(2,temp{:});
        %do FDR sig. d
        statsdata = [];
        statsdata(:,1,:) = data_temp;
        stats = Stats_testing_WTAdecoding(statsdata,1,0,'both');
        meand = nanmean(data_temp,2);    
        nsubj = size(data_temp,2);
        se = std(data_temp,[],2)/sqrt(nsubj-1);

        plot(meand,'color',mycolor(i,:),'linewidth',1.8,'linestyle','-');
        err = patch([xx fliplr(xx)],...
            [(meand-se)',fliplr([meand+se]')],...
            mycolor(i,:),'FaceAlpha',0.2, 'EdgeColor','none');
        err.Annotation.LegendInformation.IconDisplayStyle = 'off';
        %sig. time points
        b1 = scatter(xx(find(stats.RES.h_ttest_clustered(1,:))),...
            0.1*ones(length(find(stats.RES.h_ttest_clustered(1,:))),1),20,mycolor(i,:),'*');
        b1.Annotation.LegendInformation.IconDisplayStyle = 'off';
        xticks(xpoints);
        xticklabels(x)
        if run < 5
            thresh = 8; %8 for signal
            yticks([-thresh:1:thresh])
        else
            thresh = 1.5;
            yticks([-thresh:0.5:thresh])
        end
        ylim([-thresh thresh])
        xlim([0 490])
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        h = plot([0,490],[0,0],'color',[0.5 0.5 0.5],'linewidth',1,'LineStyle','--')
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        set(gca,'fontsize',14,'linewidth',1.2)
        title(lists{i},'fontweight','normal')
    end
    sgtitle(roiname,'fontsize',18)
end
if run < 5
    print(gcf,fullfile(bbpath,sprintf('Timedomain-oddball-%s-30Hz-allGroups-FDR.tiff',roiname)),'-dtiff','-r300')
    saveas(gcf,fullfile(bbpath,sprintf('Timedomain-oddball-%s-30Hz-allGroups-FDR.png',roiname)))
% %     % for noise data
%     print(gcf,fullfile(bbpath,sprintf('Timedomain-noise-oddball-%s-allGroups-FDR.tiff',roiname)),'-dtiff','-r300')
%     saveas(gcf,fullfile(bbpath,sprintf('Timedomain-noise-oddball-%s-allGroups-FDR.png',roiname)))
else
    print(gcf,fullfile(adpath,sprintf('Timedomain-oddball-%s-adultGroup-%s-FDR.tiff',roiname,trialmode)),'-dtiff','-r300')
    saveas(gcf,fullfile(adpath,sprintf('Timedomain-oddball-%s-adultGroup-%s-FDR.png',roiname,trialmode)))
end
%     print(gcf,[adpath,'Timedomain-OT+MedialROIS0-adult4Hz.tiff'],'-dtiff','-r300')

%% ODDBALLS: 4 groups x 5 conditions - OVERLAY
clear;clc;close all;
% bbpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset/';
bbpath = '/Users/xiaoqian/Desktop/arranged_dataset';
cd(bbpath)
LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68 % 63 % 50 on top
RightOT_bb = [90 91 95 96 100 99 94]; %101 on top
Posterior_bb = [71 76 75 70 69 74 82 83 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

% roiname = 'OTROI'
% rois = OTROI;
% roiname = 'OCCROI'
% rois = OCCROI;
% roiname = 'LeftOTROI'
% rois = LeftOT_bb;
roiname = 'RightOTROI'
rois = RightOT_bb;

Freqlist = [0.8572 4.286];
xlen = 1166.7;
x = 0:100:1200;
xpoints = (x*420)/1000;
sampling = 490;
step = 1/420;
odds = 1/5;
xx = 1:step:sampling/420;
lists = {'Faces','Cars','Corridors','Limbs','Characters'};

mycolor = [254,204,92;
    161,218,180;
    65,182,196;
    37,52,148]./255;

xx = 1:490;
load('New3Groups_arrange_clean_finalbldata.mat');

figure('position',[100 100 2000 500],'color','w')
for i = 4 %CONDITION HERE
    subplot_tight(1,5,i,[0.2 0.04])
    for run = 1:4
        if run == 5
            figure('position',[100 100 1800 450],'color','w')
        end
        if run == 1
            newData = groupData1;
            subjList_short = idx1;
            groupname = 'Group1-3-4';
        elseif run == 2
            newData = groupData2;
            subjList_short = idx2;
            groupname = 'Group2-4-6';
        elseif run == 3
            newData = groupData3;
            subjList_short = idx3;
            groupname = 'Group3-6-8';
        elseif run == 4
            newData = groupData4;
            subjList_short = idx4;
            groupname = 'Group4-12-14';
        elseif run == 5
            adpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset adult/arranged_dataset/4Hz';
            dfile = 'concat_arranged_data_adult4Hz.mat';
            load(fullfile(adpath,dfile));
            %         newData = cellfun(@(x) nanmean(x(:,:,[1:46]),3),cdataOut,'uni',false);
            %         trialmode = 'bbtrials';
            newData = cellfun(@(x) nanmean(x,3),cdataOut,'uni',false);
            trialmode = 'fulltrials';
            groupname = 'adults';
        end
        mVdataIn = [];
        mVdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);

        hold on;
        data = mVdataIn(:,i);
        temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
        data_temp = cat(2,temp{:});
        %do FDR sig. d
%         statsdata = [];
%         statsdata(:,1,:) = data_temp;
%         stats = Stats_testing_WTAdecoding(statsdata,1,0,'both');

        %use svdnl code for ttest and cluster-based correction
        inData = data_temp;
        [realH, realP, realT, corrH, critVal, supraTh, randDist]= ttest_permute_sstats(inData,1000,'mass');

        meand = nanmean(data_temp,2);
        nsubj = size(data_temp,2);
        se = std(data_temp,[],2)/sqrt(nsubj-1);

        plot(meand,'color',mycolor(run,:),'linewidth',1.8,'linestyle','-');
        err = patch([xx fliplr(xx)],...
            [(meand-se)',fliplr([meand+se]')],...
            mycolor(run,:),'FaceAlpha',0.2, 'EdgeColor','none');
        err.Annotation.LegendInformation.IconDisplayStyle = 'off';
        %sig. time points
        %  b1 = scatter(xx(find(stats.RES.h_ttest_clustered(1,:))),...
        %      0.1*ones(length(find(stats.RES.h_ttest_clustered(1,:))),1),20,mycolor(i,:),'*');
        hold on;
        b1 = scatter(xx(find(corrH)),0.2*run*ones(length(find(corrH)),1),5,mycolor(run,:),'*');
        b1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    xticks(xpoints);
    xticklabels(x)
    if run < 5
        thresh = 8; %8 for signal
        yticks([-thresh:2:thresh])
    else
        thresh = 1.5;
        yticks([-thresh:0.5:thresh])
    end
    ylim([-thresh thresh])
    xlim([0 490])
    xlabel('Time (ms)')
    ylabel('Amplitude (µV)')
    xtickangle(90)
    h = plot([0,490],[0,0],'color',[0.5 0.5 0.5],'linewidth',1,'LineStyle','--');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    l = legend({'3-4 months' '4-6 months' '6-8 months' '12-15 months'},'fontsize',12,'location','southwest');
    legend boxoff
    set(gca,'fontsize',14,'linewidth',1.2)
    title(lists{i},'fontweight','normal','fontsize',18)
end  %of condition
sgtitle(roiname,'fontsize',18)

if run < 5
%     print(gcf,fullfile(bbpath,sprintf('Timedomain-oddball-%s-30Hz-allGroups-FDR.tiff',roiname)),'-dtiff','-r300')
%     saveas(gcf,fullfile(bbpath,sprintf('Timedomain-oddball-%s-30Hz-allGroups-cluster-withEdgeChan.png',roiname)))
    saveas(gcf,fullfile(bbpath,sprintf('Timedomain-oddball-%s-30Hz-allGroups-cluster-rmTopChan.png',roiname)))
% %     % for noise data
%     print(gcf,fullfile(bbpath,sprintf('Timedomain-noise-oddball-%s-allGroups-FDR.tiff',roiname)),'-dtiff','-r300')
%     saveas(gcf,fullfile(bbpath,sprintf('Timedomain-noise-oddball-%s-allGroups-FDR.png',roiname)))
else
%     print(gcf,fullfile(adpath,sprintf('Timedomain-oddball-%s-adultGroup-%s-FDR.tiff',roiname,trialmode)),'-dtiff','-r300')
    saveas(gcf,fullfile(adpath,sprintf('Timedomain-oddball-%s-adultGroup-%s-FDR-cluster0.png',roiname,trialmode)))
end


%% CARRIERS: INDIVIDUAL RESPONSES - subplot
clear;clc;close all;
% addpath(genpath('/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/code2021'))
% bbpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset';
% bbpath = '/Users/xiaoqian/Desktop/arranged_dataset';
bbpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
cd(bbpath)

LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 76 75 70 69 74 82 83 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

% roiname = 'OTROI'
% rois = OTROI;
roiname = 'OCCROI'
rois = OCCROI;

Freqlist = [0.8572 4.286];
xlen = 1166.7;
x = 0:100:1200;
xpoints = (x*420)/1000;
sampling = 490;
step = 1/420;
odds = 1/5;
xx = 1:step:sampling/420;
lists = {'Faces','Cars','Corridors','Limbs','Characters'};
mycolor = [1,0,0;...
    0,0.447,0.741;...
    0.466,0.674,0.188;...
    0.929,0.694,0.125;...
    0.1,0.1,0.1];
xx = 1:98;
xxtt = round(xx*1000/420,0);

load('New3Groups_arrange_no_concat_carrier.mat');

allsubj_info_amp = [];
allsubj_info_latency = [];

for run = 1:4
    if run == 1
        newData = groupData1;
        subjList_short = idx1;
        groupname = 'Group1-3-4';
        figure('position',[10 10 1800 300*5],'color','w')
    elseif run == 2
        newData = groupData2;
        subjList_short = idx2;
        groupname = 'Group2-4-6';
        figure('position',[10 10 1800 300*3],'color','w')
    elseif run == 3
        newData = groupData3;
        subjList_short = idx3;
        groupname = 'Group3-6-8';
        figure('position',[10 10 1800 300*3],'color','w')
    elseif run == 4
        newData = groupData4;
        subjList_short = idx4;
        groupname = 'Group4-12-14';
        figure('position',[10 10 1800 300*3],'color','w')
    end

    mnewData = cellfun(@(x) squeeze(nanmean(x,3)), newData,'uni',false);
    mVdataIn = cellfun(@(x) x*10^6, mnewData, 'uni', false);

    for i = 1:size(mVdataIn,1)
        if run == 1
        subplot_tight(4,5,i)
        elseif run == 4
            subplot_tight(3,5,i)
        else
            subplot_tight(3,5,i)
        end
        hold on;
        data = mVdataIn(i);
        temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
        newtemp = reshape(temp{:},[490/5,5]);
        mtemp = nanmean(newtemp,2);
        plot(xxtt,mtemp,'color','k','linewidth',1.8,'linestyle','-');
        thresh = 12;
        yticks([-thresh:2:thresh])
        ylim([-thresh thresh])
        xlim([0 233])
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        plot([0,233],[0,0],'color',[0.5 0.5 0.5],'linewidth',1,'LineStyle','--');
        plot([65 65],[-thresh thresh],'linewidth',1,'LineStyle','--');
        plot([100 100],[-thresh thresh],'linewidth',1,'LineStyle','--');
        plot([160 160],[-thresh thresh],'linewidth',1,'LineStyle','--');
        set(gca,'fontsize',14,'linewidth',1.2)
        title(strrep(subjList_short{i},'_','-'),'fontweight','normal')
        
        %find peak and latency
        % win 1
        if run < 4
            [mval,mlatency] = timeDomain_peak_metrics(mtemp,60,90,'neg'); %65-100
        else
            [mval,mlatency] = timeDomain_peak_metrics(mtemp,60,90,'neg'); %65-90
        end
        scatter(mlatency,mval,'filled','d','MarkerFaceColor','b')
        allsubj_info_amp(run,i,1) = mval;
        allsubj_info_latency(run,i,1) = mlatency;
        
        % win 2
%         if run == 1
%             [mval,mlatency] = timeDomain_peak_metrics(mtemp,100,160,'pos');
%         elseif run < 4
%             [mval,mlatency] = timeDomain_peak_metrics(mtemp,90,120,'pos');
%         else
%             [mval,mlatency] = timeDomain_peak_metrics(mtemp,90,110,'pos');
%         end
        if run == 1
            [mval,mlatency] = timeDomain_peak_metrics(mtemp,90,160,'pos');
        else
            [mval,mlatency] = timeDomain_peak_metrics(mtemp,90,110,'pos');
        end
        scatter(mlatency,mval,'filled','d','MarkerFaceColor','r')
        allsubj_info_amp(run,i,2) = mval;
        allsubj_info_latency(run,i,2) = mlatency;

%         %win 3
%         [mval,mlatency] = timeDomain_peak_metrics(mtemp,160,233,'pos');
%         scatter(mlatency,mval,'filled','d','MarkerFaceColor','g')
%         allsubj_info_amp(run,i,3) = mval;
%         allsubj_info_latency(run,i,3) = mlatency;

    end
    sgtitle(groupname,'fontsize',18)
    figname = sprintf('TimeDomain_individual_carrier_%s_%s_subplots.png',groupname,roiname);
    saveas(gcf,fullfile(bbpath,figname));
end

% save window metrics
filename = sprintf('TimeDomain_individual_carrier_win_stats_%s_allgroups.mat',roiname);
save(fullfile(bbpath,filename),'allsubj_info_amp','allsubj_info_latency');

%% CARRIER - INDIVIDUAL -OVERLAY
clear;clc;close all;

% bbpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset';
bbpath = '/Users/xiaoqian/Desktop/arranged_dataset';
cd(bbpath)

LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 76 75 70 69 74 82 83 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

% roiname = 'OTROI'
% rois = OTROI;
roiname = 'OCCROI'
rois = OCCROI;

Freqlist = [0.8572 4.286];
xlen = 1166.7;
x = 0:100:1200;
xpoints = (x*420)/1000;
sampling = 490;
step = 1/420;
odds = 1/5;
xx = 1:step:sampling/420;
lists = {'Faces','Cars','Corridors','Limbs','Characters'};
mycolor = [1,0,0;...
    0,0.447,0.741;...
    0.466,0.674,0.188;...
    0.929,0.694,0.125;...
    0.1,0.1,0.1];
xx = 1:490/5;
ttxx = round(1000*xx/420,0)

load('New3Groups_arrange_no_concat_carrier.mat');
for run = 1:4
    if run == 1
        newData = groupData1;
        subjList_short = idx1;
        groupname = 'Group1-3-4';
    elseif run == 2
        newData = groupData2;
        subjList_short = idx2;
        groupname = 'Group2-4-6';
    elseif run == 3
        newData = groupData3;
        subjList_short = idx3;
        groupname = 'Group3-6-8';
    elseif run == 4
        newData = groupData4;
        subjList_short = idx4;
        groupname = 'Group4-12-14';

    end
    
    figure('position',[100 100 800 300*2],'color','w')
    mnewData = cellfun(@(x) squeeze(nanmean(x,3)), newData,'uni',false);
    mVdataIn = cellfun(@(x) x*10^6, mnewData, 'uni', false);
    nsubj = size(mVdataIn,1);
    colors = jet(17);

    for i = 1:nsubj
        hold on;
        data = mVdataIn(i);
        temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
        newtemp = reshape(temp{:},[490/5,5]);
        mtemp = nanmean(newtemp,2);
        plot(ttxx,mtemp,'color',colors(i,:),'linewidth',1.8,'linestyle','-');
    end
 
        thresh = 10;
        yticks([-thresh:2:thresh])
        ylim([-thresh thresh])
        xlim([0 233])
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        plot([0,490],[0,0],'color',[0.5 0.5 0.5],'linewidth',1,'LineStyle','--');
        plot([60 60],[-thresh thresh],'linewidth',1,'LineStyle','--','color',[0.5 0.5 0.5]);
        if run < 3
            plot([90 90],[-thresh thresh],'linewidth',1,'LineStyle','--','color',[0.5 0.5 0.5]);
        else
            plot([80 80],[-thresh thresh],'linewidth',1,'LineStyle','--','color',[0.5 0.5 0.5]);
        end
        plot([120 120],[-thresh thresh],'linewidth',1,'LineStyle','--','color',[0.5 0.5 0.5]);
        
        
        set(gca,'fontsize',14,'linewidth',1.2)
        title(groupname,'fontweight','normal')
    
    figname = sprintf('TimeDomain_individual_overlay_carrier_%s_%s.png',groupname,roiname);
    saveas(gcf,fullfile(bbpath,figname));
end

%% ODDBALL - INDIVIDUAL -OVERLAY
clear;clc;close all;
bbpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset/';
cd(bbpath)
LeftOT_bb = [50 57 58 59 64 65]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 101];
Posterior_bb = [66 72 71 76 75 70 69 74 82 83 84 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

roiname = 'OTROI'
rois = OTROI;
% roiname = 'OCCROI'
% rois = OCCROI;

% roiname = 'OCCROI'
% rois = OCCROI;


Freqlist = [0.8572 4.286];
xlen = 1166.7;
x = 0:100:1200;
xpoints = (x*420)/1000;
sampling = 490;
step = 1/420;
odds = 1/5;
xx = 1:490;
xxtt = round(xx*1000/420,0);

lists = {'Faces','Cars','Corridors','Limbs','Words'};
mycolor = [1,0,0;...
    0,0.447,0.741;...
    0.466,0.674,0.188;...
    0.929,0.694,0.125;...
    0.1,0.1,0.1];
xx = 1:490;

for run = 1
    if run == 5
        figure('position',[100 100 1800 450],'color','w')
    end
    if run < 5
        figure('position',[100 100 800 600],'color','w')
        load('New3Groups_arrange_clean_finalbldata.mat');
%         load('New3Groups_arrange_clean_noise_finalbldata.mat');
    end
    if run == 1
        newData = groupData1;
        subjList_short = idx1;
        groupname = 'Group1-3-4';
    elseif run == 2
        newData = groupData2;
        subjList_short = idx2;
        groupname = 'Group2-4-6';
    elseif run == 3
        newData = groupData3;
        subjList_short = idx3;
        groupname = 'Group3-6-8';
    elseif run == 4
        newData = groupData4;
        subjList_short = idx4;
        groupname = 'Group4-12-14';
    elseif run == 5
        adpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset adult/arranged_dataset/4Hz';
        dfile = 'concat_arranged_data_adult4Hz.mat';
        load(fullfile(adpath,dfile));
%         newData = cellfun(@(x) nanmean(x(:,:,[1:46]),3),cdataOut,'uni',false);
%         trialmode = 'bbtrials';
        newData = cellfun(@(x) nanmean(x,3),cdataOut,'uni',false);
        trialmode = 'fulltrials';
        groupname = 'adults';
    end
    mVdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);
    colors = jet(17);
    for i = 1
        hold on;
        data = mVdataIn(:,i);
        temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
        nsubj = size(temp,1);
        for n = 1:nsubj
            plot(xxtt,temp{n},'color',colors(n,:),'linewidth',1.8,'linestyle','-');
        end

        if run < 5
            thresh = 12; %8 for signal
            yticks([-thresh:2:thresh])
        else
            thresh = 1.5;
            yticks([-thresh:0.5:thresh])
        end
        ylim([-thresh thresh])
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        h = plot([0,490],[0,0],'color',[0.5 0.5 0.5],'linewidth',1,'LineStyle','--')
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        set(gca,'fontsize',14,'linewidth',1.2)
        title([lists{i} ' ' groupname],'fontweight','normal')
    end
    

if run < 5
%     print(gcf,fullfile(bbpath,sprintf('Timedomain-oddball-%s-30Hz-%s-overlay.tiff',roiname,groupname)),'-dtiff','-r300')
    saveas(gcf,fullfile(bbpath,sprintf('Timedomain-oddball-%s-30Hz-%s-%s-overlay.png',roiname,lists{i},groupname)))
else
%     print(gcf,fullfile(adpath,sprintf('Timedomain-oddball-%s-adultGroup-%s-FDR.tiff',roiname,trialmode)),'-dtiff','-r300')
    saveas(gcf,fullfile(adpath,sprintf('Timedomain-oddball-%s-adultGroup-%s-FDR.png',roiname,trialmode)))
end

end

%% ODDBALL - INDIVIDUAL -subplots
clear;clc;close all;
% bbpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset/';
bbpath = '/Users/xiaoqian/Desktop/arranged_dataset';
cd(bbpath)
LeftOT_bb = [58 59 64 65 63 68]; %removed 50 & 57
RightOT_bb = [90 91 95 96 94 99]; %removed 101 & 100
Posterior_bb = [71 76 75 70 69 74 82 83 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

% roiname = 'OTROI'
% rois = OTROI;
for roi = 2
    if roi == 1
        roiname = 'LeftOTROI'
        rois = LeftOT_bb;
    else
        roiname = 'RightOTROI'
        rois = RightOT_bb;
    end
    Freqlist = [0.8572 4.286];
    xlen = 1166.7;
    x = 0:100:1200;
    xpoints = (x*420)/1000;
    sampling = 490;
    step = 1/420;
    odds = 1/5;
    xx = 1:490;
    xxtt = round(xx*1000/420,0);

    lists = {'Faces','Cars','Corridors','Limbs','Characters'};
    mycolor = [1,0,0;...
        0,0.447,0.741;...
        0.466,0.674,0.188;...
        0.929,0.694,0.125;...
        0.1,0.1,0.1];
    xx = 1:490;

    allsubj_info_amp = [];
    allsubj_info_latency = [];
    for run = 1:4
        if run == 5
            figure('position',[100 100 1800 450],'color','w')
        end
        if run < 5
            load('New3Groups_arrange_clean_finalbldata.mat');
            %         load('New3Groups_arrange_clean_noise_finalbldata.mat');
        end
        if run == 1
            newData = groupData1;
            subjList_short = idx1;
            groupname = 'Group1-3-4';
        elseif run == 2
            newData = groupData2;
            subjList_short = idx2;
            groupname = 'Group2-4-6';
        elseif run == 3
            newData = groupData3;
            subjList_short = idx3;
            groupname = 'Group3-6-8';
        elseif run == 4
            newData = groupData4;
            subjList_short = idx4;
            groupname = 'Group4-12-15';
        elseif run == 5
            adpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset adult/arranged_dataset/4Hz';
            dfile = 'concat_arranged_data_adult4Hz.mat';
            load(fullfile(adpath,dfile));
            %         newData = cellfun(@(x) nanmean(x(:,:,[1:46]),3),cdataOut,'uni',false);
            %         trialmode = 'bbtrials';
            newData = cellfun(@(x) nanmean(x,3),cdataOut,'uni',false);
            trialmode = 'fulltrials';
            groupname = 'adults';
        end
        mVdataIn = [];
        mVdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);
        colors = jet(17);

        for i = 1 %condition
            if run == 4
                figure('position',[100 100 1800 350*3],'color','w')
            elseif run < 4
                figure('position',[100 100 1800 350*4],'color','w')
            end
            data = mVdataIn(:,i);
            temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
            nsubj = size(temp,1);
            ff = 0;
            for n = 1:nsubj
                ff = ff + 1;
                subplot_tight(4,4,ff)
                plot(xxtt,temp{n},'color',colors(n,:),'linewidth',1.8,'linestyle','-');
                hold on;

                % measure peak latency and amplitude metrics
                mtemp = temp{n};
                time1 = 400
                time2 = 700
                if i == 1
                    [mval,mlatency] = timeDomain_peak_metrics(mtemp,time1,time2,'pos');
                else
                    [mval,mlatency] = timeDomain_peak_metrics(mtemp,time1,time2,'neg');
                end

                scatter(mlatency,mval,'filled','d','MarkerFaceColor','r')
                allsubj_info_amp(run,i,n) = mval;
                allsubj_info_latency(run,i,n) = mlatency;

                if run == 4
                    thresh = 16; %8 for signal
                    yticks([-thresh:4:thresh])
                elseif run < 4
                    thresh = 14; %8 for signal
                    yticks([-thresh:4:thresh])
                else
                    thresh = 1.5;
                    yticks([-thresh:0.5:thresh])
                end

                xlim([0 1200])
                xticks([0:400:1200])
                ylim([-thresh thresh])
                xlabel('Time (ms)')
                ylabel('Amplitude (µV)')
                h = plot([0,1200],[0,0],'color',[0.5 0.5 0.5],'linewidth',1,'LineStyle','--')
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                set(gca,'fontsize',14,'linewidth',1.2)
                title(strrep(subjList_short{n},'_','-'),'fontweight','normal','fontsize',16)
            end
            sgtitle([groupname '-' lists{i}],'fontsize',18)

            if run < 5
                %     print(gcf,fullfile(bbpath,sprintf('Timedomain-oddball-%s-30Hz-%s-overlay.tiff',roiname,groupname)),'-dtiff','-r300')
                saveas(gcf,fullfile(bbpath,sprintf('Timedomain-oddball-%s-30Hz-%s-%s-subplot-withpeak_%d-%d.png',roiname,lists{i},groupname,time1,time2)))
            else
                %     print(gcf,fullfile(adpath,sprintf('Timedomain-oddball-%s-adultGroup-%s-FDR.tiff',roiname,trialmode)),'-dtiff','-r300')
                saveas(gcf,fullfile(adpath,sprintf('Timedomain-oddball-%s-adultGroup-%s-subplot.png',roiname,trialmode)))
            end
            %         close
        end

    end

    % save window metrics
    filename = sprintf('TimeDomain_individual_oddball_%s_win_stats_%s_allgroups_%d_%d.mat',lists{i},roiname,time1,time2);
    save(fullfile(bbpath,filename),'allsubj_info_amp','allsubj_info_latency');
    disp('done')
end
%%