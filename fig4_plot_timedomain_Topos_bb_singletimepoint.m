%% plot topographies
% bb + adult
% xy
%% CATEGORY TIME WINDOW TOPOS
clear;close all;clc;
bbpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
cd(bbpath)
splinefile = 'EGI_128.spl';
conlist = {'Faces','Cars','Corridors','Limbs','Characters'};
plotOdd = 1;
% parsing windows
odds = 1/5;
totallen = 1167;
sampling = 490;
winLenSamp = 14/2; % Temporal window length, in samples; % 33.3ms
allRunData = [];
allTimePoints = [];

for run = 5
    if run < 5
        load('New3Groups_arrange_clean_finalbldata.mat');
    end
    if run == 1
        newData = groupData1;
        subjList_short = idx1;
        groupname = '3-4 months';
    elseif run == 2
        newData = groupData2;
        subjList_short = idx2;
        groupname = '4-6 months';
    elseif run == 3
        newData = groupData3;
        subjList_short = idx3;
        groupname = '6-8 months';
    elseif run == 4
        newData = groupData4;
        subjList_short = idx4;
        groupname = '12-15 months';
    elseif run == 5
        adpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\adult_4hz\';
        cd(adpath)
        load('concat_arranged_data_oddball_adult4Hz.mat');
        newData = cellfun(@(x) nanmean(x(:,:,1:60),3),cdataOut,'uni',false);
        groupname = 'adults';
    end
    mdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);
%     timewins = 50:33:350; %ms
    if run < 5
        timewins = 500; %500 for bb and 200 for adult;
    else
        timewins = 217;
    end
    timewinpts = round(timewins*sampling/totallen,0); %time points
    for tt = 1:length(timewinpts)
        temp = timewinpts(tt);
%         tprange = [temp-winLenSamp:temp+winLenSamp];
        tprange = [temp-winLenSamp:temp];
        timepoints(tt,:) = tprange;
    end
    allTimePoints(run,:,:) = timepoints;

    %get group average
    timedata = [];
    for con = 1:5
        tempdata = mdataIn(:,con);
        avgtemp = nanmean(cat(3,tempdata{:}),3);
        for tt = 1:size(timepoints,1)
            temppts = timepoints(tt,:);
            timedata(con,tt,:) = nanmean(avgtemp(temppts,:),1);
        end
    end
    allRunData(run,:,:,:) = timedata; %run x con x nwins x chan
end
nwinpt = size(allTimePoints,3);

%% plotting 3D for infants
close all;
groupnames = {'3-4 months' '4-6 months' '6-8 months' '12-15 months'};
win = 1;
for con = 1:5    
    data = squeeze(allRunData(:,con,:,:));
    nwin = size(data,1);
    for run = 1:4
        %         ww = ww + 1;
        %         subplot_tight(1,nwin,ww,[0.2 0.04])
        figure('position',[300 300 200 250],'color','w')
        tdata = squeeze(data(run,:));
        %plot topo
        headplot_xy_oddball([0 0 0 tdata], splinefile, 'maplimits', [-5 5],...
            'colormap', jmaColors('coolhotcortex'), 'electrodes', 'on', 'view', [0 30]);
%         headplot_xy_oddball([0 0 0 tdata], splinefile, 'maplimits', [-5 5], 'title', groupnames{run},...
%             'colormap', jmaColors('coolhotcortex'), 'electrodes', 'on', 'view', [0 30]);
        titlename = sprintf('%s-%s ms',...
            num2str(round(allTimePoints(run,win,1)*1167/490,0)),...
            num2str(round(allTimePoints(run,win,nwinpt)*1167/490,0)));
%         sgtitle([conlist{con} '-' titlename],'fontsize',20)

        %save
        time1 = num2str(round(allTimePoints(run,win,1)*1167/490,0));
        time2 = num2str(round(allTimePoints(run,win,nwinpt)*1167/490,0));
        if plotOdd
            figpath = bbpath;
            figname = sprintf('3D-timedomain-definedWins-%s-%sms-topos-%s-group%d-pc.tiff',time1,time2,conlist{con},run);
        end
        %     saveas(gcf,fullfile(figpath,figname))
        print(gcf,fullfile(figpath,figname),'-dtiff','-r300')
        close
    end
end
%% plotting 3D for adults
for con = 1:5
    data = squeeze(allRunData(:,con,:,:));
    win = 1;
    run = 5
    tdata = squeeze(data(run,:));
    %plot topo  
    headplot_xy_oddball([0 0 0 tdata], splinefile, 'maplimits', [-1 1], ...
        'colormap', jmaColors('coolhotcortex'), 'electrodes', 'on', 'view', [0 30]);
    %save
    time1 = num2str(round(allTimePoints(run,win,1)*1167/490,0));
    time2 = num2str(round(allTimePoints(run,win,nwinpt)*1167/490,0));   

    XimPrintsize = 1.6;
    YimPrintsize = 2;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
    figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
    figname = sprintf('3D-timedomain-definedWins-%s-%sms-topos-%s-bbtrials.eps',time1,time2,conlist{con});
    print(gcf,fullfile(figpath,figname),'-depsc','-r600');

    close
end

%% INFANTS - CARRIERS - multiple time points
clear;close all;clc;
bbpath = ('/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset/');
cd(bbpath)
% addpath(genpath('/Users/xiaoqian/Documents/letswave7-master/'))
splinefile = 'EGI_128.spl';
plotOdd = 1;
% parsing windows
odds = 1/5;
totallen = 1167;
sampling = 490;
winLenSamp = 14/2; % Temporal window length, in samples; % 33.3ms
allRunData = [];
allTimePoints = [];
timepoints = [];
load('New3Groups_arrange_no_concat_carrier.mat');

for run = 1:4
    if run == 1
        newData = groupData1;
        subjList_short = idx1;
        groupname = 1;
    elseif run == 2
        newData = groupData2;
        subjList_short = idx2;
        groupname = 2;
    elseif run == 3
        newData = groupData3;
        subjList_short = idx3;
        groupname = 3;
    elseif run == 4
        newData = groupData4;
        subjList_short = idx4;
        groupname = 4;
    
    end
    mdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);
    meanData = cellfun(@(x) nanmean(x,3),mdataIn,'uni',false);
    comDataIn = combineCells(meanData,2,1);
%     timewins = 50:33:350; %ms
% time window and merging windows are different from oddball
% use small windows
    timewins = [135,157,176,198];
    timewinpts = round(timewins*sampling/totallen,0); %time points
    for tt = 1:length(timewinpts)
        temp = timewinpts(tt);
        tprange = [temp-winLenSamp/7:temp+winLenSamp/7];
        timepoints(tt,:) = tprange;
    end
    allTimePoints(run,:,:) = timepoints;

    %get group average
    timedata = [];
    tempdata = comDataIn;
    avgtemp = nanmean(cat(3,tempdata{:}),3);
    for tt = 1:size(timepoints,1)
        temppts = round(timepoints(tt,:),0);
        timedata(tt,:) = nanmean(avgtemp(temppts,:),1);
    end
    
    allRunData(run,:,:) = timedata; %run x con x nwins x chan
end
nwinpt = size(allTimePoints,3);

%plot 3D topos
figure('position',[300 300 1800 300*4],'color','w')
data = allRunData;
nwin = size(data,2);
ww = 0;
for run = 1:4
    for win = 1:nwin
        ww = ww + 1;
        subplot_tight(4,nwin,ww,[0.08 0.04])
        tdata = squeeze(data(run,win,:));
        thresh = 2;
        %plot topo
        if run == 1
            titlename = sprintf('%s-%s ms',...
                num2str(round(allTimePoints(run,win,1)*1167/490,0)),...
                num2str(round(allTimePoints(run,win,nwinpt)*1167/490,0)));
            headplot([0 0 0 tdata'], splinefile, 'maplimits', [-thresh thresh], 'title', titlename,...
                'colormap', jmaColors('coolhotcortex'), 'electrodes', 'on', 'view', [0 30]);
        else
            headplot([0 0 0 tdata'], splinefile, 'maplimits', [-thresh thresh],...
                'colormap', jmaColors('coolhotcortex'), 'electrodes', 'on', 'view', [0 30]);
        end
    end
end

%save
if plotOdd
    figpath = bbpath;
    figname = sprintf('3D-timedomain-definedWins-%d-%dms-topos-%s-scale%d.png',timewins(1),timewins(end),'carrier',thresh);
end
saveas(gcf,fullfile(figpath,figname))
% close

%% ADULTS - CARRIERS
clear;close all;clc;
adpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\adult_4hz\';
cd(adpath)
%         load('new_arranged_data_4Hz.mat');
load('concat_arranged_data_carrier_adult4Hz.mat');
newData = cellfun(@(x) nanmean(x(:,:,1:60),3),cdataOut,'uni',false);
groupname = 'adults';
splinefile = 'EGI_128.spl';
% parsing windows
odds = 1/5;
totallen = 1167;
sampling = 490;
winLenSamp = 14/2; % Temporal window length, in samples; % 33.3ms
allRunData = [];
allTimePoints = [];
timepoints = [];

mdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);
meanData = cellfun(@(x) nanmean(x,3),mdataIn,'uni',false);
comDataIn = combineCells(meanData,2,1);

% timewins = 50:33:350; %ms
% time window and merging windows are different from oddball
% use small windows
timewins = 150;
timewinpts = round(timewins*sampling/totallen,0); %time points
for tt = 1:length(timewinpts)
    temp = timewinpts(tt);
    tprange = [temp-winLenSamp/3.5:temp+winLenSamp/3.5];
    timepoints(tt,:) = tprange;
end
run = 5;
allTimePoints(run,:,:) = timepoints;

%get group average
timedata = [];
tempdata = comDataIn;
avgtemp = nanmean(cat(3,tempdata{:}),3);
for tt = 1:size(timepoints,1)
    temppts = round(timepoints(tt,:),0);
    timedata(tt,:) = nanmean(avgtemp(temppts,:),1);
end
allRunData(run,:,:) = timedata; %run x chan
nwinpt = size(allTimePoints,3);
%%
%plot 3D topos
figure('color','w')
data = allRunData;
win = 1;
tdata = squeeze(data(run,win,:));
thresh = 1;
%plot topo
% titlename = sprintf('%s-%s ms',...
%     num2str(round(allTimePoints(run,win,1)*1167/490,0)),...
%     num2str(round(allTimePoints(run,win,nwinpt)*1167/490,0)));
time1 = num2str(round(allTimePoints(run,win,1)*1167/490,0));
time2 = num2str(round(allTimePoints(run,win,nwinpt)*1167/490,0));
headplot_xy_carrier([0 0 0 tdata'], splinefile, 'maplimits', [-thresh thresh],...
    'colormap', jmaColors('coolhotcortex'), 'electrodes', 'on', 'view', [0 30]);
%save

XimPrintsize = 1.6;
YimPrintsize = 2;
set(gcf,'PaperUnits','centimeters');
set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
figname = sprintf('3D-timedomain-definedWins-carrier-%s-%sms-topos-[1]-bbtrials.eps',time1,time2);
print(gcf,fullfile(figpath,figname),'-depsc','-r600');

close

