%% time domain topos
% bb
% xy @ fudan
clear;close all;clc;
% addpath(genpath('D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\code2021\'))

%% category TIME WINDOW TOPOS
clear;close all;clc;
bbpath = ('D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\');
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
load('New3Groups_arrange_clean_finalbldata.mat');
for run = 1:4
    if run == 1
        newData = groupData1;
        subjList_short = idx1;
        groupname = '3-4months';
    elseif run == 2
        newData = groupData2;
        subjList_short = idx2;
        groupname = '4-6months';
    elseif run == 3
        newData = groupData3;
        subjList_short = idx3;
        groupname = '6-8months';
    elseif run == 4
        newData = groupData4;
        subjList_short = idx4;
        groupname = '12-15months';
    
    end
    mdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);
    timewins = [450+33,500]; % 483-500ms time window
    timewinpts = round(timewins*sampling/totallen,0); %time points    
    allTimePoints(run,:) = [timewinpts(1):timewinpts(2)];

    %get group average
    timedata = [];
    for con = 1:5
        tempdata = mdataIn(:,con);
        avgtemp = nanmean(cat(3,tempdata{:}),3);
        temppts = allTimePoints(run,:);
        timedata(con,:) = nanmean(avgtemp(temppts,:),1);
    end
    allRunData(run,:,:) = timedata; %run x con x chan
end


%% plot
for con = 1:5
    data = squeeze(allRunData(:,con,:,:));
    for run = 1:4
    switch run
        case 1
            groupname = '3-4months';
        case 2
            groupname = '4-6months';
        case 3
            groupname = '6-8months';
        case 4
            groupname = '12-15months';
    end
        figure('color','w')
        tdata = squeeze(data(run,:));
        %plot topo
        headplot_xy_oddball([0 0 0 tdata], splinefile, 'maplimits', [-5 5],...
            'colormap', jmaColors('coolhotcortex'), 'electrodes', 'on', 'view', [0 30]);
        %save
        figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
        figname = sprintf('fig4-3DTopo-timedomain-definedWins-%d-%dms-%s-%s.eps',timewins(1),timewins(end),conlist{con},groupname);

        XimPrintsize = 1.6;
        YimPrintsize = 2;
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
        print(gcf,fullfile(figpath,figname),'-depsc','-r300');

        close
    end
end