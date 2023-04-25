% plot topographies

clear;close all;clc;

%% ODDBALL TIME WINDOW TOPOS
clear;close all;clc;
bbpath = ('');
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
    timewins = [450+33,500];
    timewinpts = round(timewins*sampling/totallen,0); %time points
    for tt = 1:length(timewinpts)
        temp = timewinpts(tt);
        tprange = [temp-winLenSamp:temp+winLenSamp];
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
for con = 4
    figure('position',[300 300 1800 300*4],'color','w')
    data = squeeze(allRunData(:,con,:,:));

    nwin = size(data,2);
    ww = 0;
    for run = 1:4
        for win = 1:nwin
            ww = ww + 1;
            subplot_tight(4,nwin,ww,[0.08 0.04])
            tdata = squeeze(data(run,win,:));
            %plot topo
            if run == 1
                titlename = sprintf('%s-%s ms',...
                    num2str(round(allTimePoints(run,win,1)*1167/490,0)),...
                    num2str(round(allTimePoints(run,win,nwinpt)*1167/490,0)));
                headplot([0 0 0 tdata'], splinefile, 'maplimits', [-5 5], 'title', titlename,...
                    'colormap', jmaColors('coolhotcortex'), 'electrodes', 'on', 'view', [0 30]);
            else
                headplot([0 0 0 tdata'], splinefile, 'maplimits', [-5 5],...
                    'colormap', jmaColors('coolhotcortex'), 'electrodes', 'on', 'view', [0 30]);
            end
        end
    end
    sgtitle(conlist{con},'fontsize',20)
    %save
    if plotOdd
        figpath = bbpath;
        figname = sprintf('3D-timedomain-definedWins-%d-%dms-topos-%s.png',timewins(1),timewins(end),conlist{con});
    end
    saveas(gcf,fullfile(figpath,figname))
    close
end
%% plotting 3D for adults
for con = 1:5
    figure('position',[300 300 1800 300],'color','w')
    data = squeeze(allRunData(:,con,:,:));
    nwin = size(data,2);
    ww = 0;
    for run = 5
        for win = 1:nwin
            ww = ww + 1;
            subplot_tight(1,nwin,ww,[0.08 0.04])
            tdata = squeeze(data(run,win,:));
            %plot topo
            titlename = sprintf('%s-%s ms',...
                num2str(round(allTimePoints(run,win,1)*1167/490,0)),...
                num2str(round(allTimePoints(run,win,nwinpt)*1167/490,0)));
            headplot([0 0 0 tdata'], splinefile, 'maplimits', [-1 1], 'title', titlename,...
                'colormap', jmaColors('coolhotcortex'), 'electrodes', 'off', 'view', [0 30]);
        end
    end
    sgtitle(conlist{con},'fontsize',20)
    %save
    if plotOdd
        figname = sprintf('3D-timedomain-definedWins-topos-%s-bbtrials.png',conlist{con});
        figpath = adpath;
    end
    saveas(gcf,fullfile(figpath,figname))
    close
end


%% INFANTS - CARRIERS
clear;close all;clc;
bbpath = ('');
cd(bbpath)
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
