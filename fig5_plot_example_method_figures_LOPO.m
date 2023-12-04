%% plot the example figure for the LOPO method

clear;clc;close all;
bbpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
cd(bbpath)
load('New3Groups_arrange_clean_concat.mat');

run = 3
newData = groupData3;
subjList_short = idx3;
groupname = 'Group3-6-8';

%permute category labels
permidx = [1,4,3,5,2];
newData = newData(:,permidx);
mcdataOut = cellfun(@(x) squeeze(nanmean(x,3)),newData,'Uni',false);
dataIn = mcdataOut; %update here!!!
[nsubj,ncon] = size(dataIn);

% ROI
LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 76 75 70 69 74 82 83 89];
roiname = 'OT+OCC';

% temporal sliding window
Freqlist = [0.8572 4.286]; sampling = 490; odds = 1/5; tduration = 1167; %ms
Fs = sampling/tduration;
winLenSamp = 14; % Temporal window length, in samples; % 17*2 ms
winHopSamp = 14; % Temporal window hop size, in samples; %no overlapping between windows
wintime = winLenSamp*1000/420;
temp = dataIn{1,1};
[nTime, nSpace, nTrial] = size(temp); % Dimensions of input data matrix

if nTime == sampling * 10
    nTime = sampling;
end
nWins = floor((nTime - winLenSamp) / winHopSamp + 1); % # classifications
winLen = round(winLenSamp*tduration/sampling,1);
ncon = 5;
withZ = 1; % normalize or not, do zsocre across conditions for each electrode

%%%%%% LOPO %%%%%%
% get time domain data averaged across 3 ROIs and concatanete them for each
% cateogry.
for n = 1:nsubj
    testix = n;
    trainix = setdiff([1:nsubj],testix);
    train_init = dataIn(trainix,:);
    test_init = dataIn(testix,:);
    train = combineCells(train_init,1,1);
    test = combineCells(test_init,1,1);

    %z-scoring -use zscore_merge.m function
    %normalize amplitudes across conditions for each electrode
    tempTrainZ(n,:) = zscore_merge(train);
    tempTestZ(n,:) = zscore_merge(test);
    %no normalization
    tempTrainNZ(n,:) = train;
    tempTestNZ(n,:) = test;
end

% Classify the data in each time window
dataInX = tempTrainZ;
dataInY = tempTestZ;

% Classify the data in each time window
dataInX = tempTrainZ;
dataInY = tempTestZ;
subj = 5;
%windowing z-scored train and test data
xdata = cellfun(@(x) [nanmean(x(:,LeftOT_bb),2);nanmean(x(:,Posterior_bb),2);nanmean(x(:,RightOT_bb),2)], dataInX(subj,:), 'uni', false);
ydata = cellfun(@(x) [nanmean(x(:,LeftOT_bb),2);nanmean(x(:,Posterior_bb),2);nanmean(x(:,RightOT_bb),2)], dataInY(subj,:), 'uni', false);

%% plot: 
% xdata: N-1
conlists = {'Faces' 'Limbs' 'Corridors'  'Characters' 'Cars'};
for nc = 1:5
    figure('color','w')
    temp_xdata = xdata{nc};
    imagesc(temp_xdata,[-2 2])
    cmap=mrvColorMaps('coolhot'); 
    colormap(cmap);
%     set(gca,'fontsize',8,'fontname','arial')
    set(gca,'xtick',[],'ytick',[])
%     yticks([1,490,490*2,490*3])

    %save
    XimPrintsize = 0.2;
    YimPrintsize = 3;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
    figname = sprintf('example_patternfigure_kgsv_LOPO_N-1_%s.eps',conlists{nc});
    print(gcf,fullfile(figpath,figname),'-depsc','-r300');
end
close all


%% ydata: N

for nc = 1%:5
    temp_ydata = ydata{nc};
    figure('color','w')
    imagesc(temp_ydata,[-2 2])
    cmap=mrvColorMaps('coolhot'); 
    colormap(cmap);
    set(gca,'xtick',[],'ytick',[])
    %save
    XimPrintsize = 0.2;
    YimPrintsize = 3;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
    figname = sprintf('example_patternfigure_kgsv_LOPO_N_%s.eps',conlists{nc});
    print(gcf,fullfile(figpath,figname),'-depsc','-r300');
end

%% 
%make color bar
figure('position',[1200 100 200 100],'color','w')
imagesc(xdata{2},[-2 2])
cmap=mrvColorMaps('coolhot');
colormap(cmap);
cb = colorbar;
cb.FontSize = 14;
cb.Ticks = [-2 -1 0 1 2];
cb.LineWidth = 1.2;
figname = 'example_patternfig_kgs_colorbar.tiff'
print(gcf,fullfile(figpath,figname),'-dtiff','-r300')

%% plot RSM
% face, cars, corri,limb,character
% permute to face, limb,corri,character, and cars, 1,4,3,5,2
load('WTAresults_Distinctiveness_kgsv_LOPO_OT+OCC_6-8 months.mat')
figure('position',[100 100 300 300],'color','w')
ss = 6
CorrValSym = corrmat(:,:,ss);
imagesc(CorrValSym, [-0.5 0.5]);
axis('image');
cmap=mrvColorMaps('coolhot');
colormap(cmap);colorbar;
labels = {'Faces' 'Limbs' 'Corridors'  'Characters' 'Cars'};
set(gca,'Xtick', [1:1:5], 'XtickLabel',labels)
set(gca,'Ytick', [1:1:5], 'YtickLabel',labels)
set(gca,'fontsize',14,'linewidth',1.2)
axis square
xtickangle(90)
xlabel('Training')
ylabel('Testing')
title('RSM','fontsize',16)

figname = sprintf('example_patternfigure_kgsv_LOPO_RSM_subj%dGroup%d.tiff',ss,3)
print(gcf,fullfile(figpath,figname),'-dtiff','-r300')