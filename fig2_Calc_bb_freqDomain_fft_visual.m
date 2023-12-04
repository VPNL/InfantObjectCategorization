%% frequency domain 
% bb
% xy @ fudan

%% visual responses - VECTOR AVERAGING - average across all five conditions
clear;clc;close all;
addpath(genpath('D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\code2021\'))
bbpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
cd(bbpath)
load('New3Groups_arranged_EEGData_oddball.mat');%NONFILTERED DATA

LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 76 75 70 69 74 82 83 89];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

roiname = 'OCCROI';
rois = OCCROI;
% roiname = 'OTROI';
% rois = OTROI;
% roiname = 'RightOTROI';
% rois = RightOT_bb;
% roiname = 'LeftOTROI';
% rois = LeftOT_bb;

Fs = 420;
T = 1/420;
LenSignal = 490*2;
t = (0:LenSignal-1)*T;
f = (Fs*(0:(LenSignal/2))/LenSignal);

oddb_ix = [3,5,7,9,13,15,17,19]; nb_oddb = length(oddb_ix);
carri_ix = [11,21,31,41,51,61]; nb_carri = length(carri_ix);

conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Words'};
% mycolor = [254,204,92;...
%     161,218,180;...
%     65,182,196;...
%     37,52,148]./255;
mycolor = [245,210,15;
20,210,95;    
31,144,237;
37,52,148;]./255;

AllChansInfo = [];
GroupBarMean = [];
GroupBarErr = [];
GroupRSSAmp = [];
GroupNoiseAmp_short = {};
GroupNoiseAmp_long = {};
for run = 1:4
    indivData = [];
    absIndivData =[];
    fftData = [];
    fftDataAmp = [];
    switch run
        case 1
            newData = groupData1;
            subjList_short = idx1;
            groupname = '3-4 months';
            months = '3-4';
        case 2
            newData = groupData2;
            subjList_short = idx2;
            groupname = '4-6 months';
            months = '4-6';
        case 3
            newData = groupData3;
            subjList_short = idx3;
            groupname = '6-8 months';
            months = '6-8';
        case 4
            newData = groupData4;
            subjList_short = idx4;
            groupname = '12-15 months';
            months = '12-15';
    end

    mnewData = cellfun(@(x) squeeze(nanmean(x,3)),newData,'uni',false);
    comData = combineCells(mnewData,2,0);

    minEpochDurationSamples = 490*2;
    resampledDataCell =  cellfun(@(x) reshapeTrialToEpochs(x, minEpochDurationSamples), comData, ...
        'uni', false);
    mvoltData = cellfun(@(x) x*10^6, resampledDataCell,'uni',false);
    mEEGData = cellfun(@(x) squeeze(nanmean(nanmean(x,4),2)),mvoltData,'uni',false);

    % counting number of valid trials (in 1 second)
    % check if we have NaNs in the data, and we should have 0
    dataNanCounts = cellfun(@(x) sum(isnan(x(:))), mEEGData, 'Uni',false)
    size(mEEGData{1})
    % roiData = cellfun(@(x) nanmean(x(:,rois),2),mEEGData,'uni',false);
    % fftData = cellfun(@(x) fft(x),roiData,'uni',false);
    % fftDataAmp = cellfun(@(x) x/LenSignal,fftData,'uni',false);

    fftData = cellfun(@(x) fft(x),mEEGData,'uni',false);
    roiData = cellfun(@(x) nanmean(x(:,rois),2),fftData,'uni',false);
    fftDataAmp = cellfun(@(x) x/LenSignal,roiData,'uni',false);

    % this section is used to plot the fft spectra
    % getting error bars for complex numbers

    groupData_subj = combineCells(fftDataAmp,1,0);
    groupData_subj = cellfun(@(x) squeeze(x), groupData_subj,'uni',false);

    % estimate the ellipse errors
    [dA, dF, zSNR, ellipse] = fitErrorEllipse_3D_XY(groupData_subj,[2:61],roiname,groupname,'Carrier');
    nErrLen = size(dA,1);

    groupData = combineCells(fftDataAmp,1,1); %mean of the group
    absData = cellfun(@(x) abs(x),groupData,'uni',false);
    groupData = absData;

    %calculate noise bin amp based on abs(fft)-groupData
    indivData = combineCells(fftDataAmp,1,0); %mean of the group
    absIndivData = cellfun(@(x) abs(x),indivData,'uni',false);
    noisebin_short = [2:2:11];
    noisebin_long = [2:2:21];
    absIndivData = squeeze(absIndivData{:});

    noiseamp_short = absIndivData(noisebin_short,:);
    noiseamp_long = absIndivData(noisebin_long,:); % bin X nbsubjects
    GroupNoiseAmp_short{run}= noiseamp_short;
    GroupNoiseAmp_long{run} = noiseamp_long;

    %calculate projected amplitudes
    new_carri_ix = [11,21,31,41];
    projectedSubj = projectSubjectAmplitudes_XY(fftDataAmp,new_carri_ix);
    DataAmp = projectedSubj.amp;
    nsubj = size(DataAmp,2);
    runMean = mean(DataAmp,2);
    runErr = std(DataAmp,[],2)/sqrt(nsubj-1);
    %save group bar amplitude info
    GroupBarMean(:,run) = runMean;
    GroupBarErr(:,run) = runErr;
    %rssq
    GroupRSSAmpMean(run) = mean(rssq(DataAmp,1)); %1D is subj nb.
    GroupRSSAmpErr(run) = std(rssq(DataAmp,1))/sqrt(nsubj-1); %1D is subj nb.

    %plot
    figure('color','w')
    con = 1
    hold on;
    b = bar(f,groupData{con}(1:LenSignal/2+1),'FaceColor',[1 1 1],'barwidth',0.5);
    b.EdgeColor = [0.5 0.5 0.5];
    b.LineWidth = 0.5;
    b.FaceColor = 'flat';
    b.CData(oddb_ix,:) = repmat(mycolor(run,:),length(oddb_ix),1);
    b.CData(carri_ix,:) = repmat([0.2 0.2 0.2],length(carri_ix),1);

    % errorbar(f,groupData{con}(1:LenSignal/2+1),error_lower(1:LenSignal/2+1,con)', ...
    %     error_upper(1:LenSignal/2+1,con)','k','linestyle','none','CapSize',0,'linewidth',0.8);
    errorbar(f(1:nErrLen),groupData{con}(1:nErrLen),dA(:,con,1)', ...
        dA(:,con,2)','k','linestyle','none','CapSize',0,'linewidth',0.5);

    xlim([0.2 18])
    xticks(round(f([3:2:21,31:10:61]),3));
    xtickangle(90)
    ylim([-0.2 1.5])
    ylabel('Amplitude (µV)')
    set(gca,'linewidth',0.5,'fontsize',8,'fontname','arial')

    %save
    XimPrintsize = 21/4;
    YimPrintsize = 21/7;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
    figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
    figname = sprintf('bb_carrier_fft_2sepochs_%s_VectorAVG_%s.eps',roiname,groupname);
    print(gcf,fullfile(figpath,figname),'-depsc','-r300');

end
close all
%% plot the noise bin amplitudes
set(0,'DefaultAxesFontSize',8,...
    'defaultTextFontName','Arial',...
    'defaultAxesFontName','Arial');

mygroupcolor = [245,210,15;
20,210,95;    
31,144,237;
37,52,148;]./255;

grouplists = {'3-4 months' '4-6 months' '6-8 months' '12-15 months'};
figure('color','w'); hold on;
for run = 1:4
    temp_short = mean(GroupNoiseAmp_short{run},1);
    temp_long = mean(GroupNoiseAmp_long{run},1);
    mean_short = mean(temp_short);
    mean_long = mean(temp_long);
    sd_short = std(temp_short);
    sd_long = std(temp_long);   
    nsubj = size(GroupNoiseAmp_short{run},2);
    se_short = sd_short/sqrt(nsubj);
    se_long = sd_long/sqrt(nsubj);

    b = bar(run,mean_short);
    err = errorbar(run,mean_short,se_short,'LineStyle','none','Color',mygroupcolor(run,:),'linewidth',0.5);
    b.FaceColor = 'w';
    b.EdgeColor = mygroupcolor(run,:);
    
    b1 = bar(run,mean_long);
    err = errorbar(run,mean_long,se_long,'LineStyle','none','Color',mygroupcolor(run,:),'linewidth',0.5);
    b1.FaceColor = mygroupcolor(run,:);
    b1.FaceAlpha = 0.5;
    b1.EdgeColor = mygroupcolor(run,:);

end

xticks([0:1:5])
xticklabels({'' '3-4 mo' '4-6 mo' '6-8 mo' '12-15 mo' ''})
xtickangle(90)
xlim([0 5])
ylim([0 0.5])
yticks([0:0.1:0.5])
ylabel('Amplitude (µV)')
set(gca,'fontsize',7,'fontname','arial','linewidth',0.5);
fontsize(gca,7,"points")
fontsize(gcf,7,"points")
fontname(gca,'Arial')
fontname(gcf,'Arial')

% %save-skinny
% XimPrintsize = 21/7;
% YimPrintsize = 21/5;
% set(gcf,'PaperUnits','centimeters');
% set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
% figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
% figname = sprintf('Freqdomain-%s-overlay-4Groups-Carrier-noiseAmp-skinny.eps',roiname)
% print(gcf,fullfile(figpath,figname),'-depsc','-r600');

%save-chunky
XimPrintsize = 21/6;
YimPrintsize = 21/7;
set(gcf,'PaperUnits','centimeters');
set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
figname = sprintf('Freqdomain-%s-overlay-4Groups-Carrier-noiseAmp-chunky.eps',roiname)
print(gcf,fullfile(figpath,figname),'-depsc','-r600');

%% calculate the noise bin amplides across groups
load('New3Groups_arrange_bbages.mat');
 %load baby ages
    new_idx1 = cellfun(@(x) x(1:4),idx1,'uni',false);
    new_idx2 = cellfun(@(x) x(1:4),idx2,'uni',false);
    new_idx3 = cellfun(@(x) x(1:4),idx3,'uni',false);
    new_idx4 = cellfun(@(x) x(1:4),idx4,'uni',false);

    %fit lmm and plot lmm
    %model fitting - lmm: fixed intercept and fixed slope
    %get baby labels and mark down longitudinal babies
    comlabel = [new_idx1;new_idx2;new_idx3;new_idx4];
    unilabel = unique(comlabel);
    neworder = [];
    for ll = 1:length(unilabel)
        temlabel = unilabel(ll);
        tempix = strmatch(temlabel,comlabel);
        neworder(tempix) = ll;
    end

    com_Age = [group1;group2;group3;group4];
    com_Amp_short = [mean(GroupNoiseAmp_short{1},1),...
        mean(GroupNoiseAmp_short{2},1),...
        mean(GroupNoiseAmp_short{3},1),...
        mean(GroupNoiseAmp_short{4},1)];
    com_Amp_long = [mean(GroupNoiseAmp_long{1},1),...
        mean(GroupNoiseAmp_long{2},1),...
        mean(GroupNoiseAmp_long{3},1),...
        mean(GroupNoiseAmp_long{4},1)];


tbl1 = table(log10(com_Age),com_Amp_short',neworder','VariableNames',{'Age','NoiseAmp','Baby'});
lme1 = fitlme(tbl1,'NoiseAmp~ Age +(1|Baby)')

tbl2 = table(log10(com_Age),com_Amp_long',neworder','VariableNames',{'Age','NoiseAmp','Baby'});
lme2 = fitlme(tbl2,'NoiseAmp~ Age +(1|Baby)')

figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
filename =  sprintf('FreqDomain_noisebin_%s_%s_AcrossAges_lmm_stats_log10Age.mat','amp',roiname);
save(fullfile(figpath,filename),'lme1','lme2');


