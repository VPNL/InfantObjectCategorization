%% ODDBALL: INDIVIDUAL HARMONICS
clear;clc;close all;
% addpath(genpath('/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/code2021'));
% bbpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset';
% cd(bbpath)
% addpath(genpath('/Users/xiaoqian/Desktop/code2021'))
% bbpath = '/Users/xiaoqian/Desktop/arranged_dataset';
bbpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
cd(bbpath)
load('New3Groups_arranged_EEGData_oddball.mat');%NONFILTERED DATA, SO SHOULD NOT CALLED ODDBALL

LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 76 75 70 69 74 82 83 89];
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
allpvals = [];
for roi = 3
    % roiname = 'OCCROI';
    % rois = OCCROI;
    % roiname = 'OTROI';
    % rois = OTROI;
    if roi == 1
        roiname = 'RightOTROI';
        rois = RightOT_bb;
    elseif roi == 2
        roiname = 'LeftOTROI';
        rois = LeftOT_bb;
    else
        roiname = 'OCCROI';
        rois = OCCROI;
    end
    Fs = 420;
    T = 1/420;
    LenSignal = 490*2;
    t = (0:LenSignal-1)*T;
    f = (Fs*(0:(LenSignal/2))/LenSignal);
    conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Characters'};

    for run = 1:4
        switch run
            case 1
                newData = groupData1;
                subjList_short = idx1;
                groupname = '3-4 months';
            case 2
                newData = groupData2;
                subjList_short = idx2;
                groupname = '4-6 months';
            case 3
                newData = groupData3;
                subjList_short = idx3;
                groupname = '6-8 months';
            case 4
                newData = groupData4;
                subjList_short = idx4;
                groupname = '12-15 months';
        end

        mnewData = cellfun(@(x) squeeze(nanmean(x,3)),newData,'uni',false);

        minEpochDurationSamples = 490*2;
        resampledDataCell =  cellfun(@(x) reshapeTrialToEpochs(x, minEpochDurationSamples), mnewData, ...
            'uni', false);
        mvoltData = cellfun(@(x) x*10^6, resampledDataCell,'uni',false);
        mEEGData = cellfun(@(x) squeeze(nanmean(nanmean(x,4),2)),mvoltData,'uni',false);

        % counting number of valid trials (in 1 second)
        % check if we have NaNs in the data, and we should have 0
        dataNanCounts = cellfun(@(x) sum(isnan(x(:))), mEEGData, 'Uni',false)
        size(mEEGData{1})
%         roiData = cellfun(@(x) nanmean(x(:,rois),2),mEEGData,'uni',false);
%         fftData = cellfun(@(x) fft(x),roiData,'uni',false);
%         offtDataAmp = cellfun(@(x) x/LenSignal,fftData,'uni',false);

        fftData = cellfun(@(x) fft(x),mEEGData,'uni',false);
        roiData = cellfun(@(x) nanmean(x(:,rois),2),fftData,'uni',false);
        offtDataAmp = cellfun(@(x) x/LenSignal,roiData,'uni',false);

        %start stats analysis
%         oddb_ix = [3,5,7,9];
        oddb_ix = [3,5];
        FinalResults = [];
        for con = 1:5
            ii = 0
            for ix = oddb_ix
                ii = ii + 1;
                conData = offtDataAmp(:,con);
                conDataFreq = cellfun(@(x) x(ix),conData,'uni',false);
                conDataMarix = cell2mat(conDataFreq);
                output = analysecplx(conDataMarix); % use indiviudal data to caculate the significance of group data
                allpvals(roi,run,con,ii)= output.pval;
            end
        end
    end
end
size(allpvals)
%%
squeeze(allpvals(1,3:4,4,:))
uncorrect_ps = (allpvals<= .05);
figure('units','normalized','position',[0 0 1 1]);
set(gcf,'color','w')
ii = 0
roilists = {'RightOT' 'LeftOT'};
grouplists = {'3-4 months' '4-6 months' '6-8 months' '12-15 months'};
for roi = 1:2
    for run = 1:4
        %FDR correction
        ii = ii + 1;
        subplot_tight(2,4,ii)
        bar(squeeze(uncorrect_ps(roi,run,:,:)))
        xticklabels(conlists)
        xtickangle(90)
        title([roilists{roi} ' - ' grouplists{run}])
    end
end
sgtitle('uncorrected pvals','fontsize',20)

bbpath_stats = fullfile(bbpath,'stats_freq_dbaker_methods');
figname = 'non-FDR_corrected_pvals_allgroups_oddball.png'
saveas(gcf,fullfile(bbpath_stats,figname));

%%
alpha_level = 0.05;
figure('units','normalized','position',[0 0 1 1]);
set(gcf,'color','w')
ii = 0
roilists = {'RightOT' 'LeftOT'};
grouplists = {'3-4 months' '4-6 months' '6-8 months' '12-15 months'};
for roi = 1:2
    for run = 1:4
        %FDR correction
        ii = ii + 1;
        clear h;
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(squeeze(allpvals(roi,run,:,:)),alpha_level,'pdep','yes');
        all_ps(roi,run,:,:) = h;
        subplot_tight(2,4,ii)
        bar(h)
        legend({'0.857 Hz' '1.714 Hz'})
        xticklabels(conlists)
        xtickangle(90)
        title([roilists{roi} ' - ' grouplists{run}])
    end
end
sgtitle('FDR corrected pvals','fontsize',20)

bbpath_stats = fullfile(bbpath,'stats_freq_dbaker_methods');
filename = sprintf('stats_freq_%s_oddball_allgroups_%.2f.mat',roiname,alpha_level);
save(fullfile(bbpath_stats,filename),'allpvals','h','crit_p','adj_p');
figname = 'FDR_corrected_pvals_allgroups_oddball.png'
saveas(gcf,fullfile(bbpath_stats,figname));
%%
% plot the pval
conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Words'};
mycolor = [1,0,0;...
    0,0.447,0.741;...
    0.466,0.674,0.188;...
    0.929,0.694,0.125;...
    0.4,0.4,0.4];
figure('position',[100 100 1800 350*4],'color','w')
ii = 0;
for run = 1:4
    for con = 1:5
        ii = ii +1;
        subplot(4,5,ii)
        tempval = squeeze(h(run,con,:));
        b = bar(tempval);
        b.FaceColor = mycolor(con,:);
        b.FaceAlpha = 0.8;
        title(conlists{con},'fontsize',18)
        set(gca,'linewidth',1.2,'fontsize',14);
        xlim([0 5])
        xticks([0:1:5])
        xticklabels({'' '0.857' '1.714' '2.571' '3.429' '4.286'})
        xtickangle(90)
        ylim([0 1])
        hold on;
        box off
        ylabel('significant (p<.05)')
    end
end
sgtitle(roiname,'fontsize',18)
figname = sprintf('Fig_stats_freq_tcirctest_significance_oddball_%s.png',roiname);
saveas(gcf,fullfile(bbpath_stats,figname))

%% CARRIER - INDIVIDUAL HARMONICS
clear;clc;close all;
% addpath(genpath('/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/code2021'));
% bbpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset';

addpath(genpath('/Users/xiaoqian/Desktop/code2021'))
bbpath = '/Users/xiaoqian/Desktop/arranged_dataset';
cd(bbpath)
load('New3Groups_arranged_EEGData_oddball.mat');%NONFILTERED DATA, SO SHOULD NOT CALLED ODDBALL

LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 76 75 70 69 74 82 83 89];
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;

roiname = 'OCCROI';
rois = OCCROI;
% roiname = 'OTROI';
% rois = OTROI;

Fs = 420;
T = 1/420;
LenSignal = 490*2;
t = (0:LenSignal-1)*T;
f = (Fs*(0:(LenSignal/2))/LenSignal);

AllChansInfo = [];
GroupAmp = [];
allpvals = [];
for run = 1:4
    switch run
        case 1
            newData = groupData1;
            subjList_short = idx1;
            groupname = 'Group1-3-4';
            months = '3-4';
        case 2
            newData = groupData2;
            subjList_short = idx2;
            groupname = 'Group2-4-6';
            months = '4-6';
        case 3
            newData = groupData3;
            subjList_short = idx3;
            groupname = 'Group3-6-8';
            months = '6-8';
        case 4
            newData = groupData4;
            subjList_short = idx4;
            groupname = 'Group4-12-15';
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
    roiData = cellfun(@(x) nanmean(x(:,rois),2),mEEGData,'uni',false);
    fftData = cellfun(@(x) fft(x),roiData,'uni',false);
    offtDataAmp = cellfun(@(x) x/LenSignal,fftData,'uni',false);

    %start stats analysis
    carri_ix = [11,21,31,41];  
    FinalResults = [];
    for con = 1
        ii = 0
        for ix = carri_ix
            ii = ii + 1;
            conData = offtDataAmp(:,con);
            conDataFreq = cellfun(@(x) x(ix),conData,'uni',false);
            conDataMarix = cell2mat(conDataFreq);
            output = analysecplx(conDataMarix); %code from dbaker
            allpvals(run,con,ii)= output.pval;
        end
    end
end
%%
%FDR correction- do 4 levels for each age group
allpvals = squeeze(allpvals)
alpha_level = 0.01;
for run = 1:4
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(allpvals(run,:),alpha_level,'pdep','yes');
    all_hs(run,:) = h;
end
all_hs
bbpath_stats = fullfile(bbpath,'stats_freq_dbaker_methods');
mkdir(bbpath_stats)
filename = sprintf('stats_freq_%s_carri_allgroups_%d.mat',roiname,alpha_level);
save(fullfile(bbpath_stats,filename),'allpvals','h','crit_p','adj_p');
%%
% plot the pval
% to be edited: use log(10)x to convert the pval to 10-x format
groupName = {'3-4months' '4-6months' '6-8months' '12-15months'}
figure('position',[100 100 1400 350],'color','w')
for run = 1:4  
    subplot(1,4,run)
    tempval = squeeze(h(run,:));
    b = bar(tempval);
    b.FaceColor = 'k';
    b.FaceAlpha = 0.5;
    title(groupName{run},'fontsize',18)
    set(gca,'linewidth',1.2,'fontsize',14);
    xlim([0 5])
    xticks([0:1:5])
    xticklabels({'' '4.286' '8.572' '12.858' '17.144' ''})
    xtickangle(90)
    ylim([0 1])
    xtickangle(90)
    box off
    ylabel('significant (p<.01)')  
end
sgtitle(roiname,'fontsize',18)
figname = sprintf('Fig_stats_freq_tcirctest_significance_carri_%s.png',roiname);
saveas(gcf,fullfile(bbpath_stats,figname))