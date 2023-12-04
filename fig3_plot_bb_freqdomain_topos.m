% plot freq. domain topographies for each harmonic
% bb
% xy @ fudan
%% category responses
clear;clc;close all;
bbpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
cd(bbpath)
load('New3Groups_arranged_EEGData_oddball.mat');%NONFILTERED DATA
splinefile = 'EGI_128.spl';
Fs = 420;
T = 1/420;
LenSignal = 490*2;
t = (0:LenSignal-1)*T;
f = (Fs*(0:(LenSignal/2))/LenSignal);
conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Characters'};
freqlists = {'0.8572 Hz' '1.7144 Hz' '2.5716 Hz' '3.4288 Hz'};
AllChansInfo = [];
GroupAmp = [];
allpvals = [];

for con = 1:5
    if con == 1
        thresh = 2;
    else
        thresh = 1;
    end
    for run = 1:4
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
        minEpochDurationSamples = 490*2;
        resampledDataCell =  cellfun(@(x) reshapeTrialToEpochs(x, minEpochDurationSamples), mnewData, ...
            'uni', false);
        mvoltData = cellfun(@(x) x*10^6, resampledDataCell,'uni',false);
        mEEGData = cellfun(@(x) squeeze(nanmean(nanmean(x,4),2)),mvoltData,'uni',false);

        % counting number of valid trials (in 1 second)
        % check if we have NaNs in the data, and we should have 0
        dataNanCounts = cellfun(@(x) sum(isnan(x(:))), mEEGData, 'Uni',false)
        size(mEEGData{1})
        fftData = cellfun(@(x) fft(x),mEEGData,'uni',false);
        offtDataAmp = cellfun(@(x) x/LenSignal,fftData,'uni',false);
        groupData = combineCells(offtDataAmp,1,1);
        absData = cellfun(@(x) abs(x),groupData,'uni',false);

        %start stats analysis
        oddb_ix = [3,5];
        conData = absData(con);
        conDataFreq = cellfun(@(x) x(oddb_ix,:),conData,'uni',false);
        conDataMatrix = cell2mat(conDataFreq);
        %         cleanConData = conDataMarix;
        %         cleanConData(:,chan2go) = 0;
        for ii = 1:length(oddb_ix)
            figure('color','w')
            headplot_xy_oddball([0 0 0 conDataMatrix(ii,:)], splinefile, 'maplimits', [-thresh thresh], ...
                'colormap', jmaColors('coolhotcortex'), 'view', [0 30]);

            XimPrintsize = 1.6;
            YimPrintsize = 2;
            set(gcf,'PaperUnits','centimeters');
            set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
            figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
            figname = sprintf('Topos_bb_freqDomain_oddball_harm%d_group%d_%s.eps',ii,run,conlists{con});
            print(gcf,fullfile(figpath,figname),'-depsc','-r300');
            close
        end
    end
end %end of conditions
%% visual responses

clear;clc;close all;
bbpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
cd(bbpath)
set(0,'DefaultAxesFontSize',8,...
    'defaultTextFontName','Arial',...
    'defaultAxesFontName','Arial');
load('New3Groups_arranged_EEGData_oddball.mat');%NONFILTERED DATA, SO SHOULD NOT CALLED ODDBALL
splinefile = 'EGI_128.spl';
Fs = 420;
T = 1/420;
LenSignal = 490*2;
t = (0:LenSignal-1)*T;
f = (Fs*(0:(LenSignal/2))/LenSignal);
conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Characters'};
freqlists = {'4.286 Hz' '8.572 Hz' '12.858 Hz' '17.144 Hz'}
AllChansInfo = [];
GroupAmp = [];
allpvals = [];

for run = 4
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
        case 5
            adpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\adult_4hz\';
            cd(adpath)
            load('new_arranged_data_4Hz.mat');
            newData = cellfun(@(x) x(:,:,1:6),EEGData,'uni',false);
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
    conEEGData = combineCells(mEEGData,2,1);
    fftData = cellfun(@(x) fft(x),conEEGData,'uni',false);
    offtDataAmp = cellfun(@(x) x/LenSignal,fftData,'uni',false);
    groupData = combineCells(offtDataAmp,1,1);
    absData = cellfun(@(x) abs(x),groupData,'uni',false);
    
    %start stats analysis
    carri_ix = [11,21,31,41];   
    conDataFreq = cellfun(@(x) x(carri_ix,:),absData,'uni',false);
    conDataMatrix = cell2mat(conDataFreq);
    
    for ii = 1:length(carri_ix)
        figure('color','w')
        if run < 5
            headplot_xy_carrier([0 0 0 conDataMatrix(ii,:)], splinefile, 'maplimits', [-1 1],...
                'colormap', jmaColors('coolhotcortex'), 'view', [0 30]);
        elseif run == 5
            headplot_xy_carrier([0 0 0 conDataMatrix(ii,:)], splinefile, 'maplimits', [-0.5 0.5],...
                'colormap', jmaColors('coolhotcortex'), 'view', [0 30]);
        end
        XimPrintsize = 1.6;
        YimPrintsize = 2;
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
    
        figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
        figname = sprintf('Topos_bb_freqDomain_carrier_harm%d_group%d.eps',ii,run);
        print(gcf,fullfile(figpath,figname),'-depsc','-r600');
        close 
    end
end