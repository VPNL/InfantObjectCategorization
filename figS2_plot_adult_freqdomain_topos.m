%% plot frequency domain topographies for each harmonic
% adult data
% XY @ fudan 2023
%% ODDBALL RESPONSES
clear;clc;close all;
adpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\adult_4hz\';
cd(adpath);
load('new_arranged_data_4Hz.mat');
splinefile = 'EGI_128.spl';
Fs = 420;
T = 1/420;
LenSignal = 490*2;
t = (0:LenSignal-1)*T;
f = (Fs*(0:(LenSignal/2))/LenSignal);
conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Characters'};
freqlists = {'0.8572 Hz' '1.7144 Hz' '2.5716 Hz' '3.4288 Hz'};

trialmode = 'bbtrial';
trialnb = 6;
newData = cellfun(@(x) x(:,:,1:trialnb),EEGData,'uni',false);
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
for con = 1:5
    conData = absData(:,con);
    conDataFreq = cellfun(@(x) x(oddb_ix,:),conData,'uni',false);
    conDataMatrix = cell2mat(conDataFreq);
    for ii = 1:length(oddb_ix)
        figure('color','w')
        headplot_xy_oddball([0 0 0 conDataMatrix(ii,:)], splinefile, 'maplimits', [-0.15 0.15],...
            'colormap', jmaColors('coolhotcortex'), 'view', [0 30]);
        %save
        XimPrintsize = 1.6;
        YimPrintsize = 2;
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
        figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
        figname = sprintf('TOPO_adult4hz_oddb_fft_bbtrials_[0.15]_harmo%d_%s.eps',ii,conlists{con});
        print(gcf,fullfile(figpath,figname),'-depsc','-r600');
    end
end
% close all

%% CARRIER RESPONSES
clear;clc;close all;
adpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\adult_4hz\';
cd(adpath);
load('new_arranged_data_4Hz.mat');
splinefile = 'EGI_128.spl';
Fs = 420;
T = 1/420;
LenSignal = 490*2;
t = (0:LenSignal-1)*T;
f = (Fs*(0:(LenSignal/2))/LenSignal);
conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Characters'};
freqlists = {'4.286 Hz' '8.572 Hz' '12.858 Hz' '17.144 Hz'}

GroupAmp = [];
allpvals = [];
trialmode = 'bbtrial'
trialnb = 6;
newData = cellfun(@(x) x(:,:,1:trialnb),EEGData,'uni',false);
mnewData = cellfun(@(x) squeeze(nanmean(x,3)),newData,'uni',false);
comData = combineCells(mnewData,2,1);

minEpochDurationSamples = 490*2;
resampledDataCell =  cellfun(@(x) reshapeTrialToEpochs(x, minEpochDurationSamples), comData, ...
    'uni', false);
mvoltData = cellfun(@(x) x*10^6, resampledDataCell,'uni',false);
mEEGData = cellfun(@(x) squeeze(nanmean(x,2)),mvoltData,'uni',false);

% counting number of valid trials (in 1 second)
% check if we have NaNs in the data, and we should have 0
dataNanCounts = cellfun(@(x) sum(isnan(x(:))), mEEGData, 'Uni',false)
size(mEEGData{1})

fftData = cellfun(@(x) fft(x),mEEGData,'uni',false);
offtDataAmp = cellfun(@(x) x/LenSignal,fftData,'uni',false);
groupData = combineCells(offtDataAmp,1,1);
absData = cellfun(@(x) abs(x),groupData,'uni',false);

%start stats analysis
carri_ix = [11,21,31,41];

conDataFreq = cellfun(@(x) x(carri_ix,:),absData,'uni',false);
conDataMatrix = cell2mat(conDataFreq);

%use only the first half of the colors in colormap
for ii = 1:length(carri_ix)
    figure('color','w')
    headplot_xy_carrier([0 0 0 conDataMatrix(ii,:)], splinefile, 'maplimits', [-0.4 0.4], ...
        'colormap', jmaColors('coolhotcortex'), 'view', [0 30]);
    %save
    XimPrintsize = 1.6;
    YimPrintsize = 2;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
    figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
    figname = sprintf('Topos_adult_4Hz_freqDomain_Harmonic%d_bbtrials_%s_[0.4].eps',ii,'carrier');
    print(gcf,fullfile(figpath,figname),'-depsc','-r600');
end

close all;


