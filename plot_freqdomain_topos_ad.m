%plot freq. domain topographies for each harmonic
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

oddb_ix = [3,5,7,9];
jj = 0;
figure('position',[100 100 800 400*3],'color','w')
for con = 1:5
    conData = absData(:,con);
    conDataFreq = cellfun(@(x) x(oddb_ix,:),conData,'uni',false);
    conDataMatrix = cell2mat(conDataFreq);
    for ii = 1:length(oddb_ix)
        jj = jj + 1;
        subplot_tight(5,length(oddb_ix),jj,[0.04 0.04])
        headplot([0 0 0 conDataMatrix(ii,:)], splinefile, 'maplimits', [-0.1 0.1], 'title', freqlists{ii},...
            'colormap', jmaColors('coolhotcortex'), 'view', [0 30]);
    end
end

figname = 'Topos_adult_4Hz_freqDomain_perHarmonic_bbtrials_[0.1]_oddball.png';
saveas(gcf,fullfile(adpath,figname));

%% CARRIER RESPONSES
clear;clc;close all;
adpath = '';
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
trialnb = 5;
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
%     chan2go = [17,43,44,48,49,56,63,68,73,81,88,94,99,107,113,114,119];
carri_ix = [11,21,31,41];
figure('position',[100 100 800 300],'color','w')
conDataFreq = cellfun(@(x) x(carri_ix,:),absData,'uni',false);
conDataMatrix = cell2mat(conDataFreq);
for ii = 1:length(carri_ix)
    subplot_tight(1,length(carri_ix),ii,[0.12 0.04])
    headplot([0 0 0 conDataMatrix(ii,:)], splinefile, 'maplimits', [-0.5 0.5], 'title', freqlists{ii},...
        'colormap', jmaColors('coolhotcortex'), 'view', [0 30]);
end

figname = sprintf('Topos_adult_4Hz_freqDomain_perHarmonic_bbtrials_%s.png','carrier');
saveas(gcf,fullfile(adpath,figname));



