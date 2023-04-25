%% ODDBALL: INDIVIDUAL HARMONICS
clear;clc;close all;
adpath = '';
cd(adpath);
load('new_arranged_data_4Hz.mat');
LeftOT_bb = [50 57 58 59 64 65]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 101];
Posterior_bb = [66 72 71 76 75 70 69 74 82 83 84 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
parietal_limb = [52 53 60 61 67 62 77 78 85 86 92 72 62];

% roiname = 'OCCROI';
% rois = OCCROI;
% roiname = 'OTROI';
% rois = OTROI;
% roiname = 'RightOTROI';
% rois = RightOT_bb;
roiname = 'LeftOTROI';
rois = LeftOT_bb;
% rois = parietal_limb;
% roiname = 'ParietalROI';

Fs = 420;
T = 1/420;
LenSignal = 490*2;
t = (0:LenSignal-1)*T;
f = (Fs*(0:(LenSignal/2))/LenSignal);
conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Characters'};

AllChansInfo = [];
GroupAmp = [];
allpvals = [];

trialmode = 'bbtrials';
trialnb = 6;
newData = cellfun(@(x) x(:,:,1:trialnb),EEGData,'uni',false);
minEpochDurationSamples = 490*2;
resampledDataCell =  cellfun(@(x) reshapeTrialToEpochs(x, minEpochDurationSamples), newData, ...
    'uni', false);
size(resampledDataCell{1})
mvoltData = cellfun(@(x) x*10^6, resampledDataCell,'uni',false);
mEEGData = cellfun(@(x) squeeze(nanmean(nanmean(x,4),2)),mvoltData,'uni',false);

% counting number of valid trials (in 1 second)
% check if we have NaNs in the data, and we should have 0
dataNanCounts = cellfun(@(x) sum(isnan(x(:))), mEEGData, 'Uni',false)
size(mEEGData{1})
roiData = cellfun(@(x) nanmean(x(:,rois),2),mEEGData,'uni',false);
fftData = cellfun(@(x) fft(x),roiData,'uni',false);
offtDataAmp = cellfun(@(x) x/LenSignal,fftData,'uni',false);
groupData = combineCells(offtDataAmp,1,0);
groupData = cellfun(@(x) squeeze(x), groupData,'uni',false);

%start stats analysis
oddb_ix = [11,21,31,41];
FinalResults = [];
allpvals = [];
allpHs = [];
for con = 1:5
    ii = 0
    for ix = oddb_ix
        ii = ii + 1;
        conData = groupData(con);
        conDataFreq = cellfun(@(x) x(ix,:),conData,'uni',false);
        conDataMarix = cell2mat(conDataFreq);
        output = analysecplx(conDataMarix');
        allpvals(con,ii)= output.pval;
        if output.pval < .05
            allpHs(con,ii) = 1;
        else
            allpHs(con,ii) = 0;
        end
    end
end

%FDR correction
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(allpvals,0.05,'pdep','yes');
adpath_stats = fullfile(adpath,'stats_freq_dbaker_methods');mkdir(adpath_stats);
filename = sprintf('stats_freq_%s_oddball_adult.mat',roiname);
% save(fullfile(adpath_stats,filename),'allpvals','h','crit_p','adj_p');

% plot the pval
conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Characters'};
mycolor = [1,0,0;...
    0,0.447,0.741;...
    0.466,0.674,0.188;...
    0.929,0.694,0.125;...
    0.4,0.4,0.4];
figure('position',[100 100 1800 350],'color','w')
ii = 0;
for con = 1:5
    ii = ii +1;
    subplot(1,5,ii)
    tempval = squeeze(h(con,:));
%     tempval = squeeze(allpHs(con,:));
    b = bar(tempval);
    b.FaceColor = mycolor(con,:);
    b.FaceAlpha = 0.8;
    title(conlists{con},'fontsize',18)
    set(gca,'linewidth',1.2,'fontsize',14);
    xlim([0 10])
    xticks([0:1:10])
%     xticklabels({'' '0.857' '1.714' '2.571' '3.429' '4.286' '5.143' '6' '6.857' '7.714' '8.571'})
    xticklabels({'' '0.857' '1.714' '2.571' '3.429' '5.143' '6' '6.857' '7.714'})
    xtickangle(90)
    ylim([0 1])
    hold on;
    box off
    ylabel('significant (p<.05)')
end

sgtitle(roiname,'fontsize',18)
roinames = [roiname '_carrier4Hs']
figname = sprintf('Fig_stats_adult_freq_tcirctest_significance_oddball_%s.png',roiname);
saveas(gcf,fullfile(adpath_stats,figname))

%% CARRIER - INDIVIDUAL HARMONICS
clear;clc;close all;
adpath = '';
cd(adpath);
load('new_arranged_data_4Hz.mat');

LeftOT_bb = [50 57 58 59 64 65]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 101];
Posterior_bb = [66 72 71 76 75 70 69 74 82 83 84 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
old_central = [7 13 30 31 37 54 55 79 80 87 105 106 112];
central_limb = [7 106 31 55 80 30 37 105 87 54 79];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

% roiname = 'OCCROI';
% rois = OCCROI;
roiname = 'OTROI';
rois = OTROI;
% roiname = 'RightOTROI';
% rois = RightOT_bb;
% roiname = 'LeftOTROI';
% rois = LeftOT_bb;

Fs = 420;
T = 1/420;
LenSignal = 490*2;
t = (0:LenSignal-1)*T;
f = (Fs*(0:(LenSignal/2))/LenSignal);

AllChansInfo = [];
GroupAmp = [];
allpvals = [];

trialmode = 'bbtrials';
trialnb = 6;
newData = cellfun(@(x) x(:,:,1:trialnb),EEGData,'uni',false);
comData = combineCells(newData,2,1);

minEpochDurationSamples = 490*2;
resampledDataCell =  cellfun(@(x) reshapeTrialToEpochs(x, minEpochDurationSamples), comData, ...
    'uni', false);
mvoltData = cellfun(@(x) x*10^6, resampledDataCell,'uni',false);
temp_mEEGData = cellfun(@(x) squeeze(nanmean(x,2)),mvoltData,'uni',false);
mEEGData = cellfun(@(x) squeeze(nanmean(x,3)),temp_mEEGData,'uni',false);

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
        allpvals(ii)= output.pval;
    end
end

%FDR correction
allpvals = squeeze(allpvals)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(allpvals,0.05,'pdep','yes');

bbpath_stats = fullfile(adpath,'stats_adult_freq_dbaker_methods');
filename = sprintf('stats_freq_adult_%s_carri.mat',roiname);
save(fullfile(adpath,filename),'allpvals','h','crit_p','adj_p');

% plot the pval
figure('position',[100 100 400 350],'color','w')
b = bar(h);
b.FaceColor = 'k';
b.FaceAlpha = 0.5;
title(roiname,'fontsize',18)
set(gca,'linewidth',1.2,'fontsize',14);
xlim([0 5])
xticks([0:1:5])
xticklabels({'' '4.286' '8.572' '12.858' '17.144' ''})
xtickangle(90)
ylim([0 1])
xtickangle(90)
box off
ylabel('significant (p<.05)')

figname = sprintf('Fig_stats_adult_freq_tcirctest_significance_carri_%s.png',roiname);
saveas(gcf,fullfile(adpath,figname))


