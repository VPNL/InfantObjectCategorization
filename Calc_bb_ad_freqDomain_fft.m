%% SPECTRA - ODDBALL + CARRIER - VECTOR AVERAGE
clear;clc;close all;
bbpath = '';
cd(bbpath)
set(0,'DefaultAxesFontSize',14,...
    'defaultTextFontName','Calibri',...
    'defaultAxesFontName','Calibri');
load('New3Groups_arranged_EEGData_oddball.mat');%NONFILTERED DATA, SO SHOULD NOT CALLED ODDBALL

LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68 % 63, remove 50
RightOT_bb = [90 91 95 96 100 94 99]; % removed 94,99 or remove 101 for the top chan
Posterior_bb = [71 76 75 70 69 74 82 83 89]; %removed 72,66,84
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
% roiname = 'OTROI';
% rois = OTROI;
roiname = 'RightOTROI';
rois = RightOT_bb;
% roiname = 'LeftOTROI';
% rois = LeftOT_bb;

Fs = 420;T = 1/420;LenSignal = 490*2;
t = (0:LenSignal-1)*T;
f = (Fs*(0:(LenSignal/2))/LenSignal);

oddb_ix = [3,5,7,9,13,15,17,19]; nb_oddb = length(oddb_ix);
carri_ix = [11,21,31,41]; nb_carri = length(carri_ix);

AllChansInfo = [];
for ncon = 1:5
    figure('position',[10 10 400 320*5],'color','w')
    for run = 1:5
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
            case 5
                adpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\adult_4hz';
                cd(adpath)
                load("new_arranged_data_4Hz.mat")
                trialnb = 6;
                newData = cellfun(@(x) x(:,:,1:trialnb),EEGData,'uni',false);
        end

        minEpochDurationSamples = 490*2;
        resampledDataCell =  cellfun(@(x) reshapeTrialToEpochs(x, minEpochDurationSamples), newData, ...
            'uni', false);
        mvoltData = cellfun(@(x) x*10^6, resampledDataCell,'uni',false);
        bbEEGData = cellfun(@(x) squeeze(nanmean(x,2)),mvoltData,'uni',false);
        mEEGData = cellfun(@(x) squeeze(nanmean(x,3)),bbEEGData,'uni',false);

        % counting number of valid trials (in 1 second)
        % check if we have NaNs in the data, and we should have 0
        dataNanCounts = cellfun(@(x) sum(isnan(x(:))), mEEGData, 'Uni',false)
        size(mEEGData{1})
        fftData = cellfun(@(x) fft(x),mEEGData,'uni',false);
        roiData = cellfun(@(x) nanmean(x(:,rois),2),fftData,'uni',false);
        fftDataAmp = cellfun(@(x) x/LenSignal,roiData,'uni',false);

        % this section is used to plot the fft spectra
        % getting error bars for complex numbers
        groupData_subj = combineCells(fftDataAmp,1,0);
        groupData_subj = cellfun(@(x) squeeze(x), groupData_subj,'uni',false);

        %estimate ellipse error
        [dA, dF, zSNR, ellipse] = fitErrorEllipse_3D_XY(groupData_subj,[2:41],roiname,groupname,'oddb');
        %plot lillipop
        % plot(errorEllipse(:,1), errorEllipse(:,2),'k-','LineWidth',1)

        groupData = combineCells(fftDataAmp,1,1); %mean of the group
        absData = cellfun(@(x) abs(x),groupData,'uni',false);
        groupData = absData;

        conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Characters'};
        mycolor = [1,0,0;...
            0,0.447,0.741;...
            0.466,0.674,0.188;...
            0.929,0.694,0.125;...
            0.5,0.5,0.5];
        %plot
        nErrLen = size(dA,1);

        for con = ncon%:5
            subplot_tight(5,1,run,[0.1 0.15])
            hold on;
            b = bar(f,groupData{con}(1:LenSignal/2+1),'FaceColor',[1 1 1],'barwidth',0.5);
            b.EdgeColor = [0.5 0.5 0.5];
            b.LineWidth = 1;
            b.FaceColor = 'flat';
            b.CData(oddb_ix,:) = repmat(mycolor(con,:),length(oddb_ix),1);
            b.CData(carri_ix,:) = repmat([0.2 0.2 0.2],length(carri_ix),1);
            errorbar(f(1:nErrLen),groupData{con}(1:nErrLen),dA(:,con,1)', ...
                dA(:,con,2)','k','linestyle','none','CapSize',0,'linewidth',0.8);

            xlim([0.2 8.7])
            xticks(round(f([3:2:21,31:10:61]),3));
            xtickangle(90)
            if run < 5
                ylim([-0.2 2.5])
            elseif run == 5
                ylim([-0.2 0.5])
                yticks([0 0.5])
            end
            ylabel('Amplitude (µV)')
            set(gca,'fontsize',14,'linewidth',1.2)
            babyages = [months,'months'];
        end
    end
    figname = sprintf('bb+ad_fft_2sepochs_%s_%s_%dHarmonics_VectorAVG_pc.tiff',roiname,conlists{ncon},nb_oddb);
    print(gcf,'-dtiff',fullfile(bbpath,figname),'-r300')
end


%% CARRIERS - VECTOR AVERAGING - GROUP LEVEL
clear;clc;close all;
bbpath = '';
cd(bbpath)
set(0,'DefaultAxesFontSize',7,...
    'defaultTextFontName','Helvetica',...
    'defaultAxesFontName','Helvetica');
load('New3Groups_arranged_EEGData_oddball.mat');%NONFILTERED DATA, SO SHOULD NOT CALLED ODDBALL

LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 76 75 70 69 74 82 83 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
old_central = [7 13 30 31 37 54 55 79 80 87 105 106 112];
central_limb = [7 106 31 55 80 30 37 105 87 54 79];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
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
mycolor = [254,204,92;...
    161,218,180;...
    65,182,196;...
    37,52,148]./255;

%plot
figure('color','w')
AllChansInfo = [];
GroupBarMean = [];
GroupBarErr = [];
GroupRSSAmp = [];
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

    subplot_tight(4,1,run,[0.12 0.1])
    con = 1
    hold on;
    b = bar(f,groupData{con}(1:LenSignal/2+1),'FaceColor',[1 1 1],'barwidth',0.5);
    b.EdgeColor = [0.5 0.5 0.5];
    b.LineWidth = 1;
    b.FaceColor = 'flat';
    b.CData(oddb_ix,:) = repmat(mycolor(run,:),length(oddb_ix),1);
    b.CData(carri_ix,:) = repmat([0.2 0.2 0.2],length(carri_ix),1);
    errorbar(f(1:nErrLen),groupData{con}(1:nErrLen),dA(:,con,1)', ...
        dA(:,con,2)','k','linestyle','none','CapSize',0,'linewidth',0.8);

    xlim([0.2 18])
    xticks(round(f([3:2:21,31:10:61]),3));
    xtickangle(90)
    ylim([-0.2 1.5])
    ylabel('Amplitude (µV)')
    if run == 4
        xlabel('Frequency (Hz)')
    end
    set(gca,'linewidth',1.2)
    babyages = [months,'months'];

end
% sgtitle(roiname,'fontsize',7)
set(gcf, 'PaperUnits', 'centimeters');
xwidth = 8.8;
ywidth = 15;
set(gcf, 'PaperPosition', [1 1 xwidth ywidth]); % 8.5x11 paper size, half is 4.25x5.5
set(gcf, 'PaperOrientation', 'portrait');
% figname = sprintf('bb_fft_2sepochs_%s_%dHarmonics_VectorAVG_allGroups_carrier.png',roiname,nb_carri);
% saveas(gcf,fullfile(bbpath,figname));
figname = sprintf('bb_fft_2sepochs_%s_%dHarmonics_VectorAVG_allGroups_carrier_pc.eps',roiname,nb_carri);
print(gcf,'-depsc',fullfile(bbpath,figname),'-r600')


