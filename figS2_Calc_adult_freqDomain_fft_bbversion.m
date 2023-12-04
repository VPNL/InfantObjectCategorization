%% fft spectra for oddball + carrier, separate categories
% used make the frequency spectra of both visual and categorical responses
% for visual responses at 4.286 Hz, we average across all five conditions
% as we expect no response differences across conditions
% for adult data

%% category responses at 0.857Hz
clear;clc;close all;
adpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\adult_4hz\';
cd(adpath);
set(0,'DefaultAxesFontSize',8,...
    'defaultTextFontName','Arial',...
    'defaultAxesFontName','Arial');
load('new_arranged_data_4Hz.mat');

LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 76 75 70 69 74 82 83 89];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;

conlists = {'Faces' 'Limbs' 'Corridors'  'Characters' 'Cars'};

mycolor = [1,0,0;
    0.929,0.694,0.125;
    0.466,0.674,0.188;
    0.1,0.1,0.1;
    0,0.447,0.741];

Fs = 420;
T = 1/420;
LenSignal = 490*2;
t = (0:LenSignal-1)*T;
f = (Fs*(0:(LenSignal/2))/LenSignal);

oddb_ix = [3,5,7,9,13,15,17,19]; nb_oddb = length(oddb_ix);
carri_ix = [11,21,31,41,51,61]; nb_carri = length(carri_ix);

trialmode = 'bbtrial';
trialnb = 6;
if strcmp(trialmode,'bbtrial')
    newData = cellfun(@(x) x(:,:,1:trialnb),EEGData,'uni',false);
end
permidx = [1,4,3,5,2]; %reorder
newData = newData(:,permidx);
minEpochDurationSamples = 490*2;
resampledDataCell =  cellfun(@(x) reshapeTrialToEpochs(x, minEpochDurationSamples), newData, ...
    'uni', false);
mvoltData = cellfun(@(x) x*10^6, resampledDataCell,'uni',false);
temp_mEEGData = cellfun(@(x) squeeze(nanmean(x,2)),mvoltData,'uni',false);
mEEGData = cellfun(@(x) squeeze(nanmean(x,3)),temp_mEEGData,'uni',false);
% counting number of valid trials (in 1 second)
for roi = 2 %1:2
    if roi == 1
        roiname = 'LeftOTROI';
        rois = LeftOT_bb;
    else
        roiname = 'RightOTROI';
        rois = RightOT_bb;
    end
    dataNanCounts = cellfun(@(x) sum(isnan(x(:))), mEEGData, 'Uni',false)
    size(mEEGData{1})
    roiData = cellfun(@(x) nanmean(x(:,rois),2),mEEGData,'uni',false);
    fftData = cellfun(@(x) fft(x),roiData,'uni',false);

    % this section is used to plot the fft spectra
    % getting error bars for complex numbers
    fftDataAmp = cellfun(@(x) x/LenSignal,fftData,'uni',false);
    groupData_subj = combineCells(fftDataAmp,1,0);
    groupData_subj = cellfun(@(x) squeeze(x), groupData_subj,'uni',false);

    %estimate ellipse error
    [dA, dF, zSNR, ellipse] = fitErrorEllipse_3D_XY(groupData_subj,[2:61],roiname,'adult','oddb');
    %plot lillipop
    %plot(errorEllipse(:,1), errorEllipse(:,2),'k-','LineWidth',1)
    groupData = combineCells(fftDataAmp,1,1); %mean of the group
    absData = cellfun(@(x) abs(x),groupData,'uni',false);
    groupData = absData;

    %plot
    nErrLen = size(dA,1);
    ncon = 5;
    for con = 1:ncon
        %         subplot_tight(1,5,con,[0.22 0.05])
        figure('color','w')
        hold on;
        b = bar(f,groupData{con}(1:LenSignal/2+1),'FaceColor',[1 1 1],'barwidth',0.5);
        b.EdgeColor = [0.5 0.5 0.5];
        b.LineWidth = 0.8;
        b.FaceColor = 'flat';
        b.CData(oddb_ix,:) = repmat(mycolor(con,:),length(oddb_ix),1);
        b.CData(carri_ix,:) = repmat([0.2 0.2 0.2],length(carri_ix),1);
        errorbar(f(1:nErrLen),groupData{con}(1:nErrLen),dA(:,con,1)', ...
            dA(:,con,2)','k','linestyle','none','CapSize',0,'linewidth',0.5);

        xlim([0.2 8.8])
        xticks(round(f([3:2:21,31:10:61]),3));
        xtickangle(90)
        ylim([-0.1 0.5])
        yticks([0:0.1:0.5])
        if con == 1
            ylabel('Amplitude (µV)')
            h = gca;
            h.YAxis.Visible = 'on';
        end
%         if roi == 2
            xlabel('Frequency (Hz)')
%         end
        %         if roi == 1
        %             title(conlists(con),'fontsize',16)
        %         end
        %         set(gca,'fontsize',16,'linewidth',1.2)
        %         set(gca,'linewidth',1.2)
%         if con > 1
%             h = gca;
%             h.XAxis.Visible = 'off';
%         end

        fontsize(gca,7,"points")
        fontsize(gcf,7,"points")
        fontname(gca,'Arial')
        fontname(gcf,'Arial')

        %save
        XimPrintsize = 21/6;
        YimPrintsize = 21/5;
        if roi == 2
            YimPrintsize = (21/5)+0.5;
        end
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);

        figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
        figname = sprintf('adult_oddb_fft_2sepochs_%s_%dtrials_VectorAVG_%s_v2.eps',roiname,trialnb,conlists{con});
        print(gcf,fullfile(figpath,figname),'-depsc','-r600');
    end
end


%% fft spectra for carrier, merged across categories
clear;clc;close all;
adpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\adult_4hz\';
cd(adpath);
load('new_arranged_data_4Hz.mat');

LeftOT_bb = [63 68 57 58 59 64 65]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 76 75 70 69 74 82 83 89];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;

conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Characters'};
mycolor = [1,0,0;...
    0,0.447,0.741;...
    0.466,0.674,0.188;...
    0.929,0.694,0.125;...
    0.5,0.5,0.5];
Fs = 420;
T = 1/420;
LenSignal = 490*2;
t = (0:LenSignal-1)*T;
f = (Fs*(0:(LenSignal/2))/LenSignal);

oddb_ix = [3,5,7,9,13,15,17,19]; nb_oddb = length(oddb_ix);
carri_ix = [11,21,31,41,51,61]; nb_carri = length(carri_ix);

trialmode = 'bbtrial';
trialnb = 6;

if strcmp(trialmode,'bbtrial')
    newData = cellfun(@(x) x(:,:,1:trialnb),EEGData,'uni',false);
end

mnewData = cellfun(@(x) squeeze(nanmean(x,3)),newData,'uni',false);
comData = combineCells(mnewData,2,0);
minEpochDurationSamples = 490*2;
resampledDataCell =  cellfun(@(x) reshapeTrialToEpochs(x, minEpochDurationSamples), comData, ...
    'uni', false);
mvoltData = cellfun(@(x) x*10^6, resampledDataCell,'uni',false);
temp_mEEGData = cellfun(@(x) squeeze(nanmean(x,2)),mvoltData,'uni',false);
mEEGData = cellfun(@(x) squeeze(nanmean(x,3)),temp_mEEGData,'uni',false);

for roi = 1:2
    if roi == 1
        roiname = 'Occipitotemporal';
        rois = OTROI;
    else
        roiname = 'Occipital';
        rois = OCCROI;
    end
    % counting number of valid trials (in 1 second)
    % check if we have NaNs in the data, and we should have 0
    dataNanCounts = cellfun(@(x) sum(isnan(x(:))), mEEGData, 'Uni',false)
    size(mEEGData{1})
    roiData = cellfun(@(x) nanmean(x(:,rois),2),mEEGData,'uni',false);
    fftData = cellfun(@(x) fft(x),roiData,'uni',false);

    % this section is used to plot the fft spectra
    % getting error bars for complex numbers
    fftDataAmp = cellfun(@(x) x/LenSignal,fftData,'uni',false);
    groupData_subj = combineCells(fftDataAmp,1,0);
    groupData_subj = cellfun(@(x) squeeze(x), groupData_subj,'uni',false);

    %estimate ellipse error
    [dA, dF, zSNR, ellipse] = fitErrorEllipse_3D_XY(groupData_subj,[2:61],roiname,'adult','carri');
    %plot lillipop
    %plot(errorEllipse(:,1), errorEllipse(:,2),'k-','LineWidth',1)

    groupData = combineCells(fftDataAmp,1,1); %mean of the group
    absData = cellfun(@(x) abs(x),groupData,'uni',false);
    groupData = absData;


    %plot
    nErrLen = size(dA,1);
    ncon = 1;
    figure('color','w')
    for con = 1:ncon
        
        hold on;
        b = bar(f,groupData{con}(1:LenSignal/2+1),'FaceColor',[1 1 1],'barwidth',0.5);
        b.EdgeColor = [0.5 0.5 0.5];
        b.LineWidth = 0.5;
        b.FaceColor = 'flat';
        b.CData(oddb_ix,:) = repmat([182, 3, 252]./255,length(oddb_ix),1);
        b.CData(carri_ix,:) = repmat([0.2 0.2 0.2],length(carri_ix),1);
        errorbar(f(1:nErrLen),groupData{con}(1:nErrLen),dA(:,con,1)', ...
            dA(:,con,2)','k','linestyle','none','CapSize',0,'linewidth',0.5);

        xlim([0.2 20])
        xticks(round(f([3:2:21,31:10:61]),3));
        xtickangle(90)
        ylim([-0.05 0.6])
        ylabel('Amplitude (µV)')
        xlabel('Frequency (Hz)')
        set(gca,'linewidth',0.5)
        fontsize(gca,7,"points")
        fontsize(gcf,7,"points")
        fontname(gca,'Arial')
        fontname(gcf,'Arial')

        %         if roi == 2
        %             h = gca;
        %             h.YAxis.Visible = 'off';
        %         end
    end
    %save
    XimPrintsize = 21/4;
    YimPrintsize = 21/6;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);

    figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
    figname = sprintf('adult_fft_2sepochs_%s_%s_%dtrials_VectorAVG.tiff','carrier',roiname,trialnb);
    print(gcf,fullfile(figpath,figname),'-depsc','-r600');

end


