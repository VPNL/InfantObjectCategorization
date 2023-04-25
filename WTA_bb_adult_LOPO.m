% 2021 bigideas
% winner-take-all maximum correlation classification: split half of the
% trials 
% this is for oddball response
% data is partly preprocessed beforehand:
% 1. low-pass 30hz
% 2. detrend: ddataOut
%% step1: add paths
clear;clc;close all;


%% step2: read in EEG data and extra analyses
bbpath = '';
cd(bbpath)
load('New3Groups_arrange_clean_concat.mat');
set(0,'DefaultAxesFontSize',14,...
    'defaultTextFontName','Calibri',...
    'defaultAxesFontName','Calibri');
All_corrC = [];
% permute to face, limb,corri,character, and cars, 1,4,3,5,2
for run = 1
    if run == 1
        newData = groupData1;
        subjList_short = idx1;
        groupname = '3-4 months';
    elseif run == 2
        newData = groupData2;
        subjList_short = idx2;
        groupname = '4-6 months';
    elseif run == 3
        newData = groupData3;
        subjList_short = idx3;
        groupname = '6-8 months';
    elseif run == 4
        newData = groupData4;
        groupname = '12-15 months';
    elseif run == 5
       adpath = '';
        cd(adpath)
        load('concat_arranged_data_adult4Hz.mat');
        bbtrials = 60;
        newData = cellfun(@(x) x(:,:,[1:bbtrials]),cdataOut,'uni',false);
        subjList_short = subjList;
        groupname = 'Adults_4Hz';
        bbpath = adpath;
    end
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
    correctC = [];
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

    % start WTA
%     allWins = [];XwithWin = [];YwithWin = [];
    for n = 1:nsubj
        %windowing z-scored train and test data
        xdata = cellfun(@(x) [nanmean(x(:,LeftOT_bb),2);nanmean(x(:,RightOT_bb),2);nanmean(x(:,Posterior_bb),2)], dataInX(n,:), 'uni', false);
        ydata = cellfun(@(x) [nanmean(x(:,LeftOT_bb),2);nanmean(x(:,RightOT_bb),2);nanmean(x(:,Posterior_bb),2)], dataInY(n,:), 'uni', false);

        offdiag = [];
        for ii = 1:ncon %rows for train
            X = xdata{ii};
            if ~isempty(find(isnan(X)))
                disp(sprintf('missing values detected in training: subj %d',n))
            end
            XX = X;
            for jj = 1:ncon %column for test
                Y = ydata{jj};
                if ~isempty(find(isnan(Y)))
                    disp(sprintf('missing values detected in testing: subj %d',n))
                end
                YY = Y;
                [r,p] = corrcoef(XX,YY);
                corrmat(ii,jj,n) = r(1,2);
            end
        end
        %winner-take-all: corrmat
        for con = 1:5
            set1 = corrmat(con,:,n);
            wtaR = find(set1 == max(set1));
            if wtaR == con % if both predictions are correct
                correctC(con,n) = 1;
            else % both are incorrect
                correctC(con,n) = 0;
            end
        end
    end

    %save corrC
    filename = sprintf('WTAresults_Distinctiveness_kgsv_LOPO_OT+OCC_%s.mat',groupname)
    save(fullfile(bbpath,filename),'correctC','corrmat')
end

%% infants new groups- subplot plot
myColor = [1,0,0;
    0.929,0.694,0.125;
    0.466,0.674,0.188;
    0.1,0.1,0.1;
    0,0.447,0.741;];
grouplists = {'3-4 months' '4-6 months' '6-8 months' '12-15 months'};
figure('position',[100 100 320*4 360],'color','w')
for run = 1:4
    clear correctC
    load(sprintf('WTAresults_Distinctiveness_kgsv_LOPO_OT+OCC_%s.mat',grouplists{run}))
    tempCorrC = correctC;
    %get mean decoding across categories
    MeanCorrC = mean(tempCorrC,1);
    mcorrC = mean(tempCorrC,2);
    sd = std(tempCorrC,[],2);
    se = sd./sqrt(nsubj);
    subplot_tight(1,4,run,[0.28 0.06])
    for i = 1:5
        hold on;
        b = bar(i,mcorrC(i));
%         err = errorbar(i,mcorrC(i),se(i),'LineStyle','none','Color',myColor(i,:),'linewidth',1.2);
        b.FaceColor = myColor(i,:);
        b.FaceAlpha = 0.5;
        b.EdgeColor = myColor(i,:);
    end
    hold on;
    plot([0 6],[0.2 0.2],'linewidth',1.2,'linestyle','--','color',[0.5 0.5 0.5])
    xticks([0:1:6])
    xticklabels({'' 'Faces' 'Limbs' 'Corridors'  'Characters' 'Cars' ''})
    xtickangle(90)
    xlim([0 6])
    ylim([0 1])
    yticks([0:0.2:1])
    yticklabels([0:0.2:1]*100)
    ylabel('Decoding accuracy (%)')
    if run > 1
        ax = gca;
        ax.YAxis.Visible = 'off';
    end
    set(gca,'linewidth',1.2,'fontsize',16)
    title(grouplists{run},'fontsize',18)
end
figname = sprintf('TimeDomain_lopo_newgroups_OCC+OTrois_%s_new.tiff','allgroups');
print(gcf,fullfile(bbpath,figname),'-dtiff','-r300');
%% get mean decoding across categories - infant
mygroupcolor = [254,204,92;
    151,207,95;
    46,147,200;
    37,52,148;
    0,0,0]./255;
grouplists = {'3-4 months' '4-6 months' '6-8 months' '12-15 months'};
figure('position',[100 100 320 360],'color','w')
for run = 1:4
    hold on;
    clear correctC
    load(sprintf('WTAresults_Distinctiveness_kgsv_LOPO_OT+OCC_%s.mat',grouplists{run}))
    tempCorrC = correctC;
    MeanCorrC = mean(tempCorrC,1);
    [h,p,ci,stats] = ttest(MeanCorrC,0.2,'Tail','right');
    ttestrestuls(run).p = p;
    ttestrestuls(run).stats = stats;
    nsubj = length(MeanCorrC);
    sd = std(MeanCorrC,[],2);
    se = sd/sqrt(nsubj);
    b = bar(run,mean(MeanCorrC));
    err = errorbar(run,mean(MeanCorrC),se,'LineStyle','none','Color',mygroupcolor(run,:),'linewidth',1.2);
    b.FaceColor = mygroupcolor(run,:);
    b.FaceAlpha = 0.5;
    b.EdgeColor = mygroupcolor(run,:);
end
plot([0 5],[0.2 0.2],'linewidth',1.2,'linestyle','--','color',[0.5 0.5 0.5])
xticks([0:1:5])
xticklabels({'' '3-4 mo' '4-6 mo' '6-8 mo' '12-15 mo' ''})
xtickangle(90)
xlim([0 5])
ylim([0 1])
yticks([0:0.2:1])
yticklabels([0:0.2:1]*100)
ylabel('Decoding accuracy (%)')
set(gca,'linewidth',1.2,'fontsize',16)
figname = sprintf('TimeDomain_lopo_newgroups_OCC+OTrois_%s_avgCategories_new.tiff','allgroups');
print(gcf,fullfile(bbpath,figname),'-dtiff','-r300');

%% adult
adpath = '';
cd(adpath)
load('WTAresults_Distinctiveness_kgsv_LOPO_OT+OCC_Adults_4Hz.mat')
myColor = [1,0,0;
    0.929,0.694,0.125;
    0.466,0.674,0.188;
    0.1,0.1,0.1;
    0,0.447,0.741];
figure('position',[100 100 400 360],'color','w')
run = 5
tempCorrC = correctC;
nsubj = size(tempCorrC,2);
corrC_cat = mean(tempCorrC,1);
sd = std(corrC_cat);
se = sd/sqrt(nsubj-1);
mcorrC = mean(tempCorrC,2);
%plot
hold on;
b = bar(1,mean(corrC_cat));
err = errorbar(1,mean(corrC_cat),se,'LineStyle','none','Color','k','linewidth',1.2);
b.FaceColor = 'w';
b.EdgeColor = 'k';
b.LineWidth = 1.2;

for con = 1:5
    hold on;
    b = bar(con+1,mcorrC(con));
    b.FaceColor = myColor(con,:);
    b.FaceAlpha = 0.5;
    b.EdgeColor = myColor(con,:);
end
hold on;
plot([0 6],[0.2 0.2],'linewidth',1.2,'linestyle','--','color',[0.5 0.5 0.5])
xticks([0:1:6])
xticklabels({'' 'Mean' 'Faces' 'Limbs' 'Corridors'  'Characters' 'Cars'})
xtickangle(90)
xlim([0 7])
ylim([0 1])
yticks([0:0.2:1])
yticklabels([0:0.2:1]*100)
ylabel('Decoding accuracy (%)')
set(gca,'linewidth',1.2,'fontsize',16)
%     title('Adult','fontsize',18)

figname = 'TimeDomain_lopo_adult_OCC+OTrois_adults_pc_new.tiff';
print(gcf,fullfile(adpath,figname),'-dtiff','-r300');
%statistical test
[h,p,ci,stats] = ttest(corrC_cat,0.2,'Tail','right')

