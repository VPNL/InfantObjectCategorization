%% ODDBALLS: left+right each category in subplot + add rightOT LMM
% bb
% xy @ fudan
%% part 1: plot time domain

clear;clc;close all;
bbpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
cd(bbpath)
set(0,'DefaultAxesFontSize',8,...
    'defaultTextFontName','Arial',...
    'defaultAxesFontName','Arial');
LeftOT_bb = [58 59 64 65 63 68 57]; %removed chan68 % 63 % 50 as top 57
RightOT_bb = [90 91 95 96 99 94 100]; % 101 as top 100
Posterior_bb = [71 76 75 70 69 74 82 83 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

Freqlist = [0.8572 4.286];
xlen = 1166.7;
x1 = 0:100:1200;
xpoints = (x1*420)/1000;
sampling = 490;
step = 1/420;
odds = 1/5;
xx = 1:step:sampling/420;
conlists = {'Faces','Cars','Corridors','Limbs','Characters'};
mycolor = [245,210,15;
    20,210,95;
    31,144,237;
    37,52,148;]./255;
xx = 1:490;

load('New3Groups_arrange_clean_finalbldata.mat');
%load i = 1:5 seperately
for i = 4 % condition: faces, limbs, corridors, characters, cars
    for run = 1:4
        if run == 1
            newData = groupData1;
            subjList_short = idx1;
            groupname = 'Group1-3-4';
            groupname = '3-4 months';
        elseif run == 2
            newData = groupData2;
            subjList_short = idx2;
            groupname = 'Group2-4-6';
            groupname = '4-6 months';
        elseif run == 3
            newData = groupData3;
            subjList_short = idx3;
            groupname = 'Group3-6-8';
            groupname = '6-8 months';
        elseif run == 4
            newData = groupData4;
            subjList_short = idx4;
            groupname = 'Group4-12-14';
            groupname = '12-15 months';
        end
        if i == 1
            thresh = 8;
        else
            thresh = 6; %8 for signal
        end
        mVdataIn = [];
        mVdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);
        %save individual plot
        figure('color','w');hold on;
        for roi = 1:2
            if roi == 1
                roiname = 'LeftOTROI'
                rois = LeftOT_bb;
            else
                roiname = 'RightOTROI'
                rois = RightOT_bb;
            end
            data = mVdataIn(:,i);
            temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
            %             data_temp = cat(2,temp{:});
            data_temp = combineCells(temp,1,0);
            data_temp = squeeze(data_temp{:});
            inData = data_temp;
            %use svdnl code for ttest and cluster-based correction
            [realH, realP, realT, corrH, critVal, supraTh, randDist]= ttest_permute_sstats(inData,10000,'mass');

            meand = nanmean(data_temp,2);
            nsubj = size(data_temp,2);
            se = std(data_temp,[],2)/sqrt(nsubj-1);
            if roi == 1
                plot(meand,'color',[0.5 0.5 0.5],'linewidth',0.8,'linestyle','-');
                err = patch([xx fliplr(xx)],...
                    [(meand-se)',fliplr([meand+se]')],...
                    [0.5 0.5 0.5],'FaceAlpha',0.2, 'EdgeColor','none');
                err.Annotation.LegendInformation.IconDisplayStyle = 'off';
                %sig. time points
                b1 = scatter(xx(find(corrH)),1-thresh*ones(length(find(corrH)),1),2,[0.5 0.5 0.5],'*');
                b1.Annotation.LegendInformation.IconDisplayStyle = 'off';

            else
                plot(meand,'color',mycolor(run,:),'linewidth',0.8,'linestyle','-');
                err = patch([xx fliplr(xx)],...
                    [(meand-se)',fliplr([meand+se]')],...
                    mycolor(run,:),'FaceAlpha',0.2, 'EdgeColor','none');
                err.Annotation.LegendInformation.IconDisplayStyle = 'off';
                %sig. time points
                b1 = scatter(xx(find(corrH)),0.5-thresh*ones(length(find(corrH)),1),2,mycolor(run,:),'*');
                b1.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
        end

        xticks(xpoints);
        xticklabels(x1)
        xtickangle(90)

        if run == 1
            yticks([-thresh:2:thresh])
        else
            h = gca;
            h.YAxis.Visible = 'off';
        end
        ylim([-thresh thresh])
        xlim([0 490])
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        h = plot([0,490],[0,0],'color',[0.5 0.5 0.5],'linewidth',0.5,'LineStyle','--')
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        set(gca,'linewidth',0.5)
        fontsize(gca,7,"points")
        fontsize(gcf,7,"points")
        fontname(gca,'Arial')
        fontname(gcf,'Arial')
        % save
        XimPrintsize = 21/7;
        YimPrintsize = 21/6;
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
        figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
        figname = sprintf('TimeDomain_withlmm_oddb_%s_%s_%s_log10Age.eps','LeftRightOT',conlists{i},groupname);
        print(gcf,fullfile(figpath,figname),'-depsc','-r600');

    end
    close all
end

%% part 2: plot LMM for each condition
% last column: add lmm
load('New3Groups_arrange_bbages.mat')
for i = 4
    con = i;
    for roi = 2
        if roi ==1
            load(sprintf('TimeDomain_individual_oddball_%s_win_stats_LeftOTROI_allgroups_400_700.mat',conlists{con}));
            roiname = 'LeftOTROI'
        else
            load(sprintf('TimeDomain_individual_oddball_%s_win_stats_RightOTROI_allgroups_400_700.mat',conlists{con}));
            roiname = 'RightOTROI'
        end
        figure('color','w')

        metric = 'amp'  %amp or latency
        data = allsubj_info_amp(:,con,:);
        data = squeeze(data);
        %take absolute value all the amplitudes
        data = data;
        hold on;
        scatter(log10(group1'),data(1,1:length(group1)),10,mycolor(1,:),'filled');
        scatter(log10(group2'),data(2,1:length(group2)),10,mycolor(2,:),'^','filled');
        scatter(log10(group3'),data(3,1:length(group3)),10,mycolor(3,:),'*');
        scatter(log10(group4'),data(4,1:length(group4)),10,mycolor(4,:),'d','filled');
        xtickangle(90)
        ylabel('Peak amplitude (µV)')
        xlabel('Age in days (log10)')
        set(gca,'linewidth',0.5)

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
        com_Amp = [data(1,1:length(group1)),...
            data(2,1:length(group2)),...
            data(3,1:length(group3)),...
            data(4,1:length(group4))];
        max(com_Amp)
        min(com_Amp)
        if con == 1
            ylim([-5 20])
        else
            ylim([-15 5])
        end

        tbl = table(log10(com_Age),com_Amp',neworder','VariableNames',{'Age','PeakResponse','Baby'});
        lme1 = fitlme(tbl,'PeakResponse~ Age +(1|Baby)')
%         model(con).lme1 = lme1;
%         x = 1.9:0.2:2.7; y = lme1.Coefficients.Estimate(1) + (lme1.Coefficients.Estimate(2))*x;
% 
%         %predix model
%         newtbl= table(x', y', x' ,'VariableNames',{'Age','PeakResponse','Baby'});
%         [Ypred,YCI] = predict(lme1,newtbl);
%         [val order] = sort(x);
%         p1 = plot(val,y(order), 'color', 'r','linewidth',0.8);
%         p2 = patch([val, fliplr(val)],[YCI(order,1)', fliplr(YCI(order,2)')],[1 1 1],'edgecolor', 'none', 'facecolor', 'r', 'FaceAlpha', 0.2)
%         p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         hold off;
%         x1 = [1.9:0.2:2.7]; xticks(x1);
%         xticklabels({'79' '126' '200' '316' '501'})
%         xtickangle(90)
%         fontsize(gca,7,"points")
%         fontsize(gcf,7,"points")
%         fontname(gca,'Arial')
%         fontname(gcf,'Arial')
    end
    % save
    XimPrintsize = 21/6.8;
    YimPrintsize = 21/6;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
    figpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
    figname = sprintf('TimeDomain_lmm_oddb_%s_%s_%s_log10Age.eps',metric,roiname,conlists{con});
    print(gcf,fullfile(figpath,figname),'-depsc','-r600');

end
% close all