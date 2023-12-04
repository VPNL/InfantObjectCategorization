%% measure peak latency and peak response and save to .mat file
% bb
% xy @ fudan

%% CARRIERS: INDIVIDUAL RESPONSES - subplot
clear;clc;close all;
bbpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
cd(bbpath)

LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 94 99];
Posterior_bb = [71 76 75 70 69 74 82 83 89];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; 

OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

% roiname = 'OTROI'
% rois = OTROI;
roiname = 'OCCROI'
rois = OCCROI;

Freqlist = [0.8572 4.286];
xlen = 1166.7;
x = 0:100:1200;
xpoints = (x*420)/1000;
sampling = 490;
step = 1/420;
odds = 1/5;
xx = 1:step:sampling/420;
lists = {'Faces','Cars','Corridors','Limbs','Characters'};
mycolor = [1,0,0;...
    0,0.447,0.741;...
    0.466,0.674,0.188;...
    0.929,0.694,0.125;...
    0.1,0.1,0.1];
xx = 1:98;
xxtt = round(xx*1000/420,0);

load('New3Groups_arrange_no_concat_carrier.mat');

allsubj_info_amp = [];
allsubj_info_latency = [];

for run = 1:4
    if run == 1
        newData = groupData1;
        subjList_short = idx1;
        groupname = 'Group1-3-4';
        figure('position',[10 10 1800 300*5],'color','w')
    elseif run == 2
        newData = groupData2;
        subjList_short = idx2;
        groupname = 'Group2-4-6';
        figure('position',[10 10 1800 300*3],'color','w')
    elseif run == 3
        newData = groupData3;
        subjList_short = idx3;
        groupname = 'Group3-6-8';
        figure('position',[10 10 1800 300*3],'color','w')
    elseif run == 4
        newData = groupData4;
        subjList_short = idx4;
        groupname = 'Group4-12-14';
        figure('position',[10 10 1800 300*3],'color','w')
    end

    mnewData = cellfun(@(x) squeeze(nanmean(x,3)), newData,'uni',false);
    mVdataIn = cellfun(@(x) x*10^6, mnewData, 'uni', false);

    for i = 1:size(mVdataIn,1)
        if run == 1
        subplot_tight(4,5,i)
        elseif run == 4
            subplot_tight(3,5,i)
        else
            subplot_tight(3,5,i)
        end
        hold on;
        data = mVdataIn(i);
        temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
        newtemp = reshape(temp{:},[490/5,5]);
        mtemp = nanmean(newtemp,2);
        plot(xxtt,mtemp,'color','k','linewidth',1.8,'linestyle','-');
        thresh = 12;
        yticks([-thresh:2:thresh])
        ylim([-thresh thresh])
        xlim([0 233])
        xlabel('Time (ms)')
        ylabel('Amplitude (µV)')
        plot([0,233],[0,0],'color',[0.5 0.5 0.5],'linewidth',1,'LineStyle','--');
        plot([65 65],[-thresh thresh],'linewidth',1,'LineStyle','--');
        plot([100 100],[-thresh thresh],'linewidth',1,'LineStyle','--');
        plot([160 160],[-thresh thresh],'linewidth',1,'LineStyle','--');
        set(gca,'fontsize',14,'linewidth',1.2)
        title(strrep(subjList_short{i},'_','-'),'fontweight','normal')
        
        %find peak and latency
        % win 1
        if run < 4
            [mval,mlatency] = timeDomain_peak_metrics(mtemp,60,90,'neg'); %65-100
        else
            [mval,mlatency] = timeDomain_peak_metrics(mtemp,60,90,'neg'); %65-90
        end
        scatter(mlatency,mval,'filled','d','MarkerFaceColor','b')
        allsubj_info_amp(run,i,1) = mval;
        allsubj_info_latency(run,i,1) = mlatency;
        
        % win 2
%         if run == 1
%             [mval,mlatency] = timeDomain_peak_metrics(mtemp,100,160,'pos');
%         elseif run < 4
%             [mval,mlatency] = timeDomain_peak_metrics(mtemp,90,120,'pos');
%         else
%             [mval,mlatency] = timeDomain_peak_metrics(mtemp,90,110,'pos');
%         end
        if run == 1
            [mval,mlatency] = timeDomain_peak_metrics(mtemp,90,160,'pos');
        else
            [mval,mlatency] = timeDomain_peak_metrics(mtemp,90,110,'pos');
        end
        scatter(mlatency,mval,'filled','d','MarkerFaceColor','r')
        allsubj_info_amp(run,i,2) = mval;
        allsubj_info_latency(run,i,2) = mlatency;

%         %win 3
%         [mval,mlatency] = timeDomain_peak_metrics(mtemp,160,233,'pos');
%         scatter(mlatency,mval,'filled','d','MarkerFaceColor','g')
%         allsubj_info_amp(run,i,3) = mval;
%         allsubj_info_latency(run,i,3) = mlatency;

    end
    sgtitle(groupname,'fontsize',18)
    figname = sprintf('TimeDomain_individual_carrier_%s_%s_subplots.png',groupname,roiname);
    saveas(gcf,fullfile(bbpath,figname));
end

% save window metrics
filename = sprintf('TimeDomain_individual_carrier_win_stats_%s_allgroups.mat',roiname);
save(fullfile(bbpath,filename),'allsubj_info_amp','allsubj_info_latency');


%% ODDBALL - INDIVIDUAL -subplots
clear;clc;close all;
% bbpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset/';
bbpath = '/Users/xiaoqian/Desktop/arranged_dataset';
cd(bbpath)
LeftOT_bb = [58 59 64 65 63 68]; %removed 50 & 57
RightOT_bb = [90 91 95 96 94 99]; %removed 101 & 100
Posterior_bb = [71 76 75 70 69 74 82 83 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

% roiname = 'OTROI'
% rois = OTROI;
for roi = 2
    if roi == 1
        roiname = 'LeftOTROI'
        rois = LeftOT_bb;
    else
        roiname = 'RightOTROI'
        rois = RightOT_bb;
    end
    Freqlist = [0.8572 4.286];
    xlen = 1166.7;
    x = 0:100:1200;
    xpoints = (x*420)/1000;
    sampling = 490;
    step = 1/420;
    odds = 1/5;
    xx = 1:490;
    xxtt = round(xx*1000/420,0);

    lists = {'Faces','Cars','Corridors','Limbs','Characters'};
    mycolor = [1,0,0;...
        0,0.447,0.741;...
        0.466,0.674,0.188;...
        0.929,0.694,0.125;...
        0.1,0.1,0.1];
    xx = 1:490;

    allsubj_info_amp = [];
    allsubj_info_latency = [];
    for run = 1:4
        if run == 5
            figure('position',[100 100 1800 450],'color','w')
        end
        if run < 5
            load('New3Groups_arrange_clean_finalbldata.mat');
            %         load('New3Groups_arrange_clean_noise_finalbldata.mat');
        end
        if run == 1
            newData = groupData1;
            subjList_short = idx1;
            groupname = 'Group1-3-4';
        elseif run == 2
            newData = groupData2;
            subjList_short = idx2;
            groupname = 'Group2-4-6';
        elseif run == 3
            newData = groupData3;
            subjList_short = idx3;
            groupname = 'Group3-6-8';
        elseif run == 4
            newData = groupData4;
            subjList_short = idx4;
            groupname = 'Group4-12-15';
        elseif run == 5
            adpath = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset adult/arranged_dataset/4Hz';
            dfile = 'concat_arranged_data_adult4Hz.mat';
            load(fullfile(adpath,dfile));
            %         newData = cellfun(@(x) nanmean(x(:,:,[1:46]),3),cdataOut,'uni',false);
            %         trialmode = 'bbtrials';
            newData = cellfun(@(x) nanmean(x,3),cdataOut,'uni',false);
            trialmode = 'fulltrials';
            groupname = 'adults';
        end
        mVdataIn = [];
        mVdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);
        colors = jet(17);

        for i = 1 %condition
            if run == 4
                figure('position',[100 100 1800 350*3],'color','w')
            elseif run < 4
                figure('position',[100 100 1800 350*4],'color','w')
            end
            data = mVdataIn(:,i);
            temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
            nsubj = size(temp,1);
            ff = 0;
            for n = 1:nsubj
                ff = ff + 1;
                subplot_tight(4,4,ff)
                plot(xxtt,temp{n},'color',colors(n,:),'linewidth',1.8,'linestyle','-');
                hold on;

                % measure peak latency and amplitude metrics
                mtemp = temp{n};
                time1 = 400
                time2 = 700
                if i == 1
                    [mval,mlatency] = timeDomain_peak_metrics(mtemp,time1,time2,'pos');
                else
                    [mval,mlatency] = timeDomain_peak_metrics(mtemp,time1,time2,'neg');
                end

                scatter(mlatency,mval,'filled','d','MarkerFaceColor','r')
                allsubj_info_amp(run,i,n) = mval;
                allsubj_info_latency(run,i,n) = mlatency;

                if run == 4
                    thresh = 16; %8 for signal
                    yticks([-thresh:4:thresh])
                elseif run < 4
                    thresh = 14; %8 for signal
                    yticks([-thresh:4:thresh])
                else
                    thresh = 1.5;
                    yticks([-thresh:0.5:thresh])
                end

                xlim([0 1200])
                xticks([0:400:1200])
                ylim([-thresh thresh])
                xlabel('Time (ms)')
                ylabel('Amplitude (µV)')
                h = plot([0,1200],[0,0],'color',[0.5 0.5 0.5],'linewidth',1,'LineStyle','--')
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                set(gca,'fontsize',14,'linewidth',1.2)
                title(strrep(subjList_short{n},'_','-'),'fontweight','normal','fontsize',16)
            end
            sgtitle([groupname '-' lists{i}],'fontsize',18)

            if run < 5
                %     print(gcf,fullfile(bbpath,sprintf('Timedomain-oddball-%s-30Hz-%s-overlay.tiff',roiname,groupname)),'-dtiff','-r300')
                saveas(gcf,fullfile(bbpath,sprintf('Timedomain-oddball-%s-30Hz-%s-%s-subplot-withpeak_%d-%d.png',roiname,lists{i},groupname,time1,time2)))
            else
                %     print(gcf,fullfile(adpath,sprintf('Timedomain-oddball-%s-adultGroup-%s-FDR.tiff',roiname,trialmode)),'-dtiff','-r300')
                saveas(gcf,fullfile(adpath,sprintf('Timedomain-oddball-%s-adultGroup-%s-subplot.png',roiname,trialmode)))
            end
            %         close
        end

    end

    % save window metrics
    filename = sprintf('TimeDomain_individual_oddball_%s_win_stats_%s_allgroups_%d_%d.mat',lists{i},roiname,time1,time2);
    save(fullfile(bbpath,filename),'allsubj_info_amp','allsubj_info_latency');
    disp('done')
end
%%