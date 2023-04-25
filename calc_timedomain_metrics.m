% get time domain peak info
% ODDBALL
% not finished yet, don't use this one
% the time domain script calculates the peak metrics for each baby at
% certain time window, and plot it out in the time domain figure

clear;clc;close all;
bbpath = '';
cd(bbpath)
LeftOT_bb = [50 57 58 59 64 65]; %removed chan68 % 63
RightOT_bb = [90 91 95 96 100 101];
Posterior_bb = [66 72 71 76 75 70 69 74 82 83 84 89];
parietal_bb = [52 53 60 61 67 62 77 78 85 86 92];
central_bb = [30 37 36 42 41 7 106 31 80 55 105 104 103 87 93];
posROI = [LeftOT_bb,RightOT_bb,Posterior_bb]; %26 channels
cenROI = [parietal_bb,central_bb]; %26 channels
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = Posterior_bb;
pos_merge = [posROI,cenROI];

roiname = 'OTROI'
rois = OTROI;
% roiname = 'OCCROI'
% rois = OCCROI;

Freqlist = [0.8572 4.286];
xlen = 1166.7;
x = 0:100:1200;
xpoints = (x*420)/1000;
sampling = 490;
step = 1/420;
odds = 1/5;
xx = 1:490;
xxtt = round(xx*1000/420,0);

lists = {'Faces','Cars','Corridors','Limbs','Words'};
mycolor = [1,0,0;...
    0,0.447,0.741;...
    0.466,0.674,0.188;...
    0.929,0.694,0.125;...
    0.1,0.1,0.1];
xx = 1:490;
allsubj_info_amp = [];
allsubj_info_latency = [];
for run = 1:4    
    load('New3Groups_arrange_clean_finalbldata.mat');
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
    end
    mVdataIn = cellfun(@(x) x*10^6, newData, 'uni', false);
   
    for i = 1:5
        data = mVdataIn(:,i);
        temp = cellfun(@(x) nanmean(x(:,rois),2), data, 'uni', false);
        nsubj = size(temp,1);
        for s = 1:nsubj
            %measure the latency here
            data_temp = temp{s};
            if run < 4
                [mval,mlatency] = timeDomain_peak_metrics(data_temp,400,800,'pos');
            else
                [mval,mlatency] = timeDomain_peak_metrics(data_temp,300,700,'pos');
            end
            %plot and check if the peak is defined correctly
            scatter(mlatency,mval,'filled','d','MarkerFaceColor','b')
            allsubj_info_amp(run,i,s) = mval;
            allsubj_info_latency(run,i,s) = mlatency;
        end
    end
end

filename = sprintf('TimeDomain_individual_oddball_win_stats_%s_allgroups.mat',roiname);
save(fullfile(bbpath,filename),'allsubj_info_amp','allsubj_info_latency');
