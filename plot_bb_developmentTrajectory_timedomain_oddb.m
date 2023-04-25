%%  TIME DOMAIN - ODDBALL RESPONSE
% run this part always at the begining
clear;clc;close all;
bbpath = '';
cd(bbpath)
set(0,'DefaultAxesFontSize',14,...
    'defaultTextFontName','Calibri',...
    'defaultAxesFontName','Calibri');

% load('New3Groups_arrange_no_concat_carrier.mat');
load('New3Groups_arrange_bbages.mat')

% WX version adapted
LeftOT_bb = [57 58 59 64 65 63 68]; %removed chan68
RightOT_bb = [90 91 95 96 100 94 99];
OCC_bb = [71 76 75 70 69 74 82 83 89];
OTROI = [LeftOT_bb,RightOT_bb];
OCCROI = OCC_bb;
% roiname = 'OTROI';
% roiname = 'OCCROI';
mycolor = [254,204,92;...
    161,218,180;...
    65,182,196;...
    37,52,148]./255;
conlists = {'Faces' 'Cars' 'Corridors' 'Limbs' 'Characters'};

%% individual lmm - left + right ROI in subplots - log10Age
con = 2;
for roi = 2
    if roi ==1
        load(sprintf('TimeDomain_individual_oddball_%s_win_stats_LeftOTROI_allgroups_400_700.mat',conlists{con}));
        roiname = 'LeftOTROI'
    else
        load(sprintf('TimeDomain_individual_oddball_%s_win_stats_RightOTROI_allgroups_400_700.mat',conlists{con}));
        roiname = 'RightOTROI'
    end
    figure('position',[300 300 400 300],'color','w');

    metric = 'amp'  %amp or latency
    data = allsubj_info_amp(:,con,:);
    data = squeeze(data);
    %take absolute value all the amplitudes
    data = data;
    hold on;
    scatter(log10(group1'),data(1,1:length(group1)),40,mycolor(1,:),'filled');
    scatter(log10(group2'),data(2,1:length(group2)),40,mycolor(2,:),'^','filled');
    scatter(log10(group3'),data(3,1:length(group3)),45,mycolor(3,:),'s','filled');
    scatter(log10(group4'),data(4,1:length(group4)),40,mycolor(4,:),'d','filled');
    xtickangle(90)
    ylabel('peak amplitude (µV)')
    xlabel('log10[age(days)]')
    set(gca,'fontsize',14,'linewidth',1.2)
    title([conlists{con} '-' roiname],'fontsize',18)

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
    model(con).lme1 = lme1;
    x = 1.9:0.2:2.7; y = lme1.Coefficients.Estimate(1) + (lme1.Coefficients.Estimate(2))*x;

    %predix model
    newtbl= table(x', y', x' ,'VariableNames',{'Age','PeakResponse','Baby'});
    [Ypred,YCI] = predict(lme1,newtbl);
    [val order] = sort(x);
    plot(val,y(order), 'color', 'r','linewidth',1.5);
    patch([val, fliplr(val)],[YCI(order,1)', fliplr(YCI(order,2)')],[1 1 1],'edgecolor', [1 1 1], 'facecolor', 'r', 'FaceAlpha', 0.2)
    hold off;

    %doubld x-axis
    x1 = [1.9:0.2:2.7]; xticks(x1);
    x2 = 10.^[x1];
%     xrow1 = {'1.8' '2.0' '2.2' '2.4' '2.6' '2.8'};
%     xrow2 = {'log10(63)' '(100)' '(158)' '(251)' '(398)' '(631)'};
    xrow1 = {'1.9' '2.1' '2.3' '2.5' '2.7'};
    xrow2 = {'log10(79)' '(126)' '(200)' '(316)' '(501)'};
    labelArray = [xrow1; xrow2];
    tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));    
    ax = gca();
    ax.XTickLabel = tickLabels;
    xtickangle(0)
    %save
    filename =  sprintf('TimeDomain_oddb_%s_%s_AcrossAges_lmm_stats_%s_log10Age.mat',metric,roiname,conlists{con});
    save(fullfile(bbpath,filename),'lme1');
    figname = sprintf('TimeDomain_oddb_%s_%s_AcrossAges_withlmm_%s_log10Age.png',metric,roiname,conlists{con});
    saveas(gcf,fullfile(bbpath,figname));
end

%% lmm 3-way interacttion: category x hemi x age - log10Age
metric = 'latency'; %amp or latency
com_Amp = [];
com_data = [];
%load data
for con = 1:5
    clear allsubj_info_amp
    clear allsubj_info_latency
    %load left OT data
    filename = sprintf('TimeDomain_individual_oddball_%s_win_stats_LeftOTROI_allgroups_400_700.mat',conlists{con});
    load(filename);
    if strmatch(metric,'amp')
        %flip the sign, take only abs value
        com_data(:,:,1,con) = abs(squeeze(allsubj_info_amp(:,con,:)));
%         com_data(:,:,1,con) = squeeze(allsubj_info_amp(:,con,:));
    elseif strmatch(metric,'latency')
        com_data(:,:,1,con) = squeeze(allsubj_info_latency(:,con,:));
    end
    clear allsubj_info_amp
    clear allsubj_info_latency

    %load right OT data
    filename2 = sprintf('TimeDomain_individual_oddball_%s_win_stats_RightOTROI_allgroups_400_700.mat',conlists{con});
    load(filename2);
    if strmatch(metric,'amp')
        com_data(:,:,2,con) = abs(squeeze(allsubj_info_amp(:,con,:)));
%         com_data(:,:,2,con) = squeeze(allsubj_info_amp(:,con,:));
    elseif strmatch(metric,'latency')
        com_data(:,:,2,con) = squeeze(allsubj_info_latency(:,con,:));
    end
end
com_Amp = [];
temp_com_Amp = [];
temp_com_Amp_left = [com_data(1,1:length(group1),1,:),...
    com_data(2,1:length(group2),1,:),...
    com_data(3,1:length(group3),1,:),...
    com_data(4,1:length(group4),1,:)];
com_Amp = [com_Amp,temp_com_Amp_left];

temp_com_Amp_right = [com_data(1,1:length(group1),2,:),...
    com_data(2,1:length(group2),2,:),...
    com_data(3,1:length(group3),2,:),...
    com_data(4,1:length(group4),2,:)];
com_Amp = [com_Amp,temp_com_Amp_right];

com_Amp = squeeze(com_Amp);
new_com_Amp = reshape(com_Amp,size(com_Amp,1)*size(com_Amp,2),1);

hemiLen = length(temp_com_Amp_right);
com_hemi = [ones(hemiLen,1);ones(hemiLen,1)*2];
new_com_hemi = repmat(com_hemi,5,1);

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
temp_com_Age = [group1;group2;group3;group4];
com_Age = repmat(temp_com_Age,2,1);
new_com_Age = repmat(com_Age,5,1);
neworder = neworder';
com_neworder = repmat(neworder,2,1);
new_com_order = repmat(com_neworder,5,1);
com_category = ones(hemiLen*2,1);
new_com_category = [com_category;2*com_category;3*com_category;4*com_category;5*com_category];
nsubj = length(neworder);
tbl = table(log10(new_com_Age),new_com_Amp,new_com_category,new_com_order,new_com_hemi,'VariableNames',{'Age','PeakResp','Category','Baby','Hemi'});

%% left OT
hemilabel = 'leftOT';
hemis = tbl.Hemi;
hemi_idx = find(hemis == 2);
newtblhemi2 = tbl;
newtblhemi2(hemi_idx,:)= []
lme_hemi_left = fitlme(newtblhemi2,'PeakResp~ Category*Age +(1|Baby)')
%save
filename = sprintf('TimeDomain_oddball_%s_%s_AcrossAges_lmm_stats_interaction_log10Age.mat',metric,hemilabel);
save(fullfile(bbpath,filename),'lme_hemi_left');
%% find sig. age main effect in the right OT
hemilabel = 'rightOT'
hemis = tbl.Hemi
hemi_idx = find(hemis == 1);
newtblhemi = tbl;
newtblhemi(hemi_idx,:)= []
lme_hemi_right = fitlme(newtblhemi,'PeakResp~ Category*Age +(1|Baby)')
%save
filename = sprintf('TimeDomain_oddball_%s_%s_AcrossAges_lmm_stats_interaction_log10Age.mat',metric,hemilabel);
save(fullfile(bbpath,filename),'lme_hemi_right');

%%
lme3 = fitlme(tbl,'PeakResp~ Category*Age*Hemi +(1|Baby)');
lme4 = fitlme(tbl,'PeakResp~ Category*Age +(1|Baby)');
hemis = tbl.Hemi;
tbl_right = tbl;
tbl_left = tbl;
tbl_right(find(hemis==1),:) = [];
tbl_left(find(hemis==2),:) = [];
%right hemi only, find sig. age effect, not category
lme5 = fitlme(tbl_right,'PeakResp~ Category*Age +(1|Baby)') ;
%left hemi only, find only sig. category effect,vwith no interaction
lme6 = fitlme(tbl_left,'PeakResp~ Category*Age +(1|Baby)') ;

%save
% roiname = 'LROT'
% filename = sprintf('TimeDomain_oddball_%s_%s_AcrossAges_lmm_stats_3wayInteraction.mat',metric,roiname);
% save(fullfile(bbpath,filename),'lme3');


