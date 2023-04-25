%%  TIME DOMAIN - CARRIER RESPONSE
clear;clc;close all;
bbpath = '';
cd(bbpath)
set(0,'DefaultAxesFontSize',14,...
    'defaultTextFontName','Calibri',...
    'defaultAxesFontName','Calibri');
load('New3Groups_arrange_bbages.mat')
load('TimeDomain_individual_carrier_win_stats_OCCROI_allgroups.mat');
roiname = 'OCCROI';
mycolor = [254,204,92;
    151,207,95;
    46,147,200;
    37,52,148;
    0,0,0]./255;

 
%% individual lmm with log10(age)
con = 1;
metric = 'latency'; %amp or latency
%window 1
winlists = {'the first time window' 'the second time window'};
figure('position',[100 100 400*2 400],'color','w');
for win = 1:2
    subplot_tight(1,2,win,[0.24 0.1])
    if strmatch(metric,'amp')
        data = allsubj_info_amp(:,:,win);
        %take abs value
        data = abs(data);
    elseif strmatch(metric,'latency')
        data = allsubj_info_latency(:,:,win);
    end
    hold on;
    scatter(log10(group1'),data(1,1:length(group1)),30,mycolor(1,:),'filled');
    scatter(log10(group2'),data(2,1:length(group2)),40,mycolor(2,:),'^','filled');
    scatter(log10(group3'),data(3,1:length(group3)),50,mycolor(3,:),'*');
    scatter(log10(group4'),data(4,1:length(group4)),40,mycolor(4,:),'d','filled');
    xlim([1.9 2.7])
    xtickangle(90)
    if strmatch(metric,'latency')
        if win == 1
            ylabel('Peak latency (ms)')
            ylim([50 100])
        else
            ylim([80 160])
        end
    else
        ylabel('Peak amplitude (µV)')
    end
    xlabel('Age in days (log scale)')
    set(gca,'fontsize',14,'linewidth',1.2)
    title(winlists{win},'fontsize',18)

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
    %take abs value of the amplitude
    tbl1 = table(log10(com_Age),com_Amp',neworder','VariableNames',{'Age','PLatency','Baby'});
    lme1 = fitlme(tbl1,'PLatency~ Age +(1|Baby)')
    model(con).lme1 = lme1;
    x = 1.9:0.2:2.7; y = lme1.Coefficients.Estimate(1) + (lme1.Coefficients.Estimate(2))*x;

    %predix model
    newtbl= table(x', y', x' ,'VariableNames',{'Age','Amplitudes','Baby'});
    [Ypred,YCI] = predict(lme1,newtbl);
    [val order] = sort(x);
    plot(val,y(order), 'color', 'r','linewidth',1.5);
    patch([val, fliplr(val)],[YCI(order,1)', fliplr(YCI(order,2)')],[1 1 1],'edgecolor', 'none', 'facecolor', 'r', 'FaceAlpha', 0.2)
    hold off;
    
    x1 = [1.9:0.2:2.7];
    xticks(x1)
    xticklabels({'79' '126' '200' '316' '501'})
    xtickangle(0)
    %save
    filename =  sprintf('TimeDomain_carrier_%s_%s_AcrossAges_lmm_stats_win%d_log10Age_pc.mat',metric,roiname,win);
    save(fullfile(bbpath,filename),'lme1');
end
figname = sprintf('TimeDomain_carrier_%s_%s_AcrossAges_withlmm_log10Age_pc_new.tiff',metric,roiname);
% saveas(gcf,fullfile(bbpath,figname));
print(gcf,fullfile(bbpath,figname),'-dtiff','-r300')

%% lmm with interactions - log10
metric = 'latency'; %amp or latency
com_Amp = [];
nwin = 2;
for win = 1:nwin
    if strmatch(metric,'amp')
        data = allsubj_info_amp(:,:,win);
    elseif strmatch(metric,'latency')
        data = allsubj_info_latency(:,:,win);
    end
    temp_com_Amp = [];
    temp_com_Amp = [data(1,1:length(group1)),...
        data(2,1:length(group2)),...
        data(3,1:length(group3)),...
        data(4,1:length(group4))];
    com_Amp = [com_Amp,temp_com_Amp];
end

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
com_Age = repmat(temp_com_Age,nwin,1);
neworder = neworder';
com_neworder = repmat(neworder,nwin,1);
nsubj = length(neworder);
com_win = [repmat(1,nsubj,1);repmat(2,nsubj,1)];

tbl = table(log10(com_Age),com_Amp',com_neworder,com_win,'VariableNames',{'Age','PLatency','Baby','Wins'});
lme2 = fitlme(tbl,'PLatency~ Age*Wins +(1|Baby)')
