%% looking at developmental effects of visual resposnes 
% bb
% xy @ fudan


%%  TIME DOMAIN - visual RESPONSE
clear;clc;close all;
bbpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
cd(bbpath)

load('New3Groups_arrange_bbages.mat')
load('TimeDomain_individual_carrier_win_stats_OCCROI_allgroups.mat');
roiname = 'OCCROI';
mycolor = [245,210,15;
20,210,95;    
31,144,237;
37,52,148;]./255;

% individual lmm with log10(age)
con = 1;
metric = 'amp'; %amp or latency
winlists = {'the first time window' 'the second time window'};

for win = 1
    figure('color','w');
    if strmatch(metric,'amp')
        data = allsubj_info_amp(:,:,win);
        %take abs value
        data = abs(data);
    elseif strmatch(metric,'latency')
        data = allsubj_info_latency(:,:,win);
    end
    hold on;
    scatter(log10(group1'),data(1,1:length(group1)),10,mycolor(1,:),'filled');
    scatter(log10(group2'),data(2,1:length(group2)),10,mycolor(2,:),'^','filled');
    scatter(log10(group3'),data(3,1:length(group3)),10,mycolor(3,:),'*');
    scatter(log10(group4'),data(4,1:length(group4)),10,mycolor(4,:),'d','filled');
    xlim([1.9 2.7])
    xtickangle(90)
    if strmatch(metric,'latency')
        if win == 1
            ylabel('Peak latency (ms)')
            ylim([50 100])
        else
            ylim([80 160])
            ylabel('Peak latency (ms)')
        end
    else
        ylabel('Peak amplitude (µV)')
    end
    xlabel('Age in days (log scale)')
    set(gca,'fontsize',7,'linewidth',0.5,'fontname','arial');
    fontsize(gca,7,"points")
    fontsize(gcf,7,"points")
    fontname(gca,'Arial')
    fontname(gcf,'Arial')

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
    model(con).lme1 = lme1; x = 1.9:0.2:2.7; y =
    lme1.Coefficients.Estimate(1) + (lme1.Coefficients.Estimate(2))*x;

    %predix model 
    newtbl= table(x', y', x','VariableNames',{'Age','Amplitudes','Baby'}); 
    [Ypred,YCI] = predict(lme1,newtbl); 
    [val order] = sort(x); 
    plot(val,y(order),'color', 'r','linewidth',0.8); 
    patch([val,fliplr(val)],[YCI(order,1)', fliplr(YCI(order,2)')],[1 1 1], ...
        'edgecolor', 'none', 'facecolor', 'r', 'FaceAlpha', 0.2)
    hold off; 
    % single x-axis 
    x1 = [1.9:0.2:2.7]; 
    xticks(x1)
    xticklabels({'79' '126' '200' '316' '501'}) 
    xtickangle(0)
    set(gca,'fontsize',7,'linewidth',0.5,'fontname','arial');
    fontsize(gca,7,"points") 
    fontname(gca,'Arial') 
    fontname(gcf,'Arial')

    %save lme figpath =
    'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\PNAS-revisions-figures\';
    filename = sprintf('TimeDomain_carrier_%s_%s_AcrossAges_lmm_stats_win%d_log10Age.mat',metric,roiname,win);
    save(fullfile(figpath,filename),'lme1');
    
    %save-skinny figure 
    XimPrintsize = 21/7; YimPrintsize = 21/6;
    set(gcf,'PaperUnits','centimeters'); 
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]); 
    figname = sprintf('TimeDomain_carrier_lmm_%s_%s_win%d_log10Age_skinny.eps',metric,roiname,win);
    print(gcf,fullfile(figpath,figname),'-depsc','-r600');

    %save-chunky version figure 
    XimPrintsize = 21/6; YimPrintsize = 21/8;
    set(gcf,'PaperUnits','centimeters'); 
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]); 
    figname = sprintf('TimeDomain_carrier_lmm_%s_%s_win%d_log10Age_chunky.eps',metric,roiname,win);
    print(gcf,fullfile(figpath,figname),'-depsc','-r600');

end
close all

%% lmm with interactions - log10
metric = 'amp'; %amp or latency
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
%     model(con).lme1 = lme1;
%     x = 60:500; y = lme1.Coefficients.Estimate(1) + (lme1.Coefficients.Estimate(2))*x;

%     %predix model
%     newtbl= table(x', y', x' ,'VariableNames',{'Age','Amplitudes','Baby'});
%     [Ypred,YCI] = predict(lme1,newtbl);
%     [val order] = sort(x);
%     plot(val,y(order), 'color', 'r','linewidth',1.5);
%     patch([val, fliplr(val)],[YCI(order,1)', fliplr(YCI(order,2)')],[1 1 1],'edgecolor', [1 1 1], 'facecolor', 'r', 'FaceAlpha', 0.2)
%     hold off;

%save
filename = sprintf('TimeDomain_carrier_%s_%s_AcrossAges_lmm_stats_interaction_log10Age.mat',metric,roiname);
save(fullfile(bbpath,filename),'lme2');

figname = sprintf('TimeDomain_carrier_%s_%s_AcrossAges_withlmm.png',metric,roiname);
saveas(gcf,fullfile(bbpath,figname));





