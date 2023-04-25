% plot distinctive matrix
clc;clear;close all;
% bbDir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset';
% bbDir = '/Users/xiaoqian/Desktop/arranged_dataset';
% cd(bbDir)
bbpath = 'D:\Stanford_infant_EEG\infant_EEG_mac_allfiles\arranged_dataset\';
cd(bbpath)
bbDir = bbpath;
set(0,'DefaultAxesFontSize',14,...
    'defaultTextFontName','Calibri',...
    'defaultAxesFontName','Calibri');

allData = [];
grouplists = {'3-4 months' '4-6 months' '6-8 months' '12-15 months'};
data1 = sprintf('WTAresults_Distinctiveness_kgsv_noconcat_noROIavg_OT+OCC_%s.mat',grouplists{1}); 
load(data1);
len1 = size(DistValue,2);
allData = DistValue;
data2 = sprintf('WTAresults_Distinctiveness_kgsv_noconcat_noROIavg_OT+OCC_%s.mat',grouplists{2}); 
load(data2);
len2 = size(DistValue,2);
allData = cat(2,allData,DistValue);
data3 = sprintf('WTAresults_Distinctiveness_kgsv_noconcat_noROIavg_OT+OCC_%s.mat',grouplists{3}); 
load(data3);
len3 = size(DistValue,2);
allData = cat(2,allData,DistValue);
data4 = sprintf('WTAresults_Distinctiveness_kgsv_noconcat_noROIavg_OT+OCC_%s.mat',grouplists{4}); 
load(data4);
len4 = size(DistValue,2);
allData = cat(2,allData,DistValue);

roiname = 'OT+OCC';
conlists = {'Faces' 'Limbs' 'Corridors'  'Characters' 'Cars'};

%% in subplot distinctiveness plot
close all
load('New3Groups_arrange_bbages.mat');
figure('position',[100 100 1300 600],'color','w')
for con = 1:5
    subplot_tight(5,1,con,[0.13,0.06])
    imagesc(allData(con,:),[-0.5 0.5])
    cmap=mrvColorMaps('coolhot');
    colormap(cmap);
%     c = colorbar;
    yticks([0 1])
    yticklabels({''})
    
    xticks([1:1:len1+len2+len3+len4])
    xticklabels([group1;group2;group3;group4])
    xtickangle(90)
    
    set(gca,'TickLength',[0 0])
    set(gca,'fontsize',14, 'linewidth',1.2)
    title(conlists{con},'fontsize',18,'FontWeight','bold')
end
xlabel('Age (days)','fontsize',20,'FontWeight','bold')
% sgtitle('Distinctiveness values','fontsize',18)

figname = [roiname,'-DistinctivenessMatrix-allcategories-kgsv-noROIavg-pc.tiff'];
% saveas(gcf,fullfile(bbDir,figname))
print(gcf,fullfile(bbDir,figname),'-dtiff','-r300')

%% 416-499 ms time window - need to be edited
load('New3Groups_arrange_bbages.mat')
% rmidx = strmatch('bb39_3',idx1);
% idx1(rmidx) = [];
% group1(rmidx) = [];
len1 = length(idx1);
len2 = length(idx2);
len3 = length(idx3);
len4 = length(idx4);
mycolor = [254,204,92;...
    161,218,180;...
    44,127,184;...
    37,52,148]./255;

figure('position',[100 100 1300 400],'color','w')
conlists = {'Faces','Cars','Corridors','Limbs','Characters'};
data = allData;
for nc = 1:5
    subplot_tight(1,5,nc,[0.16,0.04])
    hold on;
%     scatter([1:57],mwinData(nc,:),20,myColor(nc,:),'filled')
    scatter(group1',data(nc,1:len1),40,mycolor(1,:),'filled');
    scatter(group2',data(nc,len1+1:len1+len2),40,mycolor(2,:),'^','filled');
    scatter(group3',data(nc,len1+len2+1:len1+len2+len3),45,mycolor(3,:),'s','filled');
    scatter(group4',data(nc,len1+len2+len3+1:len1+len2+len3+len4),40,mycolor(4,:),'d','filled');

    ylabel('Mean distinctiveness')
    xlabel('Age (days)')
    set(gca,'fontsize',12,'linewidth',1.2)
    title(conlists{nc},'fontsize',18)
    ylim([-1 1.5])
end

figname = ['scatter-DistinctivenessMatrix-',num2str(timewin1),'-',num2str(timewin2),'-ms',conlists{nc},'.png'];
saveas(gcf,fullfile(bbDir,figname))


%run lmm on with age & category as fixed factor
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
com_Age = repmat(temp_com_Age,5,1);
neworder = neworder';
com_neworder = repmat(neworder,5,1);
nsubj = length(neworder);
com_category = [ones(nsubj,1);2*ones(nsubj,1);3*ones(nsubj,1);4*ones(nsubj,1);5*ones(nsubj,1)];
com_Amp = [data(1,:),data(2,:),data(3,:),data(4,:),data(5,:)];

tbl = table(com_Age,com_Amp',com_neworder,com_category,'VariableNames',{'Age','Distinctiveness','Baby','Category'});

lme1 = fitlme(tbl,'Distinctiveness ~ Age*Category +(1|Baby)')

filename = sprintf('meanDistinctiveness_%s-%sms_AcrossAges_withlmm.mat',num2str(timewin1),num2str(timewin2))
save(fullfile(bbDir,filename),'tbl','lme1');


%% log10(age)
pvals = [];
load('New3Groups_arrange_bbages.mat')
% rmidx = strmatch('bb39_3',idx1);
% idx1(rmidx) = [];
% group1(rmidx) = [];
len1 = length(idx1);
len2 = length(idx2);
len3 = length(idx3);
len4 = length(idx4);
% mycolor = [254,204,92;...
%     161,218,180;...
%     46,147,200;...
%     37,52,148]./255;
mycolor = [254,204,92;
    151,207,95;
    46,147,200;
    37,52,148;
    0,0,0]./255;
data = allData;
figure('position',[100 100 1300 360],'color','w')
conlists = {'Faces' 'Limbs' 'Corridors'  'Characters' 'Cars'};
myColor = [1,0,0;0,0.447,0.741;0.466,0.674,0.188;0.929,0.694,0.125;0.1,0.1,0.1];
for nc = 1:5
%     figure('position',[100 100 500 360],'color','w')
    if nc == 1
        subplot_tight(1,5,nc,[0.25,0.05])
    else
        subplot_tight(1,5,nc,[0.25,0.04])
    end
    hold on;
    plot([1.9 2.7],[0 0],'color','k','linewidth',1);
%     scatter([1:57],mwinData(nc,:),20,myColor(nc,:),'filled')
    scatter(log10(group1'),data(nc,1:len1),30,mycolor(1,:),'filled');
    scatter(log10(group2'),data(nc,len1+1:len1+len2),40,mycolor(2,:),'^','filled');
    scatter(log10(group3'),data(nc,len1+len2+1:len1+len2+len3),50,mycolor(3,:),'*');
    scatter(log10(group4'),data(nc,len1+len2+len3+1:len1+len2+len3+len4),40,mycolor(4,:),'d','filled');
    if nc == 1
        ylabel('Distinctiveness')
    end
    set(gca,'fontsize',16,'linewidth',1.2)
    ylim([-1 1])

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
    com_Age = repmat(temp_com_Age,1,1);
    neworder = neworder';
    com_Amp = data(nc,:);
    tbl = table(log10(com_Age),com_Amp',neworder,'VariableNames',{'Age','Distinctiveness','Baby'});
    lme1 = fitlme(tbl,'Distinctiveness~ Age +(1|Baby)')
    model(nc).lme1 = lme1;
    x = 1.9:0.2:2.7; y = lme1.Coefficients.Estimate(1) + (lme1.Coefficients.Estimate(2))*x;
    %predix model
    newtbl= table(x', y', x' ,'VariableNames',{'Age','Distinctiveness','Baby'});
    [Ypred,YCI] = predict(lme1,newtbl);
    [val order] = sort(x);
    plot(val,y(order), 'color', 'r','linewidth',1.5);
    patch([val, fliplr(val)],[YCI(order,1)', fliplr(YCI(order,2)')],[1 1 1],'edgecolor', 'none', 'facecolor', 'r', 'FaceAlpha', 0.2)
    hold off;
    pvals(nc) = lme1.Coefficients.pValue(2);
    if lme1.Coefficients.pValue(2) < 0.001
        title([conlists{nc} ' ***'],'fontsize',20)
    elseif lme1.Coefficients.pValue(2) < 0.01
        title([conlists{nc} ' **'],'fontsize',20)
    elseif lme1.Coefficients.pValue(2) < 0.05
        title([conlists{nc} ' *'],'fontsize',20)
    else
        title(conlists{nc},'fontsize',20)
    end
    
    if nc == 1
%         x1 = [1.9:0.2:2.7]; xticks(x1);
%         x2 = 10.^[x1];
%         %     xrow1 = {'1.8' '2.0' '2.2' '2.4' '2.6' '2.8'};
%         %     xrow2 = {'log10(63)' '(100)' '(158)' '(251)' '(398)' '(631)'};
%         xrow1 = {'1.9' '2.1' '2.3' '2.5' '2.7'};
%         xrow2 = {'log10(79)' '(126)' '(200)' '(316)' '(501)'};
%         labelArray = [xrow1; xrow2];
%         tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
%         ax = gca();
%         ax.XTickLabel = tickLabels;
%         xtickangle(0)
%         xlabel('log10[Age (days)]')
        x1 = [1.9:0.2:2.7]; xticks(x1);
        xticklabels({'79' '126' '200' '316' '501'})
        xlabel('Age in days (log scale)')
    end
    
    if nc > 1
        h = gca;
        h.YAxis.Visible = 'off';
        h.XTick = [];  
    end
    ax = gca();
    ax.XAxis.FontSize = 13;
end

figname = ['scatter-DistinctivenessMatrix-kgsv-noROIavg-OT+OCC-log10-pc-new.tiff'];
% saveas(gcf,fullfile(bbDir,figname))
print(gcf,fullfile(bbDir,figname),'-dtiff','-r300')
%save the individual lmm
filename = 'Stats-Distinctiveness_AcrossAges_withlmm_indivCategory_noROIavg_OT+OCC-log10age-pc.mat';
save(fullfile(bbDir,filename),'model');
%%
%fdr correction
alpha_level = 0.05;
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,alpha_level,'pdep','yes');
multcomp_fdr_bh(pvals, 'alpha', 0.05)
filename = 'Stats-fdr-Distinctiveness_AcrossAges_withlmm_indivCategory_noROIavg_OT+OCC-log10age-pc.mat';
save(fullfile(bbDir,filename),'h','alpha_level','pvals');

%% run lmm on with age & category as fixed factor
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
com_Age = repmat(temp_com_Age,5,1);
neworder = neworder';
com_neworder = repmat(neworder,5,1);
nsubj = length(neworder);
com_category = [ones(nsubj,1);2*ones(nsubj,1);3*ones(nsubj,1);4*ones(nsubj,1);5*ones(nsubj,1)];
com_Amp = [data(1,:),data(2,:),data(3,:),data(4,:),data(5,:)];

tbl2 = table(log10(com_Age),com_Amp',com_neworder,com_category,'VariableNames',{'Age','Distinctiveness','Baby','Category'});

lme2 = fitlme(tbl2,'Distinctiveness ~ Age*Category +(1|Baby)')

filename = sprintf('meanDistinctiveness_AcrossAges_withlmm_noROIavg_OT+OCC_log10age.mat');
save(fullfile(bbDir,filename),'tbl2','lme2');


