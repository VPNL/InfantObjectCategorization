% plot distinctive matrix
clc;clear;close all;
% bbDir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/eeg/data/prelimresults/2021 dataset bb/arranged_dataset';
bbDir = '/Users/xiaoqian/Desktop/adult_4hz';
cd(bbDir)
allData = [];
grouplists = {'Adults_4Hz'};
data1 = sprintf('WTAresults_Distinctiveness_kgsv_OT+OCC_%s.mat',grouplists{1}); 
load(data1);
len1 = size(DistValue,2);
allData = DistValue;
roiname = 'OT+OCC';
conlists = {'Faces','Cars','Corridors','Limbs','Characters'};

%% in subplot distinctiveness plot
close all
figure('position',[100 100 650 180],'color','w')
imagesc(allData,[-1 1])
cmap=mrvColorMaps('coolhot'); 
colormap(cmap);
c = colorbar;
xticks([1:len1])
yticklabels(conlists)
set(gca,'fontsize',14, 'linewidth',1.2)
xlabel('Participants','fontsize',16)
title('Distinctiveness values','fontsize',18)

figname = [roiname,'-DistinctivenessMatrix-allcategories-kgsv.png'];
saveas(gcf,fullfile(bbDir,figname))

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

figure('position',[100 100 2000 400],'color','w')
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
mycolor = [254,204,92;...
    161,218,180;...
    44,127,184;...
    37,52,148]./255;
data = allData;
figure('position',[100 100 2000 400],'color','w')

conlists = {'Faces','Cars','Corridors','Limbs','Characters'};
myColor = [1,0,0;0,0.447,0.741;0.466,0.674,0.188;0.929,0.694,0.125;0.1,0.1,0.1];
for nc = 1:5
%     figure('position',[100 100 500 360],'color','w')
    subplot_tight(1,5,nc,[0.2,0.032])
    hold on;
%     scatter([1:57],mwinData(nc,:),20,myColor(nc,:),'filled')
    scatter(log10(group1'),data(nc,1:len1),40,mycolor(1,:),'filled');
    scatter(log10(group2'),data(nc,len1+1:len1+len2),40,mycolor(2,:),'^','filled');
    scatter(log10(group3'),data(nc,len1+len2+1:len1+len2+len3),45,mycolor(3,:),'s','filled');
    scatter(log10(group4'),data(nc,len1+len2+len3+1:len1+len2+len3+len4),40,mycolor(4,:),'d','filled');
    if nc == 1
        ylabel('Mean distinctiveness')
    end
    xlabel('log10[Age (days)]')
    set(gca,'fontsize',14,'linewidth',1.2)
    ylim([-1 1.5])

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
    patch([val, fliplr(val)],[YCI(order,1)', fliplr(YCI(order,2)')],[1 1 1],'edgecolor', [1 1 1], 'facecolor', 'r', 'FaceAlpha', 0.2)
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
end

figname = ['scatter-DistinctivenessMatrix-kgsv-OT+OCC-log10.png'];
saveas(gcf,fullfile(bbDir,figname))
%save the individual lmm
filename = 'Stats-Distinctiveness_AcrossAges_withlmm_indivCategory_OT+OCC-log10age.mat';
save(fullfile(bbDir,filename),'model');

%fdr correction
pvals = [7.9531e-06 0.00089091 0.0375 0.00017818 0.072441];
alpha_level = 0.05;
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,alpha_level,'pdep','yes');
h
filename = 'Stats-fdr-Distinctiveness_AcrossAges_withlmm_indivCategory_OT+OCC-log10age.mat';
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

filename = 'meanDistinctiveness_AcrossAges_withlmm_OT+OCC_log10age.mat'
save(fullfile(bbDir,filename),'tbl2','lme2');


