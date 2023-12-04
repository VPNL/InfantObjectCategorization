%% fig1- image low-level properties analysis 
% following Stigliani et al (2015)
% calculate image luminance
% calculate image contrast
% calculate image similarity
% calculate the mean image
% calculate image spatial frequency

% xiaoqian yan 2023 july at fudan

%% image luminance
clc;close all;clear
imgpath = 'D:\Stanford_infant_EEG\fLoc-master\stimuli\';
cd(imgpath);

categories = {'adult' 'limb' 'corridor' 'word' 'car'};
ncat = length(categories);
AllImgs = [];
AllImgs_stats = [];
lums = []; sds = [];
for n = 1:ncat
    newpath = fullfile(imgpath,categories{n});
    cd(newpath)
    files = dir('*.jpg');
    fLen = size(files,1);
    for f = 1:fLen
        temp = (imread(files(f).name))./255;
        lum_temp = mean2(temp);
        AllImgs(f) = lum_temp;
        AllImgs_stats(f,n) = lum_temp;
    end
    lums(n) = mean(AllImgs);
    sds(n) = std(AllImgs);
end
%% stats
[p,t,st] = anova1(AllImgs_stats(:,[1:5]))

%% box plot
mycolor = [1,0,0;...
    0.929,0.694,0.125;...
    0.466,0.674,0.188;...
    0.4,0.4,0.4;...
    0,0.447,0.741];

mycolor_inv = flipud(mycolor);
d2plot = (AllImgs_stats);

figure('color','w')
bb = boxplot(d2plot,'symbol','k+','OutlierSize',4,'notch','on','Colors',mycolor); % datapoits in x rows and group in columns

h = findobj(gca,'Tag', 'Box');
for i = 1:length(h)
    boxEdges = findobj(h(i), 'type', 'line', 'Tag', 'Box');
    set(boxEdges, 'Color', mycolor_inv(i,:));
    set(boxEdges, 'LineWidth', 1);
end

outliers = findobj(bb, 'Tag', 'Outliers');
% Set the color of the outliers
for i = 1:length(outliers)
    set(outliers(i), 'MarkerEdgeColor', mycolor(i,:));
end

set(bb,'LineWidth',0.6) %0.4
ylim([0.6 0.7])
yticks([0.6:0.1:0.7])
xlim([0 6])
% ylabel('Normalized luminance')
box off

% % Apply colors
% h = findobj(gca,'Tag','Box');
% for jj=1:length(h)
%    patch(get(h(jj),'XData') ,get(h(jj),'YData'),mycolor_inv(jj,:));
% end

set(gca,'LineWidth',0.5,'xtick',[])
fontsize(gca,7,'points')
%save
XimPrintsize = 2.6;
YimPrintsize = 2.6;
set(gcf,'PaperUnits','centimeters');
set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
figpath = 'D:\Stanford_infant_EEG\fLoc-analysis\';
figname = 'luminance-boxplot-v2';
% saveas(gcf, fullfile(figpath,figname))
print(gcf,fullfile(figpath,figname),'-depsc','-r600');

%% contrast
clc;close all;clear

imgpath = 'D:\Stanford_infant_EEG\fLoc-master\stimuli\';
cd(imgpath);

categories = {'adult' 'limb' 'corridor' 'word' 'car'};
ncat = length(categories);
allcontr_stats = [];

for n = 1:ncat
    newpath = fullfile(imgpath,categories{n});
    cd(newpath)
    files = dir('*.jpg');
    fLen = size(files,1);
    for f = 1:fLen
        temp = imread(files(f).name)./255;
%         temp_max = max(temp,[],'all');
%         temp_min = min(temp,[],'all');
%         contr(f) = (temp_max - temp_min)/(temp_max + temp_min);
        contr(f) = std2(single(temp));
        allcontr_stats(f,n) = std2(single(temp));
    end
    allCon(n) = mean(contr);
    allCon_sd(n) = std(single(contr)); %all images have the same max and min so that we cannot calculate the SD.
end
%% stats
[p,t,st] = anova1(allcontr_stats(:,[2:5]))

%% box plot 
mycolor = [1,0,0;...
    0.929,0.694,0.125;...
    0.466,0.674,0.188;...
    0.4,0.4,0.4;...
    0,0.447,0.741];
mycolor_inv = flipud(mycolor);
d2plot = (allcontr_stats);
figure('color','w')
% bb = boxplot(d2plot,'symbol','k+','OutlierSize',4,'notch','on'); % datapoits in x rows and group in columns

bb = boxplot(d2plot,'symbol','k+','OutlierSize',4,'notch','on','Colors',mycolor); % datapoits in x rows and group in columns

h = findobj(gca,'Tag', 'Box');
for i = 1:length(h)
    boxEdges = findobj(h(i), 'type', 'line', 'Tag', 'Box');
    set(boxEdges, 'Color', mycolor_inv(i,:));
    set(boxEdges, 'LineWidth', 1);
end

outliers = findobj(bb, 'Tag', 'Outliers');
% Set the color of the outliers
for i = 1:length(outliers)
    set(outliers(i), 'MarkerEdgeColor', mycolor(i,:));
end

set(bb,'LineWidth',0.6) %0.4
ylim([0.4 0.5])
yticks([0.4 0.5])
xlim([0 6])
box off

% % Apply colors
% h = findobj(gca,'Tag','Box');
% for jj=1:length(h)
%    patch(get(h(jj),'XData') ,get(h(jj),'YData'),mycolor_inv(jj,:));
% end

set(gca,'LineWidth',0.5,'xtick',[])
fontsize(gca,7,'points')
%save
XimPrintsize = 2.6;
YimPrintsize = 2.6;
set(gcf,'PaperUnits','centimeters');
set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);

figpath = 'D:\Stanford_infant_EEG\fLoc-analysis\';
figname = 'contrast-boxplot-v2';
% saveas(gcf, fullfile(figpath,figname))
print(gcf,fullfile(figpath,figname),'-depsc','-r600');



%% similarity
clc;close all;clear

imgpath = 'D:\Stanford_infant_EEG\fLoc-master\stimuli\';
cd(imgpath);
categories = {'adult' 'limb' 'corridor' 'word' 'car'};
ncat = length(categories);
for n = 1:ncat
    newpath = fullfile(imgpath,categories{n});
    cd(newpath)
    files = dir('*.jpg');
    fLen = size(files,1);
    allpos = nchoosek([1:fLen],2);
    for i = 1:size(allpos,1)
        name_img1 = [categories{n} '-' num2str(allpos(i,1)) '.jpg'];
        name_img2 = [categories{n} '-' num2str(allpos(i,2)) '.jpg'];
        img1 = single(imread(name_img1))./255;
        img2 = single(imread(name_img2))./255;
        im1 = reshape(img1,[1,size(img1,1)*size(img1,2)]);
        im2 = reshape(img2,[1,size(img2,1)*size(img2,2)]);
        c1 = im1-im2;
        c2 = sum(c1.^2)/length(im1);
        Simis(i) = 1-c2;
    end
allSimis(n,:) = Simis;
clear Simis;
end
%%
%stats
data2test = allSimis([1,2,3,5],:);
[p,t,st] = anova1(data2test')

[p,t,st] = anova1(allSimis')
datadir = 'D:\Stanford_infant_EEG\fLoc-master\';
save(fullfile(datadir,'image_similarity_values.mat'),'allSimis');

%% boxplot
datadir = 'D:\Stanford_infant_EEG\fLoc-master\';
load(fullfile(datadir,'image_similarity_values.mat'));

mycolor = [1,0,0;...
    0.929,0.694,0.125;...
    0.466,0.674,0.188;...
    0.4,0.4,0.4;...
    0,0.447,0.741];
mycolor_inv = flipud(mycolor);

figure('color','w')
% bb = boxplot(allSimis','symbol','k+','OutlierSize',4,'notch','on'); % datapoits in x rows and group in columns

bb = boxplot(allSimis','symbol','k+','OutlierSize',4,'notch','on','Colors',mycolor); % datapoits in x rows and group in columns

h = findobj(gca,'Tag', 'Box');
for i = 1:length(h)
    boxEdges = findobj(h(i), 'type', 'line', 'Tag', 'Box');
    set(boxEdges, 'Color', mycolor_inv(i,:));
    set(boxEdges, 'LineWidth', 1);
end

outliers = findobj(bb, 'Tag', 'Outliers');
% Set the color of the outliers
for i = 1:length(outliers)
    set(outliers(i), 'MarkerEdgeColor', mycolor(i,:));
end

set(bb,'LineWidth',0.6) %0.4

ylim([0.9 1])
yticks([0.9 1])
xlim([0 6])
box off

% % Apply colors
% h = findobj(gca,'Tag','Box');
% for jj=1:length(h)
%    patch(get(h(jj),'XData') ,get(h(jj),'YData'),mycolor_inv(jj,:));
% end
set(gca,'LineWidth',0.5,'xtick',[])
fontsize(gca,7,'points')
%save
XimPrintsize = 2.6;
YimPrintsize = 2.6;
set(gcf,'PaperUnits','centimeters');
set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
figpath = 'D:\Stanford_infant_EEG\fLoc-analysis\';
figname = 'ImageSimilarity-boxplot-v2.eps';
print(gcf,fullfile(figpath,figname),'-depsc','-r600');


%% mean image - not used
clc;close all;clear

imgpath = 'D:\Stanford_infant_EEG\fLoc-master\stimuli\';
cd(imgpath);

categories = {'adult' 'limb' 'corridor' 'word' 'car'};
ncat = length(categories);

for n = 1:ncat
    newpath = fullfile(imgpath,categories{n});
    cd(newpath)
    files = dir('*.jpg');
    fLen = size(files,1);
    for f = 1:fLen
        AllImgs(f,:,:) = imread(files(f).name);
    end
    mImgs = squeeze(mean(AllImgs,1));
    figure('color','w')
    imagesc(mImgs/255,[0,1])
    colormap('gray')
    axis off
    axis square

    XimPrintsize = 2;
    YimPrintsize = 2;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
    figpath = 'D:\Stanford_infant_EEG\fLoc-analysis\';
    figname = ['meanImg_' categories{n} '.eps'];
    print(gcf,fullfile(figpath,figname),'-depsc','-r600');
%     close all;
end

%save colorbar
XimPrintsize = 3;
YimPrintsize = 3;
set(gcf,'PaperUnits','centimeters');
set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
figpath = 'D:\Stanford_infant_EEG\fLoc-analysis\';
figname = ['meanImg_colorbar.eps'];
print(gcf,fullfile(figpath,figname),'-depsc','-r600');

%% spatial frequency
close all;clc;clear
imgpath = 'D:\Stanford_infant_EEG\fLoc-master\stimuli\';
cd(imgpath);

categories = {'adult' 'limb' 'corridor' 'word' 'car'};
ncat = length(categories);
AmpSpec2D = [];
storeLogAmpSpec = [];
storeLogAmpSpec_circ = [];
for n = 1:ncat
    newpath = fullfile(imgpath,categories{n});
    cd(newpath)
    files = dir('*.jpg');
    fLen = size(files,1);
    for f = 1:fLen
%         AllImgs(f,:,:) = imread(files(f).name);
        im = imread(files(f).name);
        %amplitude spectrum + normalize it
        AmpSpec = fftshift(abs(fft2(im)) ./ size(im,1)^2);        
        AmpSpec2D(:,:,f) = AmpSpec;
        %log10 of amp spectrum
        logAmpSpec_circ = log10(circavg(AmpSpec));
        logAmpSpec = log10(AmpSpec);
        storeLogAmpSpec(:,:,f) = logAmpSpec;
        storeLogAmpSpec_circ(:,f,n) = squeeze(logAmpSpec_circ);
    end
    %mean log-amp-spec
    meanLogAmpSpec = mean(storeLogAmpSpec,3);
    figure('color','w')
    imagesc(meanLogAmpSpec,[-3 0]);
    colormap(jet)
    axis off
    axis square
    %save
    XimPrintsize = 2;
    YimPrintsize = 2;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
    figpath = 'D:\Stanford_infant_EEG\fLoc-analysis\';
    figname = ['meanLogAmpSpectrum_' categories{n} '.eps'];
    print(gcf,fullfile(figpath,figname),'-depsc','-r600');
    c = colorbar;

%     close all;
end

%save colorbar
c = colorbar;
c.FontSize = 8;
c.LineWidth = 0.5;
XimPrintsize = 3;
YimPrintsize = 3;
set(gcf,'PaperUnits','centimeters');
set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
figpath = 'D:\Stanford_infant_EEG\fLoc-analysis\';
figname = ['meanLogAmpSpectrum_colorbar.eps'];
print(gcf,fullfile(figpath,figname),'-depsc','-r600');

%% spatial frequency amplitude plot
close all;clc;clear
imgpath = 'D:\Stanford_infant_EEG\fLoc-master\stimuli\';
cd(imgpath);
categories = {'adult' 'limb' 'corridor' 'word' 'car'};
ncat = length(categories);
storeLogAmpSpec = [];
storeLogAmpSpec_circ = [];
for n = 1:ncat
    newpath = fullfile(imgpath,categories{n});
    cd(newpath)
    files = dir('*.jpg');
    fLen = size(files,1);
    for f = 1:fLen
%         AllImgs(f,:,:) = imread(files(f).name);
        im = imread(files(f).name);
        %amplitude spectrum + normalize it
        AmpSpec = fftshift(abs(fft2(im)) ./ size(im,1)^2);        
        AmpSpec2D(:,:,f) = AmpSpec;
        %log10 of amp spectrum
        logAmpSpec_circ = log10(circavg(AmpSpec));
        storeLogAmpSpec_circ(:,f,n) = squeeze(logAmpSpec_circ);
    end
end

data = storeLogAmpSpec_circ;
mm = squeeze(mean(data,2));
ss = squeeze(std(data,[],2));
mycolor = [1,0,0;...
    0.929,0.694,0.125;...
    0.466,0.674,0.188;...
    0.4,0.4,0.4;...
    0,0.447,0.741];

figure('color','w');
hold on;
for n = 1:ncat
%     plot(log10([1:512]),mm(:,n),'color',mycolor(n,:),'linewidth',1.5);
    patch([log10([1:512]) fliplr(log10([1:512]))],[mm(:,n)-ss(:,n); flipud(mm(:,n)+ss(:,n))]',...
        mycolor(n,:),'FaceAlpha',0.5,'EdgeColor','none');
end

% xlabel('log(frequency)')
% ylabel('log(amplitude)')

set(gca,'LineWidth',0.6,'xtick',[])
fontsize(gca,7,'points')
%save
XimPrintsize = 2.6;
YimPrintsize = 2.6;
set(gcf,'PaperUnits','centimeters');
set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
figpath = 'D:\Stanford_infant_EEG\fLoc-analysis\';
figname = 'SpatialFrequency-lineplot.eps';
print(gcf,fullfile(figpath,figname),'-depsc','-r600');

%%
data = storeLogAmpSpec_circ;
mm = squeeze(mean(data,1));

allpos = nchoosek([1:5],2);
nrow = size(allpos,1);
h = [];
p = [];
ksstat = [];
for n = 1:nrow
    ix1 = allpos(n,1);
    ix2 = allpos(n,2);
    [h(n), p(n), ksstat(n)] = kstest2(mm(:,ix1), mm(:,ix2));
end
thresh = 0.001/nrow;
p < thresh

%% image sequence- not used
close all;clc;clear
imgpath = 'D:\Stanford_infant_EEG\fLoc-master\stimuli\';
cd(imgpath);

categories = {'adult' 'limb' 'corridor' 'word' 'car'};
ncat = length(categories);
nimgs = 5;
imgdir = 'D:\Stanford_infant_EEG\fLoc-master\';
for n = 2:ncat
    newpath = fullfile(imgpath,categories{n});
    cd(newpath)
    files = dir('*.jpg');
    fLen = size(files,1);
    idx = datasample([1:fLen],nimgs);
    for ii = 1:nimgs
        tempname = files(idx(ii)).name;
        tempimg = imread(tempname);
        figure
        imagesc(tempimg);
        colormap('gray')
        set(gca,'xtick',[],'ytick',[],'LineWidth',0.3)
        axis square
        XimPrintsize = 1;
        YimPrintsize = 1;
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'paperposition',[10 10 XimPrintsize YimPrintsize]);
        figpath = 'D:\Stanford_infant_EEG\fLoc-analysis\';
        figname = [tempname(1:end-4) '.eps'];
        print(gcf,fullfile(figpath,figname),'-depsc','-r300');
    end
end
        
close all
    




        
        
