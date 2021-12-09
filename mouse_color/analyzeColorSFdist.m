function analyzeColorSFdist(sftcAll_0,sftcAll_45,sftcAll_90,sftcAll_135,pSthresh)


Mpk =  max(sftcAll_0');
Spk =  max(sftcAll_90');
pSAll = Spk./(Spk+Mpk);

%Use all data points to plot %M histogram
histdom = linspace(0,1,13);
figure,histogram(1-pSAll,histdom,'FaceColor',[.5 .5 .5]), 
set(gca,'XTick',[0 1/3 2/3 1],'Tickdir','out')
xlabel('%M')


%Now cull the population to includes balanced S and M
idB = find(pSAll>pSthresh(1) & pSAll<pSthresh(2)); %Id cells that have input from both S and M
sftcAll_0 = sftcAll_0(idB,:);
sftcAll_90 = sftcAll_90(idB,:);
sftcAll_45 = sftcAll_45(idB,:);
sftcAll_135 = sftcAll_135(idB,:);

Mpk =  max(sftcAll_0');
Spk =  max(sftcAll_90');
Colpk =  max(sftcAll_135');
Lumpk =  max(sftcAll_45');


%[pColorLowSFAll pColorHiSFAll] = organizeF1vectors(rst45_135,rst0_90)


pColorAllSF =  sftcAll_135./(sftcAll_45+sftcAll_135);

pColorLowSFAll = pColorAllSF(:,1);
pColorMidSFAll = pColorAllSF(:,3);
pColorHiSFAll = pColorAllSF(:,5);
pColorPeakAll = Colpk./(Colpk+Lumpk);


figure, 

subplot(4,1,1),
histogram(pColorLowSFAll,histdom,'FaceColor',[.5 .5 .5])
set(gca,'XTick',[0 1/3 2/3 1],'Tickdir','out')
xlabel('percent color; LowSF')
med = median(pColorLowSFAll)
hold on
plot([med med],[-0 10],'b')
xlim([0 1])

subplot(4,1,2),
histogram(pColorMidSFAll,histdom,'FaceColor',[.5 .5 .5])
set(gca,'XTick',[0 1/3 2/3 1],'Tickdir','out')
xlabel('percent color; MidSF')
med = median(pColorMidSFAll)
hold on
plot([med med],[-0 10],'r')
xlim([0 1])

subplot(4,1,3),
histogram(pColorHiSFAll,histdom,'FaceColor',[.5 .5 .5])
set(gca,'XTick',[0 1/3 2/3 1],'Tickdir','out')
xlabel('percent color; HiSF')
med = median(pColorHiSFAll)
hold on
plot([med med],[-0 10],'r')
xlim([0 1])

subplot(4,1,4),
histogram(pColorPeakAll,histdom,'FaceColor',[.5 .5 .5])
set(gca,'XTick',[0 1/3 2/3 1],'Tickdir','out')
xlabel('percent color; PeakSF')
med = mean(pColorPeakAll);
hold on
plot([med med],[-0 10],'r')
xlim([0 1])


%Plot mean tuning curve and mean %color tuning curve
ma = max([sftcAll_45 sftcAll_135],[],2);
%ma = max([sftcAll_45],[],2);
ma = ones(size(ma));
%ma = sftcAll_135(:,1);
sftcAll_45norm = sftcAll_45./(ma*ones(1,6));
sftcAll_135norm = sftcAll_135./(ma*ones(1,6));

sfdom = [.0125 .025 .05 .1 .2 .4];

figure,
subplot(2,1,1)
semilogx(sfdom,nanmean(sftcAll_45norm),'o-k'), hold on, semilogx(sfdom,nanmean(sftcAll_135norm),'o-r')
xlabel('Spatial frequency (c/deg)')
ylabel('dF/F')
set(gca,'XTick',sfdom)
xlim([.01 .45])

subplot(2,1,2)
semilogx(sfdom,nanmedian(pColorAllSF),'o-k')
xlabel('Spatial frequency (c/deg)')
ylabel('Color selectivity')
set(gca,'XTick',sfdom)
xlim([.01 .45])
hold on
plot([sfdom(1) sfdom(end)],[0.5 0.5],'--k')




id = find(pColorPeakAll>0 & pColorPeakAll<1);
pColorPeakAll = pColorPeakAll(id);
lumID = find(pColorPeakAll<.33);
colorID = find(pColorPeakAll>.66);
lumcolorID = find(pColorPeakAll>.33 & pColorPeakAll<.66);

length(lumID)/length(pColorPeakAll)
length(lumcolorID)/length(pColorPeakAll)
length(colorID)/length(pColorPeakAll)

sftcAll_45 = sftcAll_45(id,:);
sftcAll_135 = sftcAll_135(id,:);

figure,
subplot(2,1,1)
sfdum = nanmean(sftcAll_45(lumID,:));
plot(sfdum/max(sfdum),'k')
hold on
sfdum = nanmean(sftcAll_45(lumcolorID,:));
plot(sfdum/max(sfdum),'b')
hold on
sfdum = nanmean(sftcAll_45(colorID,:));
plot(sfdum/max(sfdum),'r')
legend('lum','color-lum','color')
title('lum gratings')

subplot(2,1,2)
sfdum = nanmean(sftcAll_135(lumID,:));
plot(sfdum/max(sfdum),'k')
hold on
sfdum = nanmean(sftcAll_135(lumcolorID,:));
plot(sfdum/max(sfdum),'b')
hold on
sfdum = nanmean(sftcAll_135(colorID,:));
plot(sfdum/max(sfdum),'r')
legend('lum','color-lum','color')
title('color gratings')
