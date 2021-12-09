function plotSFparamVColor(tcP,pColordum,pSdum,pSlohi)


%%
idB = find(~isnan(tcP.sfprefLum.*tcP.sfprefColor.*tcP.sfCoMLum.*tcP.sfCoMColor.*tcP.sfhcoLum.*tcP.sfhcoColor.*tcP.sfBPLum.*tcP.sfBPColor));
idx = find(pSdum(idB)>pSlohi(1) & pSdum(idB)<pSlohi(2));
%idx = find(pSAll(idB)<.2);
idB = idB(idx);

%Cull the poplation
pColordum = pColordum(idB);
vars = fields(tcP);
for i = 1:length(vars)
    eval(['tcP.' vars{i} ' = tcP.' vars{i} '(idB);'])
end


idLum = find(pColordum>-10 & pColordum<1/3);
idColorLum = find(pColordum>1/3 & pColordum<2/3);
idColor = find(pColordum>2/3 & pColordum<10);


%idB = find(pSAll<.1 | pSAll>.6);


sfdom = logspace(log10(.0125),log10(.4),5);

figure,
subplot(5,2,1)
loglog([sfdom(1) sfdom(end)],[sfdom(1) sfdom(end)],'k'), hold on
loglog(tcP.sfprefLum(idColorLum),tcP.sfprefColor(idColorLum),'.r'),
loglog(tcP.sfprefLum(idLum),tcP.sfprefColor(idLum),'.k'), 
loglog(tcP.sfprefLum(idColor),tcP.sfprefColor(idColor),'or'),
%loglog(tcP.sfprefLum(idB),tcP.sfprefColor(idB),'.k'),
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf peak loc; Lum ')
ylabel('Sf peak loc; Color ')
axis square
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])
id = find(~isnan(tcP.sfprefLum.*tcP.sfprefColor));
[r p] = corrcoef(log2(tcP.sfprefLum),log2(tcP.sfprefColor));
title(['r=' num2str(round(r(1,2)*100)/100) '; p=' num2str(p(1,2))])

subplot(5,2,2)
histogram(log2(tcP.sfprefLum./tcP.sfprefColor),[-4:.25:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(tcP.sfprefLum./tcP.sfprefColor));
[h p] = ttest(log2(tcP.sfprefLum./tcP.sfprefColor));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(5,2,3)
loglog([sfdom(1) sfdom(end)],[sfdom(1) sfdom(end)],'k'), hold on
loglog(tcP.sfCoMLum(idColorLum),tcP.sfCoMColor(idColorLum),'.r'), hold on, 
loglog(tcP.sfCoMLum(idLum),tcP.sfCoMColor(idLum),'.k'), 
loglog(tcP.sfCoMLum(idColor),tcP.sfCoMColor(idColor),'or'), 
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf CoM; Lum ')
ylabel('Sf CoM; Color ')
axis square
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])
id = find(~isnan(tcP.sfCoMLum.*tcP.sfCoMColor));
[r p] = corrcoef(log2(tcP.sfCoMLum),log2(tcP.sfCoMColor));
title(['r=' num2str(round(r(1,2)*100)/100) '; p=' num2str(p(1,2))])

subplot(5,2,4)
histogram(log2(tcP.sfCoMLum./tcP.sfCoMColor),[-3:.25:3],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(tcP.sfCoMLum./tcP.sfCoMColor));
[hyp p] = ttest(log2(tcP.sfCoMLum./tcP.sfCoMColor));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(5,2,5)
loglog([sfdom(1) sfdom(end) ],[sfdom(1) sfdom(end)],'k'), hold on,
n1 = .05*randn(size(idColorLum)); n2 = .05*randn(size(idColorLum));
loglog(2.^(log2(tcP.sfhcoLum(idColorLum))+n1),2.^(log2(tcP.sfhcoColor(idColorLum))+n2),'.r'),
n1 = .05*randn(size(idLum)); n2 = .05*randn(size(idLum));
loglog(2.^(log2(tcP.sfhcoLum(idLum))+n1),2.^(log2(tcP.sfhcoColor(idLum))+n2),'.k')
n1 = .05*randn(size(idColor)); n2 = .05*randn(size(idColor));
loglog(2.^(log2(tcP.sfhcoLum(idColor))+n1),2.^(log2(tcP.sfhcoColor(idColor))+n2),'or')
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf high pass cutoff; Lum ')
ylabel('Sf high pass cutoff; Color ')
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])
axis square
id = find(~isnan(tcP.sfhcoLum.*tcP.sfhcoColor));
[r p] = corrcoef(log2(tcP.sfhcoLum),log2(tcP.sfhcoColor));
title(['r=' num2str(round(r(1,2)*100)/100) '; p=' num2str(p(1,2))])

subplot(5,2,6)
histogram(log2(tcP.sfhcoLum./tcP.sfhcoColor),[-4:.5:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(tcP.sfhcoLum./tcP.sfhcoColor));
[hyp p] = ttest(log2(tcP.sfhcoLum./tcP.sfhcoColor));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(5,2,7)
plot([0 1],[0 1],'k'), hold on
scatter(tcP.sfBPLum(idColorLum),tcP.sfBPColor(idColorLum),'.r')
scatter(tcP.sfBPLum(idLum),tcP.sfBPColor(idLum),'.k')
scatter(tcP.sfBPLum(idColor),tcP.sfBPColor(idColor),'or')
xlabel('Sf bandpass; Lum ')
ylabel('Sf bandpass; Color ')
axis square
id = find(~isnan(tcP.sfBPLum.*tcP.sfBPColor));
[r p] = corrcoef(tcP.sfBPLum,tcP.sfBPColor);
title(['r=' num2str(round(r(1,2)*100)/100) '; p=' num2str(p(1,2))])

subplot(5,2,8)
h  = histogram(tcP.sfBPLum - tcP.sfBPColor,[-1:.1:1],'FaceColor',[1 1 1]); 
%plot(h.BinEdges(1:end-1),cumsum(h.Values/sum(h.Values)))
set(gca,'TickDir','out')
mu = nanmedian(tcP.sfBPLum - tcP.sfBPColor);
[hyp p] = ttest(tcP.sfBPLum - tcP.sfBPColor);
title(['mu=' num2str(mu) ' p=' num2str(p)])
xlabel('x-y')


subplot(5,2,9)
loglog([sfdom(1) sfdom(end) ],[sfdom(1) sfdom(end)],'k'), hold on
loglog(tcP.sfBWLum(idColorLum),tcP.sfBWColor(idColorLum),'.r')
loglog(tcP.sfBWLum(idLum),tcP.sfBWColor(idLum),'.k')
loglog(tcP.sfBWLum(idColor),tcP.sfBWColor(idColor),'or')
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf bandwidth; Lum ')
ylabel('Sf bandwidth; Color ')
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])
axis square
id = find(~isnan(tcP.sfBWLum.*tcP.sfBWColor));
[r p] = corrcoef(log2(tcP.sfBWLum),log2(tcP.sfBWColor));
title(['r=' num2str(round(r(1,2)*100)/100) '; p=' num2str(p(1,2))])


legend('unity','color-lum','lum','color')

subplot(5,2,10)
histogram(log2(tcP.sfBWLum./tcP.sfBWColor),[-4:.5:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(tcP.sfBWLum./tcP.sfBWColor));
[hyp p] = ttest(log2(tcP.sfBWLum./tcP.sfBWColor));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')


%SFfit_percYield = length(idB)/length(tcP.sfprefLum)