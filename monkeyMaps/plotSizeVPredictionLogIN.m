function plotSizeVPredictionLogIN(Rsize,Rsizepred,scatterdom)

%%
%hdom = logspace(log10(1/(2*getparam('max_sf'))),log10(1/(2*getparam('min_sf'))),2*getparam('n_sfreq'));
id = find(~isnan(Rsize.*Rsizepred));
Rsize = Rsize(id);
Rsizepred = Rsizepred(id);

hdom = linspace(scatterdom(1),scatterdom(end),2*length(scatterdom)-1);

figure, 
subplot(2,2,1)
loglog(Rsizepred,Rsize,'.k'), xlim([scatterdom(1) scatterdom(end)]), ylim([scatterdom(1) scatterdom(end)])
set(gca,'XTick',scatterdom);
set(gca,'YTick',scatterdom);

% [param ffit dom varaccount] = LinetotalLS(log2(Rsizepred),log2(Rsize),[1 0]);
% param
% dom = log2(scatterdom);
% hold on
% loglog(2.^dom,2.^(dom*param(1)+param(2)))


id = find(~isnan(Rsize.*Rsizepred));
[r p] = corrcoef(Rsize(id),Rsizepred(id));
%xlim([-1 3]), ylim([-1 3])
ylabel('actual Rsize'), xlabel('Predicted Rsize')
title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2))])
axis square
hold on, plot([scatterdom(1) scatterdom(end)],[scatterdom(1) scatterdom(end)],'k')

subplot(2,2,4), 
predratio = Rsize-Rsizepred;
[h] = histogram(predratio,10,'FaceColor',[1 1 1]);
%axis square
set(gca,'TickDir','out')
xtck = h.BinEdges(1:2:end);
set(gca,'XTick',h.BinEdges(1:2:end)); 
set(gca,'XTickLabel',round(2.^xtck*10)/10);
%xlim([-2 2])
%xlabel('log2(Actual RF Size) - log2(Prediction of RF size)')
idisn = find(~isnan(predratio));
[h p] = ttest(predratio(idisn));
xlabel('actualWidth/predictedWidth')
title(['geomean(ratio) = ' num2str(nanmean(predratio)) '; p=' num2str(p)])

axis square

%%
subplot(2,2,3), 
h = histogram(Rsize,(hdom),'FaceColor',[1 0 0]);
axis square %If you put this at the end it screws up axes
set(gca,'TickDir','out')
xtck = h.BinEdges(1:2:end);
set(gca,'XTick',h.BinEdges(1:2:end)); 
set(gca,'XTickLabel',round(xtck*100)/100);
xlabel('actual width (2sig)')
title(['gmu = ' num2str(nanmean(Rsize)) '; sig = ' num2str(nanstd(Rsize)) 'oct'])
xlim([hdom(1) hdom(end)])

hold on
histogram(Rsizepred,(hdom),'FaceColor',[0 0 1]);
xtck = h.BinEdges(1:2:end);
set(gca,'XTick',h.BinEdges(1:2:end)); 
set(gca,'XTickLabel',round(xtck*100)/100);
xlabel('predicted width (2sig) = 1/(pi*BW/2)')
title(['gmu = ' num2str(nanmean(Rsizepred)) '; sig = ' num2str(nanstd(Rsizepred)) 'oct'])
set(gca,'TickDir','out')
axis square
