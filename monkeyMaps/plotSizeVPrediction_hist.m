function plotSizeVPrediction_hist(Rsize,Rsizepred)

%%

Rsize = log2(Rsize); Rsizepred = log2(Rsizepred);
 
predratio = Rsize-Rsizepred;
[h] = histogram(predratio,10,'FaceColor',[1 1 1]);
%axis square
set(gca,'TickDir','out')
xtck = h.BinEdges(1:2:end);
set(gca,'XTick',h.BinEdges(1:2:end)); 
set(gca,'XTickLabel',round(2.^xtck*10)/10);
%xlim([-2 2])
%xlabel('log2(Actual RF Size) - log2(Prediction of RF size)')
xlabel('actualWidth/predictedWidth')
title(['geomean(ratio) = ' num2str(2.^nanmean(predratio))])
