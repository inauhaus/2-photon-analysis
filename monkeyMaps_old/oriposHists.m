
figure, 

s = min([xsize ysize],[],2);
id = find(~isnan(s));
subplot(2,1,1)
plotBarGraph(s(id),0,[0:.1:2])
%set(gca,'Xtick',[0 90 180])
xlim([.3 1.8])
xlabel('RF size; 2sigma; minor axis')
axis square

id = find(~isnan(BWdiff));
subplot(2,1,2)
plotBarGraph(BWdiff(id),1,[-1.5:.15:1.5])
xlim([-1.5 1.5])
set(gca,'Xtick',[-1:.5:1])
xlabel('log2')
axis square
