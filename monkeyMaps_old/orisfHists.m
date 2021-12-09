
figure, 

id = find(~isnan(opref));
subplot(3,1,1)
plotBarGraph(opref(id),0,[0:10:180])
set(gca,'Xtick',[0 90 180])
xlim([-10 190])
xlabel('ori preference (deg)')
axis square

id = find(~isnan(sfpref_oct));
subplot(3,1,2)
plotBarGraph(sfpref_oct(id),1,[-1:.25:3])
xlim([-1.5 3])
xlabel('sf pref (octaves)')
axis square

id = find(~isnan(F1F0));
subplot(3,1,3)
plotBarGraph(F1F0(id),0,[0:.15:3])
set(gca,'Xtick',[0 1 2])
xlim([-.1 2.1])
xlabel('sf pref (octaves)')
axis square