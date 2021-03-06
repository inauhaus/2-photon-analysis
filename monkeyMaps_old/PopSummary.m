%This was made because it specifically is what I want for the supplemental
%in the paper

%% Cell-by-cell

sfBWhold = sfBW;
thresh = 1/sqrt(2);
id = find(LPness>thresh | isinf(LPness));  %bandwidth set to 0.7 of height
sfBW(id) = NaN;
%sfpref(id) = NaN;

orisighold = orisig;
id = find(orisig>60);
orisig(id) = NaN;
sfBW(id) = NaN;
sfpref(id) = NaN;

%id = find(round(sfpref*100)/100 == .5);
sfBWLin2 = (fhi-sfpref)*2;
%sfBWLin2 = sfBWLin;
id = find(fhi == sfpref);
sfBW(id) = NaN;
sfBWLin2(id) = NaN;  %high pass cells
sfpref(id) = NaN;

orisig(id) = NaN;
sfBW(id) = NaN;
sfpref(id) = NaN;

%sfBWLin2(id) = sfBWLin2(id)+.5;
%sfBWLin2 = 1./RFsize;

sfpref_oct = log2(sfpref);

figure

id = find(~isnan(orisig) & ~isnan(sfBW));
subplot(3,3,4), scatter(orisig(id),sfBW(id),'.k')
xlabel('ori fitted sigma'),ylabel('sf bandwidth (octaves)')
[R p] = corrcoef(orisig(id),sfBW(id));
title(['R = ' num2str(R(1,2)) '  p = ' num2str(p(1,2)) '; N = ' num2str(length(id))])
xlim([0 70]),ylim([0 3.5])
axis square

id = find(~isnan(sfpref_oct) & ~isnan(sfBW));
subplot(3,3,5), scatter(sfpref_oct(id),sfBW(id),'.k')
xlabel('sfpref (octaves)'),ylabel('sf bandwidth (octaves)')
[R p] = corrcoef(sfpref_oct(id),sfBW(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2)) '; N = ' num2str(length(id))])
xlim([-1 4]), ylim([0 3.5]), axis square


id = find(~isnan(sfBW));
subplot(3,3,6),
plotBarGraph(sfBW(id),1,[0:.2:4])
xlabel('sf BW (octaves)')
xlim(([0 3.5]))
axis square

id = find(~isnan(sfBWLin2));
subplot(3,3,9),
plotBarGraph(sfBWLin2(id),0,[0:.3:6])
xlabel('sf BW (cyc/deg)')
title(['mean = ' num2str(mean(sfBWLin2(id))) '  sig = ' num2str(std(sfBWLin2(id))) '; N = ' num2str(length(id))])
xlim([0 6])
axis square

id = find(~isnan(orisig) & ~isnan(sfBWLin2));
subplot(3,3,7), scatter(orisig(id),sfBWLin2(id),'.k')
xlabel('ori fitted sigma'),ylabel('sf bandwidth (cyc/deg)')
[R p] = corrcoef(orisig(id),sfBWLin2(id));
title(['R = ' num2str(R(1,2)) '  p = ' num2str(p(1,2)) '; N = ' num2str(length(id))])
xlim([0 70]),ylim([0 6])
axis square

id = find(~isnan(sfpref_oct) & ~isnan(sfBWLin2));
subplot(3,3,8), scatter(2.^sfpref_oct(id),sfBWLin2(id),'.k')
xlabel('sfpref (cyc/deg)'),ylabel('sf bandwidth (cyc/deg)')
[R p] = corrcoef(2.^sfpref_oct(id),sfBWLin2(id));
title(['R = ' num2str(R(1,2)) '  p = ' num2str(p(1,2)) '; N = ' num2str(length(id))])
xlim([0 7]), ylim([0 6]), axis square

id = find(~isnan(orisig));
subplot(3,3,1),
plotBarGraph(orisig(id),0,0:3:70)
xlabel('ori sig (deg)')
title(['mean = ' num2str(mean(orisig(id))) '  sig = ' num2str(std(orisig(id))) '; N = ' num2str(length(id))])
xlim([0 70])
axis square

id = find(~isnan(sfpref_oct));
subplot(3,3,2), 
plotBarGraph(sfpref_oct(id),1,[-1:.25:3])
xlabel('sf pref (octaves)')
axis square

%% oripref vs. sfpref
figure
id = find(~isnan(sfpref_oct) & ~isnan(opref));
scatter(opref(id),sfpref_oct(id),'.k')
ylabel('sfpref (octaves)'),xlabel('orientation preference (degrees)')
[R p] = corrcoef(opref(id),sfpref_oct(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2)) '; N = ' num2str(length(id))])
hold on

odom = 0:45:180; do = odom(2)-odom(1);
clear prc50 prc25 prc75
for i = 1:length(odom)
    o1 = odom(i)-do/2; o2 = odom(i)+do/2;
    
    if o1>0 & o2<180
        id2 = find(opref(id)<o2 & opref(id)>o1);
    end
    
    if o1<0
        o1 = o1+180;
        id2 = find(opref(id)<o2 | opref(id)>o1);
    end
    if o2>180
        o2 = o2-180;
        id2 = find(opref(id)<o2 | opref(id)>o1);
    end       
    
    sfdist{i} = sfpref_oct(id(id2));
    

    %prc50(i) = mean(sfpref_oct(id(id2)));
    prc50(i) = prctile(sfdist{i},50);
    prc25(i) = prctile(sfdist{i},25);
    prc75(i) = prctile(sfdist{i},75);    
end

plot(odom,prc50,'o-b')
plot(odom,prc25,'o--b')
plot(odom,prc75,'o--b')
%xlim([-1 4]), ylim([0 3.5]), axis square

h = zeros(length(odom)-1,length(odom)-1)*NaN; p = zeros(length(odom)-1,length(odom)-1)*NaN;
for i = 1:length(odom)-1
    for j = i+1:length(odom)-1
        [h(i,j) p(i,j)] = ttest2(sfdist{i},sfdist{j});
    end
end

%% Compare selectivity vs. map

figure
id = find(~isnan(orisig) & orisig<60 & ~isnan(oriCV_Pop));
subplot(3,2,1), scatter(oriCV_Pop(id),orisig(id),'.k')
xlabel('orimap CV'),ylabel('ori cell Sigma')
[R p] = corrcoef(oriCV_Pop(id),orisig(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])
axis square

id = find(~isnan(orisig) & ~isnan(oriCV_pks));
subplot(3,2,2), scatter(sqrt(oriCV_pks(id)),orisig(id),'.k')
xlabel('orimap circular sig (peaks)'),ylabel('ori cell Sigma')
[R p] = corrcoef(sqrt(oriCV_pks(id)),orisig(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])
axis square

%as a function of sfmap "variance"
id = find(~isnan(sfvar_Pop) & ~isnan(sfBW));
subplot(3,2,3), scatter(sqrt(sfvar_Pop(id)),sfBW(id),'.k')
xlabel('sf map stdev'),ylabel('sf cell bandwidth')
[R p] = corrcoef(sqrt(sfvar_Pop(id)),sfBW(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])
ylim([0 3.5]), xlim([.5 1.5])
axis square

%as a function of sfmap "variance" (peaks)
id = find(~isnan(sfvar_pks) & ~isnan(sfBW));
subplot(3,2,4), scatter(sqrt(sfvar_pks(id)),sfBW(id),'.k')
xlabel('sf map stddev (peaks)'),ylabel('sf cell bandwidth')
[R p] = corrcoef(sqrt(sfvar_pks(id)),sfBW(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])
ylim([0 3.5])
axis square

id = find(~isnan(sfvar_Pop) & ~isnan(oriCV_Pop));
subplot(3,2,5), scatter(sqrt(sfvar_Pop(id)),oriCV_Pop(id),'.k')
ylabel('ori map CV'),xlabel('sf map stddev')
[R p] = corrcoef(oriCV_Pop(id),sqrt(sfvar_Pop(id)));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])
xlim([.5 1.5])
axis square

id = find(~isnan(sfvar_pks) & ~isnan(oriCV_pks));
subplot(3,2,6), scatter(sqrt(sfvar_pks(id)),sqrt(oriCV_pks(id)),'.k')
ylabel('ori map circular sig (peaks)'),xlabel('sf map stddev (peaks)')
[R p] = corrcoef(oriCV_pks(id),sqrt(sfvar_pks(id)));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])
axis square

%%

figure

id = find(~isnan(sfpref_oct) & ~isnan(F1F0));
subplot(2,2,2), scatter(sfpref_oct(id),F1F0(id),'.k')
xlabel('sf preference (octaves)'),ylabel('F1/F0')
[R p] = corrcoef(sfpref_oct(id),F1F0(id));
title(['R = ' num2str(R(1,2)) '  p = ' num2str(p(1,2)) '; N = ' num2str(length(id))])
xlim([-1.3 3.5]),ylim([0 3.5])
set(gca,'YTick',[0 1 2 3])
set(gca,'XTick',-1:3)
axis square

subplot(2,2,1),
plotBarGraph(F1F0(id),0,[-.0:.3:3.4])
xlabel('F1/F0')
set(gca,'XTick',[0 1 2 3])
xlim(([-.2 3.4]))
axis square

subplot(2,2,4),
plotBarGraph(sfpref_oct(id),0,[-1:.35:3.4])
xlabel('sfpref')
set(gca,'XTick',-1:3)
xlim(([-1.3 3.5]))
axis square


%% LMS tuning
% 
% figure
% for i = 1:3
%     N = length(find(~isnan(orisig{i})))
%     id = find(~isnan(orisig{i}));
%     subplot(3,2,i*2-1),
%     plotBarGraph(orisig{i}(id),0,0:3:70)
%     xlabel('ori sig (deg)')
%     title(['median = ' num2str(median(orisig{i}(id))) '; mean = ' num2str(mean(orisig{i}(id))) '  sig = ' num2str(std(orisig{i}(id))) '; N = ' num2str(length(id))])
%     xlim([0 70])
%     axis square
% end
% 
% 
% 
% for i = 1:3
%     id = find(~isnan(sfpref{i}));
%     subplot(3,2,i*2)
%      plotBarGraph(sfpref{i}(id),1,[-1:.25:5])    
%     xlabel('sf pref (octaves)')
%     title(['geomean = ' num2str(geomean(sfpref{i}(id))) '  sig = ' num2str(std(sfpref{i}(id))) ' ; N = ' num2str(length(id))])
%     %xlim([0 70])
%     axis square
% end