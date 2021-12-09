
figure, 

id = find(~isnan(LMphaseDiff{1}));
subplot(2,1,1)
plotBarGraph(LMphaseDiff{1}(id),0,[0:10:180])
set(gca,'Xtick',[0 90 180])
xlim([-10 190])
xlabel('phase difference (deg)')
axis square

id = find(~isnan(lumpref{1}));
subplot(2,1,2)
plotBarGraph(lumpref{1}(id),1,[-.5:.05:.5])
set(gca,'Xtick',[-.5 0 .5])
xlim([-.5 .5])
if strcmp(getparam('colorspace'),'LMS')
    xlabel('log_2(L/M)')
elseif strcmp(getparam('colorspace'),'DKL')
    xlabel('log_2((L-M)/(L+M))')
end

axis square