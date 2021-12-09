function plotSizeVPrediction_scatter(Rsize,Rsizepred,varargin)

%%
if isempty(varargin)
    hdom = logspace(log10(1/(2*getparam('max_sf'))),log10(1/(2*getparam('min_sf'))),2*getparam('n_sfreq'));
    hdom = log2(hdom);
else
    hdom = varargin{1};
    hdom = log2(hdom);
end

loglog(Rsizepred,Rsize,'.k'), xlim(2.^[hdom(1) hdom(end)]), ylim(2.^[hdom(1) hdom(end)])
Rsize = log2(Rsize); Rsizepred = log2(Rsizepred);
id = find(~isnan(Rsize.*Rsizepred));
[r p] = corrcoef(Rsize(id),Rsizepred(id));
%xlim([-1 3]), ylim([-1 3])
ylabel('actual Rsize'), xlabel('Predicted Rsize')
title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2))])
axis square
hold on, plot([.01 8],[.01 8],'k')
