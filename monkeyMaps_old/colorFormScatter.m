

figure,

subplot(1,2,1)
scatter(lumpref{1},log2(sfpref{2}),'.')
id = find(~isnan(lumpref{1}.*log2(sfpref{2})));
[r p] = corrcoef(lumpref{1}(id),log2(sfpref{2}(id)));
title(['r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2))])
xlabel('log_2[(L-M)/(L+M)]')
ylabel('sfpref (octaves)  (L-M)')

subplot(1,2,2)
scatter(lumpref{1},log2(sfpref{3}),'.')
id = find(~isnan(lumpref{1}.*log2(sfpref{3})));
[r p] = corrcoef(lumpref{1}(id),log2(sfpref{3}(id)));
title(['r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2))])
xlabel('log_2[(L-M)/(L+M)]')
ylabel('sfpref (octaves)  (L+M)')