function [flo fpeak fhi] = getDoGparams(param,thresh)

dom = linspace(0,30,1000);
ffit1 = exp(-dom.^2/(2*param(1).^2))*param(2);
ffit2 = exp(-dom.^2/(2*param(3).^2))*param(4);
ffit = ffit1 - ffit2 + param(5);

[ma id] = max(ffit);
fpeak = dom(id);

[ma idma] = max(ffit); 
mi = ffit(end);
idmi = length(ffit);

thresh = (ma-mi)*thresh + mi;
[dum fhidum] = min(abs(ffit(idma:idmi) - thresh));
fhi = dom(fhidum+idma-1);  %high cutoff sfreq
[dum flodum] = min(abs(ffit(1:idma) - thresh));
flo = dom(flodum);  %low cutoff sfreq