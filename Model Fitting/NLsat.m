function [param ffit varaccount,domu] = NLsat(domain,f,G0)

%2 allows 'domain' vs. 'f' to be a scatter plot (i.e. domain is not
%monotonic)

global yy xx


%%%search%%%
xx = domain;
yy = f;
param = NLsatfitter(G0);
%%%%%%%%%%%

domu = unique(xx);

%domuPad = linspace(domu(1)-range(domu),domu(end)+range(domu),100);


%ffit = param(1)*domu;
ffit = domu-param(2);
ffit(find(ffit>param(1))) = param(1);  %saturation level
% 
% h = fspecial('gaussian',size(ffit),3);
% h = h/sum(h);
% ffit = ifft(fft(ffit).*abs(fft(h)));
% 
% id = find(domuPad>domu(1) & domuPad<domu(end));
% ffit = ffit(id);
% domu = domuPad(id);
% 
% 



%%%%
id = find(isnan(yy.*xx));
yy(id) = [];
xx(id) = [];

%expect = param(1)*xx;
expect = xx;
expect(find(expect>param(1))) = param(1);  %saturation leve


%expect = param(1)*phi(xx-param(2)).^param(3);
varaccount = (var(yy)-var(yy-expect))/var(yy);
%%%%
