function [param ffit varaccount,domu] = NLsoftsat(domain,f,G0)

%2 allows 'domain' vs. 'f' to be a scatter plot (i.e. domain is not
%monotonic)

global yy xx


%%%search%%%
xx = domain;
yy = f;
param = NLsoftsatfitter(G0);
%%%%%%%%%%%

domu = unique(xx); 

%domuPad = linspace(domu(1)-range(domu),domu(end)+range(domu),100);


ffit = param(1)*domu;
%ffit = ffit-param(2);
ffit(find(ffit>param(2))) = param(2);  %saturation level


ffit = smoothFit(domu,ffit,param,.1);


%%%%
id = find(isnan(yy.*xx));
yy(id) = [];
xx(id) = [];

expect = param(1)*xx;
%expect = xx;
expect(find(expect>param(2))) = param(2);  %saturation level


%expect = param(1)*phi(xx-param(2)).^param(3);
varaccount = (var(yy)-var(yy-expect))/var(yy);
%%%%


