function [param ffit varaccount] = Expfit4(domain,f)

%3 is an extension of 2, which only varies the decay constant and
%amplitude, while assuming it asymptotes to 'base'.

global yy xx


%%%search%%%
xx = domain;
yy = f;
param = expfitter4;

%%%%%%%%%%%

domu = unique(xx);

ffit = exp(-param(2)*(domu).^param(1));


%%%%
id = find(isnan(yy.*xx));
yy(id) = [];
xx(id) = [];
expect = exp(-param(2)*(xx).^param(1)) + param(1);
varaccount = (var(yy)-var(yy-expect))/var(yy);
%%%%