function [param ffit varaccount domu] = Decay_fit(domain,f,G0)

%This is used to model RF size vs. SF preference.  Standard model is 4sig =
%1/sf.  This allows for sig and sf to not asymptote to infinity.  e.g. size
%can have a zero crossing at really low values of sf.

global yy xx


%%%search%%%
xx = domain;
yy = f;
param = Decayfitter(G0);
%%%%%%%%%%%

domu = unique(xx);

ffit = 1./(domu*param(3)+param(1)) + param(2);


%%%%
id = find(isnan(yy.*xx));
yy(id) = [];
xx(id) = [];
expect = 1./(xx*param(3)+param(1)) + param(2);
varaccount = (var(yy)-var(yy-expect))/var(yy);
%%%%
