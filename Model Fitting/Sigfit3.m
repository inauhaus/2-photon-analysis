function [param ffit varaccount] = Sigfit3(domain,f,varargin)

global yy xx initguess


if ~isempty(varargin)
    initguess = varargin{1};
else 
    initguess = [];
end

%%%search%%%
xx = domain;
yy = f;
param = sigfitter3;
%%%%%%%%%%%

xc = param(1);
C = param(2);
B = param(3);
d = xx-xc;
ffit = 5*d./(C+abs(d)) + B;

%%%%
% id = find(isnan(yy.*xx));
% yy(id) = [];
% xx(id) = [];
% expect = param(1)./(1 + exp(-xx*param(2)));
varaccount = (var(yy)-var(yy-ffit))/var(yy);

%%%%
