function [param ffit varaccount] = Sigfit2(domain,f,varargin)

global yy xx initguess

if ~isempty(varargin)
    initguess = varargin{1};
else 
    initguess = [];
end

%%%search%%%
xx = domain;
yy = f;
param = sigfitter2;
%%%%%%%%%%%

xc = param(1);
sig = param(2);
A = 1;
B = 0;
d = xx-xc;
ffit = A./(1 + exp(-d*sig)) + B;

%%%%
% id = find(isnan(yy.*xx));
% yy(id) = [];
% xx(id) = [];
% expect = param(1)./(1 + exp(-xx*param(2)));
varaccount = (var(yy)-var(yy-ffit))/var(yy);

%%%%
