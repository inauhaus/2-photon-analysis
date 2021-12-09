function [x,f,g] = expfitter

global RF

x0 = expfitguess;
dim = length(RF);

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('expfitter_handle',x0,options);
[x,f] = fminsearch('expfitter_handle',x0);

B = x(1);
alp = x(2);

A = x(3);

xx = 0:dim-1;

g = B*exp(-alp*xx)+A;
