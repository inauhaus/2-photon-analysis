function [x,f] = gaussfitter2Drot

global RF

x0 = gaussfitguess2Drot;

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle2Drot',x0,options);
[x,f] = fminsearch('gaussfitter_handle2Drot',x0);
