function [x,f,g] = gaussfitter2

global RF dom

x0 = gaussfitguess2;


% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);

[x,f] = fminsearch('gaussfitter_handle2',x0);

% lb = x0-[1.5 1 .1 .1];
% ub = x0+[1.5 1 .2 .1];
% [x,f] = fmincon('gaussfitter_handle',x0,[],[],[],[],lb,ub);

xc = x(1);
sig = x(2);

A = x(3);
B = x(4);

d = (dom-xc).^2;
g = A*exp(-d/(2*sig^2))+B;
