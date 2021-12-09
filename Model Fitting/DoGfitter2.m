function x = DoGfitter2

x0 = DoGfitguess2;

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);
[x,f] = fminsearch('DoGfitter2_handle',x0);

