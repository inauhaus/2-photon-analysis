function param = expfitter4

global yy xx

%x0 = expfitguess1;
%ffit = exp(-param(2)*(domu).^param(1));
x0 = [1 0];

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('expfitter_handle',x0,options);
param = fminsearch('expfitter4_handle',x0);


