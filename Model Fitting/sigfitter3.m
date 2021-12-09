function [x,f,g] = sigfitter3

global RF xx initguess

if ~isempty(initguess)
    x0 = initguess;
else
    x0 = [0 1 0];
end
    
dim = length(RF);

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);

[x,f] = fminsearch('sigfitter_handle3',x0);

% lb = x0-[1.5 1 .1 .1];
% ub = x0+[1.5 1 .2 .1];
% [x,f] = fmincon('gaussfitter_handle',x0,[],[],[],[],lb,ub);

xc = x(1);
C = x(2);
B = x(3);
d = abs(xx-xc);
g = 5*d./(C+d) + B;


