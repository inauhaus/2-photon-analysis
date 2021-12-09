function [x fit varacc] = ColorPlane(x,y,f)

global data domx domy

%%%%%%%%%%%
% mi = min(f);
% f = f-mi;
% E = max(f);
% f = f/E;
%%%%%%%%%%%%

data = f;
domx = x;
domy = y;

%% Initial guess

H = [abs(x(:)) abs(y(:))];
x0 = inv(H'*H)*H'*f(:);
%x0 = [1 1]

%% Fit the parameters

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);

[x,f] = fminsearch('ColorPlanefitter_handle',x0);

alpha = x(1);
beta = x(2);

fit = abs(alpha*domx + beta*domy);

varacc = (var(data(:))-var(data(:)-fit(:)))/var(data(:));






