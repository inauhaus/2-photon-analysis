function [param ffit varacc ffitI domI] = GaussTLineFit(domain,f)

%Gaussian times a line.  i.e. the FT of a Gaussian derivative model

%%%%Normalize%%%%%%%
E = max(f);
f = f/E;

%% Compute initial guess %%

[dum,idx] = max(f);
muGuess = domain(idx);
sigGuess = range(domain)/4;
%sig = length(f)/8;

AmpGuess = max(f);
DCguess = 0;

%G = [idx length(f)/6 max(RF)-min(RF) min(RF)];
x0 = [muGuess sigGuess AmpGuess DCguess]; %initial guess for paramter search
%G = [length(f)/2 1 0 0];

%% Search %%%
global dom RF
dom = domain;
RF = f;

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);


[param,dum] = fminsearch('GaussTLineFitter_handle',x0);

%%%% Compute the fit %%%%%%%

ffit = exp(-(domain-param(1)).^2/(2*param(2).^2));
ffit = param(3)*ffit.*(domain+param(4));

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

ffit = ffit*E;
param(3) = param(3)*E;

%%%% Compute the interpolated fit %%%%%%%
N = 50;

domI = logspace(log10(domain(1)),log10(domain(end)),N);
ffitI = exp(-(domI-param(1)).^2/(2*param(2).^2));
ffitI = param(3)*ffitI.*(domI+param(4));


