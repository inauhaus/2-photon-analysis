function [param ffit varacc sigma] = Gaussfit2(x,f)

global RF dom

%%%%%%%%%%%
mi = min(f);
f = f-mi;
E = max(f);
f = f/E;
%%%%%%%%%%%%

%%%search%%%
RF = f;
dom = x;
param = gaussfitter2;

%%%%%%%%%%%

ffit = exp(-(dom-param(1)).^2/(2*param(2).^2));
ffit = param(3)*ffit + param(4);

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

ffit = ffit*E + mi;
param(3) = param(3)*E;
param(4) = param(4)*E+mi;


%% get sigma as 61% of max

xxI = linspace(0,180,1001);
xxI = xxI(1:end-1);

ffitI = exp((-xxI.^2)/(2*param(2)^2));

ffitI = ffitI-min(ffitI);
ffitI = ffitI/max(ffitI);

thresh = exp(-1/2);
thresh = 1/sqrt(2);
[dum id] = min(abs(ffitI-thresh));

sigma = (xxI(2)-xxI(1))*(id-1);



%%


