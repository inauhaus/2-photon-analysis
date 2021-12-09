function [params ffitII domI domII pfit varaccfit repFitQuality] = analyzeLOGtuning(tc,dom,cothresh,varargin)

%INPUTS
%tc - tuning curve
%dom - domain.  e.g. spatial frequency 
%cothresh is percent of max for high pass cut off.  e.g. cothresh = 0.75; 

%OUTPUTS
%params - structure containing all the computed parameters
%ffitII - interpolated fit, using the domain in domII
%domI - slightly interpolated dom for improved Gaussian fit
%domII - highly interpolated dom to generate Gaussian fit
%pfit - p-value of correlation coefficient between fit and data
%%varaccfit - variance accounted for in the fit.

InterpNfit = 2; %for domI

InterpNplot = 5; %for domII

domI = logspace(log10(dom(1)),log10(dom(end)),length(dom)*InterpNfit-1);
domII = logspace(log10(dom(1)),log10(dom(end)),length(dom)*InterpNplot-1); %strictly for plotting a smooth curve

tcI = interp1(log2(dom),tc,log2(domI),'linear ')'; %Interpolate tuning curve before fitting

%[Gparam ffitI varaccfit] = Gfit(log2(domI),tcI');
[Gparam ffitI varaccfit] = Gfit_base0(log2(domI),tcI');
if isnan(ffitI(1))
    ffitI = zeros(size(ffitI));
end
GparamRaw = Gparam;
Gparam(1) = min([log2(dom(end)) Gparam(1)]); %If the peak is at > highest sf, use highest sf instead
Gparam(1) = max([log2(dom(1)) Gparam(1)]);  %If the peak is at < lowest sf, use lowest sf instead
BWoct = Gparam(2);  %bandwidth in octaves (sigma).  N.B. taking 2.^BW does not give you linear BW.  i.e. ~= sfhi-sflo

pref = 2.^Gparam(1); %peak location in Gaussian fit, in cyc/deg

%Compute linear bandwidth (sfhi-sflo) from log-based parameters
hcoE = 2.^Gparam(1)*(2.^BWoct); %Should get the same as hco below, but only if the fit ids the peak... i.e. if its bandpass. 
lcoE = 2.^Gparam(1)/(2.^BWoct); %Should get the same as lco below, but only if the fit ids the peak... i.e. if its bandpass. 
BW = hcoE-lcoE; 

%BW = 2*sqrt((2^(BWoct^2)-1)*2^(2*Gparam(1) + BWoct^2)); %sigma of log-normal distribution

if pref<0
    BW = NaN; 
end

ffitII = interp1(log2(domI),ffitI,log2(domII),'spline')';  %High resolution interpolation of the fit for plotting purposes

[rfit pfit] = corrcoef(tcI',ffitI);
pfit = pfit(1,2);

[lco hco] = gethcolco(domII,ffitII,cothresh);

%Get bandpass factor
BP = (max(tc)-phi(tc(1)))/(max(tc)+phi(tc(1))); %Bandpass factor

%Compute center-of-mass of tuning curve
tcnorm = tc - min(tc);
tcnorm = tcnorm/sum(tcnorm);
CoM = 2.^sum(tcnorm.*log2(dom)); %center of mass

%Save parameters to structure
params.hco = hco; %high pass cutoff
params.BW = BW;  %linear bandwidth. This only makes sense if t ?
params.BWoct = BWoct;  %octave bandwidth.
params.BP = (max(tc)-phi(tc(1)))/(max(tc)+phi(tc(1))); %Bandpass factor
params.CoM = CoM;  %Center of mass.
params.pref = pref; %location of the peak, based on Gaussian fit.
params.raw = GparamRaw;

%Use the repeats to assess fit quality
if ~isempty(varargin)
    kern_rep = varargin{1}; %response matrix of log parameter vs. repeat
    ffit = exp(-(log2(dom)-params.raw(1)).^2/(2*params.raw(2).^2));
    %ffit = params.raw(3)*ffit + params.raw(4);
    ffit = params.raw(3)*ffit;
    ffitMat = ffit(:)*ones(1,getnorepeats(1));
    Repvaracc = (var(kern_rep(:))-var(kern_rep(:)-ffitMat(:)))/var(kern_rep(:));
    [rdum pfit] = corrcoef(kern_rep(:),ffitMat(:));
    repFitQuality.p = pfit(1,2);
    repFitQuality.r = rdum(1,2);    
    repFitQuality.varacc = Repvaracc;
else
    repFitQuality = [];
end

function [param ffit varacc] = Gfit(domain,f)

%%%%Normalize%%%%%%%
mi = min(f);
f = f-mi;
E = max(f);
f = f/E;

%% Compute initial guess %%

[dum,idx] = max(f);
muGuess = domain(idx);
sigGuess = range(domain)/4;
%sig = length(f)/8;

baseGuess = prctile(f,.1);
AmpGuess = max(f)-baseGuess;

%G = [idx length(f)/6 max(RF)-min(RF) min(RF)];
x0 = [muGuess sigGuess AmpGuess baseGuess]; %initial guess for paramter search
%G = [length(f)/2 1 0 0];

%% Search %%%
global dom RF
dom = domain;
RF = f;

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);


[param,dum] = fminsearch('gaussfitter_handle2',x0);


%%%% Compute the fit %%%%%%%

ffit = exp(-(domain-param(1)).^2/(2*param(2).^2));
ffit = param(3)*ffit + param(4);

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

%% Undo normalization

ffit = ffit*E + mi;
param(3) = param(3)*E;
param(4) = param(4)*E+mi;


function [param ffit varacc] = Gfit_base0(domain,f)

%%%%Normalize%%%%%%%
E = max(f);
f = f/E;

%% Compute initial guess %%

[dum,idx] = max(f);
muGuess = domain(idx);
sigGuess = range(domain)/4;
%sig = length(f)/8;

AmpGuess = max(f);

%G = [idx length(f)/6 max(RF)-min(RF) min(RF)];
x0 = [muGuess sigGuess AmpGuess]; %initial guess for paramter search
%G = [length(f)/2 1 0 0];

%% Search %%%
global dom RF
dom = domain;
RF = f;

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);


[param,dum] = fminsearch('gaussfitter_handle2',x0);


%%%% Compute the fit %%%%%%%

ffit = exp(-(domain-param(1)).^2/(2*param(2).^2));
ffit = param(3)*ffit;

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

%% Undo normalization

ffit = ffit*E;
param(3) = param(3)*E;



