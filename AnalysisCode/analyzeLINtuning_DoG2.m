function [params ffitII domI domII pfit varaccfit repFitQuality] = analyzeLINtuning_DoG2(tc,dom,cothresh,varargin)

%Got rid of all bandpass outputs and only output the lco and hco values.
%It was creating too much confusion 

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

InterpNplot = 20; %for domII

domI = logspace(log10(dom(1)),log10(dom(end)),length(dom)*InterpNfit-1);
domII = logspace(log10(dom(1)),log10(dom(end)),length(dom)*InterpNplot-1); %strictly for plotting a smooth curve

tcI = interp1(dom,tc,domI,'spline')'; %Interpolate tuning curve before fitting

[Gparam ffitI varaccfit] = Gfit(domI,tcI');
%[Gparam ffitI varaccfit] = Gfit_base0(domI,tcI');
if isnan(ffitI(1))
    ffitI = zeros(size(ffitI));
end
GparamRaw = Gparam;

%hco = Gparam(2) + Gparam(1);


ffitII = interp1(domI,ffitI,domII,'spline')';  %High resolution interpolation of the fit for plotting purposes

[rfit pfit] = corrcoef(tcI',ffitI);
pfit = pfit(1,2);

[dum pref] = max(ffitII);
pref = domII(pref);


[lco hco] = gethcolco(domII,ffitII,cothresh);
if isnan(lco)
    lco = min(dom);
end


%Compute center-of-mass of tuning curve
tcnorm = tc - min(tc);
tcnorm = tcnorm/sum(tcnorm);
CoM = sum(tcnorm.*dom); %center of mass

%Save parameters to structure
params.hco = hco; %high pass cutoff
params.lco = lco;
params.BP = (max(tc)-phi(tc(1)))/(max(tc)+phi(tc(1))); %Bandpass factor
params.CoM = CoM;  %Center of mass.
params.pref = pref; %location of the peak, based on Gaussian fit.
params.raw = GparamRaw;


%Use the repeats to assess fit quality
if ~isempty(varargin)
    kern_rep = varargin{1}; %response matrix of lin parameter vs. repeat
    
    
    ffit1 = params.raw(2)*exp(-(dom-params.raw(5)).^2./(2*params.raw(1)^2));
    ffit2 = params.raw(4)*exp(-dom.^2./(2*params.raw(3)^2));
    ffit = ffit1-ffit2;
    
%     ffit = exp(-((dom)-params.raw(1)).^2/(2*params.raw(2).^2));
%     ffit = params.raw(3)*ffit + params.raw(4);
    
    
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

E = norm(f);
f = f/E;


%% Search %%%
global logdom RF
logdom = domain; %not really a log domain
RF = f;

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);

param = DoGfitter2;

%[param,dum] = fminsearch('DoGfitter_handle',x0);
%[param,dum] = fmincon('SoGfitter_handle2',x0,[],[],[],[],[.25 0 -inf -.1 0],[inf 4 4 .2 4]);


%%%% Compute the fit %%%%%%%

%figure, plot(exp(-(domain-param(1)).^2/(2*param(2).^2)))
%shold on, plot(param(4)*exp(-(domain).^2/(2*param(2).^2)))

sx1 = param(1);
A1 = param(2);

sx2 = param(3);
A2 = param(4);

mu = param(5);

%base = param(5);

ffit1 = param(2)*exp(-(logdom-param(5)).^2./(2*param(1)^2));
ffit2 = param(4)*exp(-logdom.^2./(2*param(3)^2));
ffit = ffit1-ffit2;

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

%% Undo normalization

ffit = ffit*E;
param(3) = param(3)*E;
%param(4) = param(4)*E+mi;
