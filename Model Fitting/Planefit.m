function [param ffit varacc] = Planefit(domx,domy,f,varargin)

%Ian Nauhaus

%This fits a plane to a circular variable in the axis domain

global RF xpts ypts

RF = f; xpts = domx; ypts = domy;

%% make the guess
if isempty(varargin)
    
    %very important to get an initial estimate of the variance
    angpref = angle(nansum(exp(1i*f*2*pi/180)));
    angpref = angpref + pi*(1-sign(angpref+eps));  %0 to 2pi
    angpref = angpref*90/pi;               %0 to 180
    
    fdum = angle(exp(1i*(f-angpref)*2*pi/180))*90/pi;
    
    xslopeguess = nanstd(fdum(:))/nanstd(xpts);
    yslopeguess = nanstd(fdum(:))/nanstd(ypts);
    baseguess = nanmean(f(:));
else
    xslopeguess = varargin{1}(1);
    yslopeguess = varargin{1}(2);
    baseguess = varargin{1}(3);
end

x0 = [xslopeguess yslopeguess baseguess];


%% search
options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.000004,'TolX',.000004);
%[param,f] = fminsearch('planefitter_handle',x0,options);
x0
[param,f] = fminsearch('planefitter_handle',x0);

%f = angle(exp(1i*f*pi/180))*180/pi;



%%%%%%%%%%%
%%

ffit = param(1)*domx + param(2)*domy + param(3);

ffit = angle(exp(1i*ffit*pi/180))*180/pi; %wrap it
id = find(ffit(:)<0);
ffit(id) = ffit(id)+360;
ffit = ffit/2;

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

