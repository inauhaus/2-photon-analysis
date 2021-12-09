function [param ffit varaccount] = GaborFdomainFit(f,varargin)

global RF errall
errall = [];

% f = f-nanmedian(f(:));
f = f/max(abs(f(:)));

if ~isempty(varargin{1})
   pixperdeg = varargin{1}; 
else
    pixperdeg = 34;
end


%%%search%%%
RF = f;
global RF 

%% Get initial guess

G = getGaborGuess(RF,pixperdeg);
x0 = [G.ymu G.xmu G.ysig G.xsig G.ori G.sf G.phase];
x0 = [] 
x0 = GaborFdomainGuess;

%% Search

options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',1e-8,'TolX',.000004);  %this really helps for this function
param = fminsearch('Gaborfitter2Drot_handle',x0,options);

%% Build the fitted function

imWy = length(f(:,1));
imWx = length(f(1,:));
domx = -floor(imWx/2):ceil(imWx/2)-1;
domy = ceil(imWy/2)-1:-1:-floor(imWy/2);

[x y] = meshgrid(domx,domy);

xp = abs(x*cos(param(4)*pi/180) + y*sin(param(4)*pi/180));
yp = y*cos(param(4)*pi/180) - x*sin(param(4)*pi/180); 

ffitD1 = exp(-yp.^2/(2*param(2).^2));
ffitD2 = exp(-(xp-param(1)).^2/(2*param(3).^2));  %param(1,2) is spatial freq preference and bandwidth;
Env = ffitD1.*ffitD2;

dim = size(xp);

ffit = Env.*carr;

varaccount = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

%Put the means back on the original axis

% p2 = param(1)*cos(-param(4)*pi/180) + 0*sin(-param(4)*pi/180);
% p1 = 0*cos(-param(4)*pi/180) - param(1)*sin(-param(4)*pi/180);
% param(1) = p2;  %x location
% param(4) = p1;  %y location

%Convert units to degrees
pixperim = size(ffit,1);
param(1:4) = param(1:4)/pixperdeg;
param(6) = param(6)/pixperim*pixperdeg; %Convert cycles/image to cyc/degree


%figure,plot(errall)