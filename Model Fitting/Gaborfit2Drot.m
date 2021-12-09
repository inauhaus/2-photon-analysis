function [param ffit varaccount] = Gaborfit2Drot(f,varargin)

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
x0 = [G.ymu G.xmu G.ysig G.xsig/2 G.ori G.sf G.phase];

%% Search

options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',1e-8,'TolX',.000004);  %this really helps for this function
param = fminsearch('Gaborfitter2Drot_handle',x0,options);

%% Build the fitted function

imWy = length(f(:,1));
imWx = length(f(1,:));
domx = -floor(imWx/2):ceil(imWx/2)-1;
domy = ceil(imWy/2)-1:-1:-floor(imWy/2);

[x y] = meshgrid(domx,domy);

xp = x*cos(param(5)*pi/180) + y*sin(param(5)*pi/180);
yp = y*cos(param(5)*pi/180) - x*sin(param(5)*pi/180); 


% ffitD1 = (cos(2*pi*(yp-param(1))/(4*param(3)))+1)/2;
% ffitD2 = (cos(2*pi*(xp-param(2))/(4*param(4)))+1)/2;
% ffitD1(find(abs(yp-param(1)) > param(3)*2)) = 0;
% ffitD2(find(abs(xp-param(2)) > param(4)*2)) = 0;

ffitD1 = exp(-(yp-param(1)).^2/(2*param(3).^2));
ffitD2 = exp(-(xp-param(2)).^2/(2*param(4).^2));
Env = ffitD1.*ffitD2;

dim = size(xp);
carr = cos(xp/dim(1)*2*pi*param(6) + param(7)*pi/180);

ffit = Env.*carr;

varaccount = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

%Put the means back on the original axis

p2 = param(2)*cos(-param(5)*pi/180) + param(1)*sin(-param(5)*pi/180);
p1 = param(1)*cos(-param(5)*pi/180) - param(2)*sin(-param(5)*pi/180);
param(1) = p1;
param(2) = p2;


%Convert units to degrees
pixperim = size(ffit,1);
param(1:4) = param(1:4)/pixperdeg;
param(6) = param(6)/pixperim*pixperdeg; %Convert cycles/image to cyc/degree


%figure,plot(errall)