function [param ffit varaccount] = Gaborfit2Drot_nophase(f,varargin)

global RF errall
errall = [];

%f = f-nanmedian(f(:));
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
x0 = [G.ysig/5 G.xsig/5 G.ori G.sf];

%% Search

options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',1e-8,'TolX',.000004);  %this really helps for this function
param = fminsearch('Gaborfitter2Drot_handle_nophase',x0,options);

%% Build the fitted function

imWy = length(f(:,1));
imWx = length(f(1,:));
domx = -floor(imWx/2):ceil(imWx/2)-1;
domy = ceil(imWy/2)-1:-1:-floor(imWy/2);

[x y] = meshgrid(domx,domy);

xp = x*cos(param(3)*pi/180) + y*sin(param(3)*pi/180);
yp = y*cos(param(3)*pi/180) - x*sin(param(3)*pi/180); 

ffitD1 = exp(-(yp).^2/(2*param(1).^2));
ffitD2 = exp(-(xp).^2/(2*param(2).^2));

Env = ffitD1.*ffitD2;

dim = size(xp);

ffit = Env.*cos(xp/dim(1)*2*pi*param(4));
ffit90 = Env.*sin(xp/dim(1)*2*pi*param(4));

%ffit = ffit*param(5);

fw = abs(fft2(f));
%ffitw = abs(fft2(ffit));
ffitw = sqrt(abs(fft2(ffit)).^2 + abs(fft2(ffit90)).^2);

varaccount = (var(fw(:))-var(fw(:)-ffitw(:)))/var(fw(:));



%Put the means back on the original axis


%Convert units to degrees
pixperim = size(ffit,1);
param(1:2) = param(1:2)/pixperdeg;
param(4) = param(4)/pixperim*pixperdeg; %Convert cycles/image to cyc/degree


%figure,plot(errall)