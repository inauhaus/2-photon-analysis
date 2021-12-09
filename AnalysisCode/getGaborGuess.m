function G = getGaborGuess(im,pixperdeg)

imW = length(im(:,1));
domx = -floor(imW/2):ceil(imW/2)-1;
domy = ceil(imW/2)-1:-1:-floor(imW/2);

% ori = 0;
% xp = x*cos(ori*pi/180) + y*sin(ori*pi/180);  %clockwise
% yp = -x*sin(ori*pi/180) + y*cos(ori*pi/180);
% phase = 90*pi/180;
% im = cos(xp*2*pi/5 - phase);


%% Get carrier params

%figure,imagesc(im), axis image

[dumy dumx imC] = getGaborCenter(im,5); %This also multiplies it by a hann window

rfw = fft2((imC-mean(imC(:))));
rfwMag = abs(rfw);
rfwMag = fftshift(fftshift(rfwMag,1),2);


%hh = fspecial('gaussian',size(rfw),.5);
%magsmooth = ifft2(fft2(rfwMag).*abs(fft2(hh)));
magsmooth = rfwMag;


dumSmooth = magsmooth(:,(size(magsmooth,2)/2+1):end);

[idy idx] = find(dumSmooth == max(dumSmooth(:)));
idy = idy(1); idx = idx(1);
pkx = domx(idx) - domx(1);
pky = domy(idy)+1;
origuess = atan2(pky,pkx)*180/pi;


sfguess = sqrt(pkx.^2+pky.^2);  %cycles per image

% regenerate the plot

[x y] = meshgrid(domx,domy);

xp = x*cos(origuess*pi/180) + y*sin(origuess*pi/180);  %clockwise
%yp = -x*sin(origuess*pi/180) + y*cos(origuess*pi/180);

dim = size(xp);
imcomplex = cos(xp/dim(2)*2*pi*sfguess) + 1i*sin(xp/dim(2)*2*pi*sfguess) ;

phaseguess = angle(imcomplex(:)'*im(:)) * 180/pi;

Carrierguess = cos(xp/dim(2)*2*pi*sfguess + phaseguess*pi/180);

%figure,imagesc(Carrierguess), axis image

%% Get envelope params

[x y] = meshgrid(domx,domy);

xp = x*cos(origuess*pi/180) + y*sin(origuess*pi/180);
yp = y*cos(origuess*pi/180) - x*sin(origuess*pi/180);

[yg xg] = getGaborCenter(im,5);

xg = domx(xg);
yg = domy(yg);


xgp = xg*cos(origuess*pi/180) + yg*sin(origuess*pi/180); %clockwise
ygp = -xg*sin(origuess*pi/180) + yg*cos(origuess*pi/180);

sigx = .2*pixperdeg;
sigy = .2*pixperdeg;
ffitD1 = exp(-(yp-ygp).^2/(2*sigy.^2));
ffitD2 = exp(-(xp-xgp).^2/(2*sigx.^2));
ffit = ffitD1.*ffitD2.*Carrierguess;

%figure,imagesc(ffit)
axis image

G.ymu = ygp;
G.xmu = xgp;
G.ysig = sigy;
G.xsig = sigx;
G.ori = origuess;
G.sf = sfguess;
G.phase = phaseguess;


