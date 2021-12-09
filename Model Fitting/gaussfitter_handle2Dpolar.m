function err = gaussfitter_handle2Dpolar(param)

global RF;

dim = size(RF);

xc1 = param(1);
sx1 = param(2);
xc2 = param(3);
sx2 = param(4);

base = param(5);
A = param(6);

%% Make continuous at the peak. Center at 0 deg
[xx yy] = meshgrid(1:dim(2),1:dim(1));
xx = xx-ceil(dim(2)/2);
yy = yy-ceil(dim(1)/2);
yy = flipud(yy);

tt = atan2(yy,xx)*180/pi;

%Rotate the domain so that its continuous around the peak:
tt = tt + (90-xc1);
tt = angle(exp(1i*tt*pi/180))*180/pi;

id = find(sign(tt) == -1);
tt(id) = tt(id) + 180;
tt = tt-90;
%%

rr = sqrt(yy.^2 + xx.^2);

%tt = tt-90;
rr = rr-xc2;

img = A*exp(-(rr.^2)./(2*sx2^2)).*exp(-(tt.^2)./(2*sx1^2))  +  base;


err = sum((img(:)-RF(:)).^2);