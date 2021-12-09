function [param ffit MSE] = Gaussfit2Dpolar(f)

global RF

%%%search%%%
RF = f;
param = gaussfitter2Dpolar
%%%%%%%%%%%

cols = length(f(1,:));
rows = length(f(:,1));
[xx yy] = meshgrid(1:cols,1:rows);
xx = xx-ceil(cols/2);
yy = yy-ceil(rows/2);
yy = flipud(yy);



%Rotate the domain so that its continuous around the peak:
tt = atan2(yy,xx)*180/pi;
tt = tt + (90-param(1));
tt = angle(exp(1i*tt*pi/180))*180/pi;
id = find(sign(tt) == -1);
tt(id) = tt(id) + 180;
tt = tt-90;

rr = sqrt(yy.^2 + xx.^2);


d1 = tt.^2;  %ori domain
d2 = (rr-param(3)).^2;  %sf domain
ffit = param(6)*exp(-d2./(2*param(4)^2)).*exp(-d1./(2*param(2)^2))  +  param(5);

MSE = mean((ffit(:)-f(:)).*(ffit(:)-f(:)));


