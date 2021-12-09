function [x,f,g] = gaussfitter2Dpolar

global RF

x0 = gaussfitguess2Dpolar;
dim = size(RF);

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);
[x,f] = fminsearch('gaussfitter_handle2Dpolar',x0);

xc1 = x(1);
sig1 = x(2);

xc2 = x(3);
sig2 = x(4);

B = x(5);
A = x(6);

[xx yy] = meshgrid(1:dim(2),1:dim(1));
xx = xx-ceil(dim(2)/2);
yy = yy-ceil(dim(1)/2);
yy = flipud(yy);

%Rotate the domain so that its continuous around the peak:
tt = atan2(yy,xx)*180/pi;
tt = tt + (90-xc1);
tt = angle(exp(1i*tt*pi/180))*180/pi;
id = find(sign(tt) == -1);
tt(id) = tt(id) + 180;
tt = tt-90; %Centered at zero

%tt(:,1:floor(dim(2)/2)) = 180-fliplr(tt(:,ceil(dim(2)/2)+1:end));
rr = sqrt(yy.^2 + xx.^2);

d1 = tt.^2; %already centered above
d2 = (rr-xc2).^2;
g = A*exp(-d1/(2*sig1^2)).*exp(-d2/(2*sig2^2))+B;


