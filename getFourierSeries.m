function [Rfit F1 E] = getFourierSeries(R,Nh,phi)    

%R - raw signal
%phi - phase of F1 (same length as R)
%Nh - number of harmonics to use in reconstruction

R = R-mean(R);

%polyorder = 2; LPsig = 3;
%R = processTcourse(R,LPsig,polyorder);

harmvec = 1:Nh;
Rfit = zeros(size(R));

for h = 1:length(harmvec) %loop each harmonic
    dum = R.*hann(length(phi))';
    F1x = 2*dum*exp(1i*harmvec(h)*phi(:))/length(phi);
    Rfit = Rfit + abs(F1x)*cos(harmvec(h)*phi-angle(F1x));
end

F1 = R*exp(1i*phi(:));
F1 = 2*F1/sum(hann(length(phi)));

%SNR varies between 



function y = processTcourse(y,LpSig,polyorder)

%Called from Ggetrandposkernel2 and Ggetrevcorrkernel2

id = find(isnan(y));
y(id) = nanmean(y);

%mu = mean(y);
%First subtract a low-order polynomial fit:
id = find(y>median(y)+2*std(y));
ytrunc = y;
ytrunc(id) = median(y); %Remove outliers (large transient responses) before fitting
yfit = polyfitter(ytrunc,polyorder);
y = y-yfit';

%Linear smoothing
if ~isempty(LpSig)
    hh = abs(fft(fspecial('gaussian',size(y),LpSig)));
    y = ifft(fft(y).*hh);
end

%Nonlinear highpass filter to subtract estimated baseline

% h = fspecial('gaussian', [1 length(y)], 10); %Use stats from heavily smoothed version
% ydum = ifft(fft(y).*abs(fft(h)));
% y = y - ordfilt2(ydum, 20, ones(1,600));

% id = find(y<prctile(y,50) & y>prctile(y,0));
% y = y-median(y); 
% y = y/std(y(id));

%y = zscore(y);
%y = y.*sign(y);

function xfit = polyfitter(x,order)

dom = (0:length(x)-1)';

H = ones(length(dom),order+1);  %last column for DC
for i = 1:order
    H(:,i) = dom.^i;
end

p = inv(H'*H)*H'*x';
xfit = H*p;
