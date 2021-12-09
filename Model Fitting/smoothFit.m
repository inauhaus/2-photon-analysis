function ffit = smoothFit(xx,ffit,param,sig)

%Used by NLsoftsat and NLsat_softfitter_handle

%Make padded and interpolated domain
xpad = range(xx);
xI = linspace(min(xx)-xpad,max(xx)+xpad,100000);
dx = xI(2)-xI(1);

%Make padded and interpolated fit
pad1 = xI(1)*param(1);
pad2 = xI(end)*param(1);
if pad2>param(2)
    pad2 = param(2);
end
ffitpad = [pad1; ffit; pad2];
xxpad = [xI(1); xx; xI(end)];
ffitsharp = interp1(xxpad,ffitpad,xI);

%Smooth the fit
sig = sig/dx;
h = fspecial('gaussian',size(ffitsharp),sig);
ffitI = ifft(abs(fft(h)).*fft(ffitsharp));

%Remove the padding
id = find(xI>min(xx) & xI<(max(xx)));
xI = xI(id);
ffitI = ffitI(id);

%Look up the fit at each raw x location
clear ffit
for i = 1:length(xx)
    [dum id] = min(abs(xx(i)-xI));
    ffit(i) = ffitI(id);
end


