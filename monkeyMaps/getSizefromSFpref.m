function xsig = getSizefromSFpref(sf2sigModel,sfpref)

global alpha

switch sf2sigModel
    
    case 'Albrecht'
        
        bwmesh = -sfpref/6 + 2;  %Albrecht et al BW vs. SF (octaves)
        xsig = 1./(pi*2.^bwmesh)/2;  %convert to degrees
        
    case 'classic'
        
        
        xsig = (1./(sfpref))/alpha;  %classic model: sig = 1/(4sf)
        
    case 'exp decay'
        
        
%         logxsig = -.12*(-sfpref).^2 + 0.54*(-sfpref) + 0.189;
%         xsig = (2.^logxsig)/4;
%         
        param = [1.76 0.57 0.197];
        
        xsig = param(1)*exp(-sfpref*param(2)) + param(3);
        xsig = xsig/4;  %model was fit using 4sig.
        
    case 'NL sat'
        
       % param = [-.8422 -.5086];
                
%         a = log2(1./(2*sfpref));
%         %b = log2(2*xsig);        
%         b = a-param(2);
%         b(find(b>param(1))) = param(1);  %saturation level        
%         xsig = (2.^b)/2;
%         %sfpref = 1./(2.^a)/2; xsig = (2.^b)/2
        
       
        param = [1.5238 0.3417];
        
        a = (1./(alpha*sfpref));
        b = a*param(1);
        %b = log2(xsig);
          
        b(find(b>param(2))) = param(2);  %saturation level
        
        sig = .1;
        b = smoothFit(a(:),b(:),param,sig);

        b = reshape(b,size(a));
        
        xsig = (b);
        


end


function ffit = smoothFit(xxraw,ffit,param,sig)

[xx id idu] = unique(xxraw);
ffit = ffit(id);

%Make padded and interpolated domain
xpad = range(xx);
xI = linspace(min(xx)-xpad,max(xx)+xpad,10000); %This needs a bunch of interpolation for some reason
%xI = logspace(log10(min(xx)-xpad),log10(max(xx)+xpad),5000);
dx = xI(2)-xI(1);

%Make padded and interpolated fit
pad1 = xI(1)*param(1);
pad2 = xI(end)*param(1);
if pad2>param(2)
    pad2 = param(2);
end

% [xxsort id] = sort(xx);
% ffit = ffit(id);

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
ffit = zeros(1,length(xxraw));
for i = 1:length(xx)
    
    [dum id] = min(abs(xx(i)-xI));
    ffitdum(i) = ffitI(id);
    
    %id = find(xx(i) == xxraw);
    %ffit(id) = ffitdum;
    
end

ffit = ffitdum(idu); %This is way faster than doing it in the 'for' loop


