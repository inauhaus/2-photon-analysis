function y = processTcourse_NaN(y,hhi,hlo,polyorder,sp)

%Called from Ggetrandposkernel2 and Ggetrevcorrkernel2


%mu = mean(y);
%First subtract a low-order polynomial fit:
id = find(y<nanmedian(y)-2*nanstd(y));
ytrunc = y;
ytrunc(id) = NaN; %Remove outliers (large transient responses) before fitting
yfit = polyfitter(ytrunc,polyorder);
y = y-yfit';

%Linear Bandpass Filter from GUI
ydum = y;

if ~isempty(hlo)
    
    %hhMat = toeplitz([hh(1) fliplr(hh(2:end))], hh);
    
    counter = zeros(size(hlo));
    counter(find(~isnan(y))) = 1;
    
    y(find(isnan(y))) = 0;
    yh = ifft(fft(y).*hlo);
    ch = ifft(fft(counter).*hlo);
    y = yh./ch;
    y(find(1-counter)) = NaN;

end

if ~isempty(hhi)
    
    %hhMat = toeplitz([hh(1) fliplr(hh(2:end))], hh);
    
    counter = zeros(size(hhi));
    counter(find(~isnan(y))) = 1;
    
    y(find(isnan(y))) = 0;
    yh = ifft(fft(y).*hhi);
    ch = ifft(fft(counter).*hhi);
    yh = yh./ch;
    
    y = y-yh;
    y(find(1-counter)) = NaN;
end




%[y noise] = wiener2(y, [1 round(300/sp)]);
%y = zscore(y);

%figure,plot(abs(fft(y)))

%Get rid of any peaks around the breathing rate
% fdom = linspace(0,1/(sp/1000),length(y)+1);
% fdom = fdom(1:end-1);
% atten = [.4 .15 .1 .15 .4]; %Notch filter
% idbreath = find(fdom<2.5 & fdom>.25); %Look in this band for peaks
% Np = 0;
% Npeaks = 2;
% for i = 1:Npeaks
% 
%     yf = fft(y);
%     [ma id] = max(abs(yf(idbreath)));
%     if ma > 3*std(abs(yf(idbreath))) + median(abs(yf(idbreath)));
%         idpeak = id+idbreath(1)-1;
%         yf(idpeak-2:idpeak+2) = yf(idpeak-2:idpeak+2).*atten;
% 
%         yf = fliplr(yf);
%         [dum id] = max(abs(yf(idbreath)));
%         idpeak = id+idbreath(1)-1;
%         yf(idpeak-2:idpeak+2) = yf(idpeak-2:idpeak+2).*atten;
%         yf = fliplr(yf);
%         y = real(ifft(yf));
%         Np = Np+1;
%     else
%         break
%     end
% 
% end




function xfit = polyfitter(x,order)

%This one accounts for NaNs

dom = (0:length(x)-1)';

H = ones(length(dom),order+1);  %last column for DC
for i = 1:order
    H(:,i) = dom.^i;
end
Hfull = H;

id = find(isnan(x));
H(id,:) = [];
x(id) = [];


p =  inv(H'*H)*H'*x';


xfit = Hfull*p(:);