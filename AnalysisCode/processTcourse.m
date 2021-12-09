function y = processTcourse(y,hh,polyorder,sp)

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

%Linear Bandpass Filter from GUI
if ~isempty(hh)
    
    dL = length(y)-length(hh);
    if dL<0
        y = [y zeros(1,-dL)];
    elseif dL>0
        y = y(1:end-dL);
    end
    
    y = ifft(fft(y).*hh);
end
%[y noise] = wiener2(y, [1 round(300/sp)]);
%y = zscore(y);

%figure,plot(abs(fft(y)))

%Get rid of any peaks around the breathing rate
fdom = linspace(0,1/(sp/1000),length(y)+1);
fdom = fdom(1:end-1);
atten = [.4 .15 .1 .15 .4]; %Notch filter
idbreath = find(fdom<2.5 & fdom>.25); %Look in this band for peaks
Np = 0;
Npeaks = 2;
for i = 1:Npeaks

    yf = fft(y);
    [ma id] = max(abs(yf(idbreath)));
    if ma > 3*std(abs(yf(idbreath))) + median(abs(yf(idbreath)));
        idpeak = id+idbreath(1)-1;
        yf(idpeak-2:idpeak+2) = yf(idpeak-2:idpeak+2).*atten;

        yf = fliplr(yf);
        [dum id] = max(abs(yf(idbreath)));
        idpeak = id+idbreath(1)-1;
        yf(idpeak-2:idpeak+2) = yf(idpeak-2:idpeak+2).*atten;
        yf = fliplr(yf);
        y = real(ifft(yf));
        Np = Np+1;
    else
        break
    end

end

%hold on,
%plot(abs(fft(y)),'r')
%title(num2str(Np))

%Nonlinear highpass filter to subtract estimated baseline

h = fspecial('gaussian', [1 length(y)], 10); %Use stats from heavily smoothed version
ydum = ifft(fft(y).*abs(fft(h)));
y = y - ordfilt2(ydum, 20, ones(1,600));

id = find(y<prctile(y,50) & y>prctile(y,0));
%y = y-median(y); 
%y = y/std(y(id));
%y = y/sqrt(mean(y(id).^2));

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