function F1 = getF1(tResp)

global ACQinfo cellS Analyzer

Flim = getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning
Flim(2) = Flim(1) + 1000*getparam('stim_time');

tdom = (0:length(tResp)-1)*ACQinfo.msPerLine*ACQinfo.linesPerFrame;
[dum Flim(1)] = min(abs(tdom-Flim(1)));
[dum Flim(2)] = min(abs(tdom-Flim(2)));


try
    STIMfrate = Analyzer.framerate; %ms/period
catch
    STIMfrate = 60;
    'Frame rate does not exist in Analyzer. Setting to 60Hz'
end

T = getParamVal('t_period')/STIMfrate*1000; %ms/cycle

tdom = tdom - tdom(Flim(1));  %Need phase to start at zero

maxCycle = max(unique(floor(tdom(1:Flim(2))/T)));
maxID = find(tdom/T>=maxCycle);
maxID = maxID(1);
Flim(2) = maxID;

tdom = tdom(Flim(1):Flim(2));


phi = 2*pi*tdom/T;

%%

tResp = processTcourse(tResp,1,2)'; %smooth and subtract low-order polynomial
tResp = tResp(Flim(1):Flim(2));

dum = tResp(:)';
F1 = 2*dum*exp(1i*phi(:))/length(phi);

%fit = abs(F1)*cos(phi-angle(F1));
%varacc = (var(tResp(:))-var(tResp(:)-fit(:)))/var(tResp(:));


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

