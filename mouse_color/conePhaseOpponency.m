function conePhaseOpponency

%getCellStats  

%%

global ACQinfo cellS Analyzer

Flim = getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning
Flim(2) = Flim(1) + 1000*getparam('stim_time');

tdom = (0:size(cellS.cellMat{1},2)-1)*ACQinfo.msPerLine*ACQinfo.linesPerFrame;
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
colordom = getdomain('theta');
sdir = find(colordom == 90);
mdir = find(colordom == 0);

k = 1;
Ncell = length(cellS.mukern);
vthresh = -20; %threshold for variance accounted for.
figure
for i = 2:Ncell 

   oritc = mean(cellS.mukern{i},1);
   
   [dum oribest] = max(oritc);
   %colortc = squeeze(cellS.mukern{i}(id,:,:));
   
   %oritc = squeeze(mean(kern{i},1));
   tcourse1 = squeeze(cellS.mukernTime{i}(mdir,oribest,:));
   tcourse2 = squeeze(cellS.mukernTime{i}(sdir,oribest,:));
   
   tcourse1 = tcourse1(Flim(1):Flim(2)); %truncate to stimulus window.
   tcourse2 = tcourse2(Flim(1):Flim(2));
   
   tcourse1 = processTcourse(tcourse1',1,2)'; %smooth and subtract low-order polynomial
   tcourse2 = processTcourse(tcourse2',1,2)';
   
   dum = tcourse1(:)';
   F1_1 = 2*dum*exp(1i*phi(:))/length(phi);   
   fit_1 = abs(F1_1)*cos(phi-angle(F1_1));
   
   dum = tcourse2(:)';
   F1_2 = 2*dum*exp(1i*phi(:))/length(phi);
   fit_2 = abs(F1_2)*cos(phi-angle(F1_2));
   
   varacc_1(i) = (var(tcourse1(:))-var(tcourse1(:)-fit_1(:)))/var(tcourse1(:));
   varacc_2(i) = (var(tcourse2(:))-var(tcourse2(:)-fit_2(:)))/var(tcourse2(:));
   
   phasediff(i) = abs(angle(exp(1i*(angle(F1_2)-angle(F1_1)))));
   
   if varacc_1(i)>vthresh & varacc_2(i)>vthresh
       
       subplot(8,9,k)
       
       plot(tdom/1000,tcourse1,'b'),
       hold on,
       %plot(tdom(1:3:end)/1000,fit_1(1:3:end),'.b') %overlay the fit
       
       plot(tdom/1000,tcourse2,'r')
       %hold on,
       %plot(tdom(1:3:end)/1000,fit_2(1:3:end),'.r')
       
       axis off
       title(num2str(round(phasediff(i)*180/pi)))
       k = k+1;
       
       ylim([min([tcourse1; tcourse2]) max([tcourse1; tcourse2])])
   end
   
end

id = find(varacc_1<vthresh | varacc_2<vthresh);
phasediff(id) = NaN;
figure,
hist(abs(phasediff*180/pi))
xlabel('phase differential')

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

