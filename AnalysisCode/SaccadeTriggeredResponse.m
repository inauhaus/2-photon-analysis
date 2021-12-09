function SaccadeTriggeredResponse

global cellS ACQinfo TimingInfo

rThresh = 0.95;

nT = length(TimingInfo.EyeLink_xy);

fp_Eye = median(diff(TimingInfo.VBLTimestamp{1}))*1000;
taudom = -3000:fp_Eye:3000;



fp_Brain = ACQinfo.msPerLine*ACQinfo.linesPerFrame; %ms
figure
for t = 1:nT
    
    [c r] = getcondrep(t);
    
    xpos = TimingInfo.EyeLink_xy{t}(1,:);
    ypos = TimingInfo.EyeLink_xy{t}(2,:);
    
    %id lost track times
    dx = [0 abs(diff(xpos))];
    dy = [0 abs(diff(ypos))];
    dr = sqrt(dx.^2 + dy.^2);  
    idNoTrack = find(dr>5000 | dr == 0);   
    %%%%%%%%%%%%%%%
    
    xpos(idNoTrack) = NaN;
    ypos(idNoTrack) = NaN;
    
    medfiltN = 21;
    yposF = medfilt1(ypos,medfiltN);
    xposF = medfilt1(xpos,medfiltN);
    
%     figure,
%     plot(((ypos)))
%     hold on
%     plot(((yposF)))
   
    dx = [0 (diff(xposF))];
    dy = [0 (diff(yposF))];
    
    ythresh = nanstd(dx)*3;
    xthresh = nanstd(dy)*3;
    
    dr = sqrt(dx.^2 + dy.^2);
    
    idnotSaccade = find(dx<xthresh & dx>-xthresh);
    dx = abs(dx);
    dx(idnotSaccade) = 0;
    
    idnotSaccade = find(dy<ythresh & dy>-ythresh);
    dy = abs(dy);
    dy(idnotSaccade) = 0;
    
    drCleaned = sqrt(dx.^2 + dy.^2);
    
    thrsh = 50;
    drCleaned = (sign(drCleaned-thrsh)+1)/2;

%     if t == 21
%        
%         'hi'
%         
%     end
    
    %drCleaned = sign(drCleaned);
    
%     figure,plot(dr)
%     hold on 
%     plot(sign(drCleaned)*200)
    
   
    tdomEye = 0:fp_Eye:fp_Eye*(length(xpos)-1);

    %dr = dr-mean(dr);
    %dr = dr/std(dr);
    
    Nbrain = size(cellS.cellMat{c},2); %number of samples in trial
    tdomBrain = 0:fp_Brain:fp_Brain*(Nbrain-1);
    
    
    if length(size(cellS.cellMat{c})) > 2
        Rpop = squeeze(cellS.cellMat{c}(:,:,r))';
    else 
        Rpop = squeeze(cellS.cellMat{c})';
    end
    
    skewthresh = 0;
    
    idSkew = find(skewness(Rpop)>skewthresh);
    %Rpop = Rpop(:,idSkew);
    
    hh = makeTemporalfilter(Nbrain);
    
    for p = 1:length(Rpop(1,:))
        %Rpop(:,p) = processTcourse(Rpop(:,p)',hh,1,fp_Brain); 
        %Rpop(:,p) = zscore(Rpop(:,p));
        
        Rpop(:,p) = Rpop(:,p) - nanmean(Rpop(:,p));
        %Rpop(:,p) = Rpop(:,p)/nanstd(Rpop(:,p));
    end
    
    %Rpop = Rpop.*(sign(Rpop)+1)/2;
    
    Rpop = nanmean(Rpop,2)';
    %Rpop = Rpop(:,15);
    Rpop = Rpop/nanstd(Rpop);
   % Rpop = Rpop - nanmean(Rpop);
    %Rpop = Rpop/nanstd(Rpop);
    %Rpop = zscore(Rpop);
    
    figure(20),
    
    subplot(ceil(sqrt(nT)),ceil(sqrt(nT)),t)
    plot(Rpop)
    
    
    
    motionVec = cellS.motionInfo.rTemp{t}';
    
    LD = length(Rpop) - length(motionVec);
    if LD>0
       motionVec = [motionVec(:)' zeros(1,LD)];
    end

    RpopI = interp1(tdomBrain,Rpop,tdomEye);
    %motionVecI = interp1(tdomBrain,motionVec,tdomEye);
    
%     
%     rThresh = prctile(motionVecI,10);
%     idMotion = find(motionVecI<rThresh);
%     RpopI(idMotion) = NaN;
%     RpopI = RpopI-nanmean(RpopI);
    
    drCleaned = drCleaned-nanmean(drCleaned);
    
    XC(t,:) = myXcorr(drCleaned,RpopI,round(taudom/fp_Eye));
    
%     x = randn(size(RpopI));
%     y = circshift(x,[0 20]);
%     XC(t,:) = myXcorr(y,x,round(taudom/fp_Eye));
    
    figure(22)
    subplot(ceil(sqrt(nT)),ceil(sqrt(nT)),t)
    
    plot(taudom,XC(t,:))
    %plot(RpopI)
    
    
end

figure,plot(taudom,nanmean(XC)), xlabel('ms after saccade')
%%

function Rxy = myXcorr(x,y,shiftdom)

for i = 1:length(shiftdom)
       
    dn = shiftdom(i);
    if dn>=0
        Rxy(i) = nansum(x(1:end-dn).*y(dn+1:end));
    else        
        dn = -dn;
        Rxy(i) = nansum(x(dn+1:end).*y(1:end-dn));
    end
end
    
    
function hh = makeTemporalfilter(N)

global G_RChandles ACQinfo

togstateHP = get(G_RChandles.HPflag,'Value');
togstateLP = get(G_RChandles.LPflag,'Value');

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 

% if ~isempty(varargin)
%     N = getTrialLength(varargin{1});
% else
%     N = getTrialLength(1); %Just use the first trial
% end

if togstateHP == 1
    Hwidth = str2double(get(G_RChandles.Hwidth,'string'));
    Hwidth = round(Hwidth/acqPeriod); %convert to samples
    ind = get(G_RChandles.HPWind,'value');

    switch ind
        case 1 %Gaussian
            dom = (1:N)-round(N/2);
            H = exp(-dom.^2/(2*Hwidth^2));
            H = -H/sum(H);
            H(round(N/2)) = 1+H(round(N/2));
        case 2 %Hann
            H = zeros(1,N);
            Hd = hann(Hwidth);
            Hd = -Hd./sum(Hd(:));
            Hd(round(Hwidth/2)) = 1+Hd(round(Hwidth/2));
            H(1:Hwidth) = Hd;
        case 3 %Disc
            H = zeros(1,N);
            Hd = -ones(1,Hwidth)/Hwidth;
            Hd(round(Hwidth/2)) = 1+Hd(round(Hwidth/2));
            H(1:Hwidth) = Hd;
    end
    if togstateLP == 0
        hh = abs(fft(H));   %Eliminate phase information
    end
end

if togstateLP == 1
    Lwidth = str2double(get(G_RChandles.Lwidth,'string'));
    ind = get(G_RChandles.LPWind,'value');

    switch ind
        case 1            
            Lwidth = Lwidth/acqPeriod; %convert to samples (can be decimal for Gaussian)
            dom = (1:N)-round(N/2);
            L = exp(-dom.^2/(2*Lwidth^2));
            L = L/sum(L);
        case 2            
            Lwidth = round(Lwidth/acqPeriod); %convert to samples
            L = zeros(1,N);
            Ld = hann(Lwidth);
            Ld = Ld./sum(Ld(:));
            L(1:Lwidth) = Ld;            
        case 3            
            Lwidth = round(Lwidth/acqPeriod); %convert to samples
            L = zeros(1,N);
            Ld = ones(1,Lwidth)/Lwidth;
            L(1:Lwidth) = Ld;            
    end
    if togstateHP == 0
        hh = abs(fft(L(:)))';   %Eliminate phase information
    else
        hh = abs(fft(L(:)).*fft(H(:)))';   %Take mag because phase gives a slight shift.
    end
end

if ~or(togstateLP,togstateHP)
    hh = [];
end


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

% h = fspecial('gaussian', [1 length(y)], 10); %Use stats from heavily smoothed version
% ydum = ifft(fft(y).*abs(fft(h)));
% y = y - ordfilt2(ydum, 20, ones(1,600));
% 
% id = find(y<prctile(y,50) & y>prctile(y,0));
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




