function kern = revCorrLS9(cellS,trialdom)

%9 is an extension of 5 which shows the prediction

global ACQinfo maskS Analyzer G_RChandles DM

hh = makeTemporalfilter;

DC = 0;
polyorder = 1;
blnkflag = 0;

cellMat = cellS.cellMat;
%synctimes = cellS.synctimes;
%frameid = cellS.frameid;
clear cellS

nID = getNeuronMask;  %get the index values for the neurons
masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
Ncell = length(nID);

ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];
% try
%     load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt],'frate')
% catch
%     load(['e:\2p_data\' Analyzer.M.anim '\log_files\' expt],'frate')
% end
frate = 60;
Tf = 1000/frate;  %Frame period in ms (frate obtained from log file) 

[domains seqs] = getSeqInfo(trialdom);

%%%%%%%%%%%%%%%%%%%

oridom = DM.oridom;
sfdom = DM.sfdom;

paramP = length(oridom)*length(sfdom);

tauP = round(166/acqPeriod)+1;

taudom = 0:acqPeriod:((tauP-1)*acqPeriod);
clear kern
covMat = 0;
xCorr = cell(1,Ncell);
covMat = cell(1,Ncell);
for i = 1:Ncell
    xCorr{i} = 0;
    covMat{i} = 0;
end

ARord = 1;
ARdel = 0;

%pdom = [4 7 12];
pdom = 1:18;
for pdum = 1:length(pdom)
    p = pdom(pdum);
    pID = nID(p);
    [idcelly idcellx] = find(masklabel == celldom(pID));
    CoM = [mean(idcelly) mean(idcellx)];  %center of mass
    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);
    tauConst(p) = NaN;
    
%     [param ffit fitdom varaccount] = Saturationfit2(yhat_all{p},y_all{p},[1.5 .7 0]);
%     [param2 ffit2 varaccount2 fitdom2] = Expfit2(yhat_all{p},y_all{p},[1.5 -.7 0]);
%     
%     if varaccount2>varaccount
%         param = param2;
%         ffit = ffit2;
%         fitdom = fitdom2;
%     end
    
%     figure, scatter(yhat_all{p},y_all{p},'.k')
%     hold on, plot(fitdom,ffit,'r')

    for trialid = 1:length(trialdom)

        T = trialdom(trialid);
        [cond rep] = getcondrep(T);
        
        hper = gethper(cond);         
        Tupdate = Tf*hper;   

        y = squeeze(cellMat{cond}(pID,:,rep)); 
        
        tdom = (0:length(y)-1)*acqPeriod;
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain
        %of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;                
        tdom_pix = tdom_pix + acqPeriod/2;
        
        %y = squeeze(mean(cellMat{cond}(:,:,rep),1));
              
        y = processTcourse(y,hh,polyorder,acqPeriod);
        
        
        %%%%
%         y = zscore(y);
%         acov = fftshift(xcov(y,'unbiased'));
%         noise = (acov(end) - acov(end-1))*length(acov); %estimate of noise variance; assuming its white, additive, and independent from everything else
%         N = noise*ones(size(y));
%         S = conj(fft(y)).*fft(y);
%         Ctrx = exp(-tdom/800);
%         H = fft(Ctrx);
%         G = (conj(H).*S)./(S.*abs(H).^2 + N);
%         y = ifft(G.*fft(y));
        %y = ifft(fft(y)./fft(Ctrx));
        %%%%
%         
%         if varaccount2>varaccount
%             y = -log((y-param(3))/param(1))/param(2);
%         else
%             y = -log(.01 + phi(1 - (y-param(3))/param(1)))/param(2);
%         end       
%         y = real(y);
        

        HHdum = zeros(length(y),length(oridom)*length(sfdom));
        pid = 1;
        for ori = 1:length(oridom)
            for sf = 1:length(sfdom)

                id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq == sfdom(sf));
                %id = find(seqs{T}.oriseq == oridom(ori));

                if ~isempty(id)

                    stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times (ms)
                    %stimes = synctimes{cond,rep}(id)*1000;

                    for i = 1:length(stimes)
                        [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                        HHdum(idx1,pid) = 1;
                    end
                end
                pid = pid+1;
            end
        end
        
        if blnkflag

            blnkvec = zeros(length(HHdum(:,1)),1);
            blnkid = find(isnan(seqs{T}.oriseq));
            stimes = synctimes{cond,rep}(blnkid)*1000;
            for i = 1:length(stimes)
                [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                blnkvec(idx1) = 1;
            end
            HHdum = [HHdum blnkvec];
            
        end
        
        %Truncate to timepoints with stimulation
        sid = sum(HHdum,2);
        id = find(sid);
        HHdum = HHdum(id(1):id(end),:);
        y = y(id(1):id(end));
        
        %Create stimulus matrix
        HH = zeros(length(HHdum(:,1)),length(HHdum(1,:))*tauP); %preallocate
        HHdum = [.0*randn(tauP-1,length(HHdum(1,:))); HHdum];  %Pad with zeros
        for z = tauP:length(HHdum(:,1))
            chunk = squeeze(HHdum(z-tauP+1:z,:))';
            chunk = chunk(:)';
            HH(z-tauP+1,:) = chunk;
        end
        %HH = zscore(HH);
        
        if DC
            HH = [HH ones(length(HH(:,1)),1)];
        end       
        
        if ARord
            %dydt = y;
            dydt = [(y(2:end)-y(1:end-1)) 0];
            %dydt = [0 (y(3:end)-y(1:end-2))/2 0];            
            
            %HH = [HH [0 y(1:end-1)]'];
            HH = [HH [-(y(2:end)+y(1:end-1))/2 0]'];
            %HH = [HH -y'];
        end
        
        covMat{p} = covMat{p} + HH'*HH;        
        
        xCorr{p} = HH'*dydt(:) + xCorr{p};


    end

    %covMat{p} = covMat{p}.*eye(size(covMat{p}));
    params{p} = (covMat{p})\xCorr{p};
    %params{p} = xCorr{p};
    
    paramsdum = params{p};

    alpha = params{p}(end-ARord+1:end);
    if ~ARord
        alpha = 0;
    end

    %tauConst(p) = acqPeriod/-log(alpha(1));
    tauConst(p) = 1/alpha;
    
    paramsdum = params{p}(1:end-ARord);
    xCorr{p} = xCorr{p}(1:end-ARord);

    if DC
        xCorr{p} = xCorr{p}(1:end-1);
        paramsdum = paramsdum(1:end-1);  %Get rid of DC shift
    end

    kerndum = reshape(paramsdum,paramP+blnkflag,tauP);
    if blnkflag
        kerndum = kerndum(1:end-1,:);
        blnkResp = kerndum(end,:);
    end
    
   
    for i = 1:length(kerndum(1,:))  %loop through each time point
        dum = kerndum(:,i);
        kern{p}(:,:,i) = reshape(dum,length(sfdom),length(oridom))';  %first dimension will be sf
        blnkResp(i) = mean(kern{p}(:,end,i),1); %Mean over orientation at highest sfreq        
        
        %if blnkflag
            
            %kern{p}(:,:,i) = kern{p}(:,:,i) - blnkResp(i);  %Use highest sp freq as blank
        %end
    end

    %allResp(:,p) = mean(kerndum)-blnkResp; %Mean over orientation at highest sfreq
    
    %subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

    %plotkerns(kern{p},oridom,sfdom,taudom,acqPeriod,blnkResp,alpha)
    %figure, plot(blnkResp)
    %imagesc(fliplr(Xcorrkern))

    %plot(kern(6,:))
end

%%

trialdom = 12;
for p = 16:16

    pID = nID(p);
    [idcelly idcellx] = find(masklabel == celldom(p));
    CoM = [mean(idcelly) mean(idcellx)];  %center of mass
    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);
    tauConst(p) = NaN;    

    for trialid = 1:length(trialdom)

        T = trialdom(trialid);
        [cond rep] = getcondrep(T);
        
        hper = gethper(cond);         
        Tupdate = Tf*hper;   

        y = squeeze(cellMat{cond}(pID,:,rep)); 
        
        tdom = (0:length(y)-1)*acqPeriod;
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain
        %of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;        

        %y = squeeze(mean(cellMat{cond}(:,:,rep),1));
               
        y = processTcourse(y,hh,polyorder,acqPeriod);
        
        
        %%%%
%         y = zscore(y);
%         acov = fftshift(xcov(y,'unbiased'));
%         noise = (acov(end) - acov(end-1))*length(acov); %estimate of noise variance; assuming its white, additive, and independent from everything else
%         N = noise*ones(size(y));
%         S = conj(fft(y)).*fft(y);
%         Ctrx = exp(-tdom/800);
%         H = fft(Ctrx);
%         G = (conj(H).*S)./(S.*abs(H).^2 + N);
%         y = ifft(G.*fft(y));
        %y = ifft(fft(y)./fft(Ctrx));
        %%%%
%         

        HHdum = zeros(length(y),length(oridom)*length(sfdom));
        pid = 1;
        for ori = 1:length(oridom)
            for sf = 1:length(sfdom)

                id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq == sfdom(sf));
                %id = find(seqs{T}.oriseq == oridom(ori));

                if ~isempty(id)

                    stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times (ms)
                    %stimes = synctimes{cond,rep}(id)*1000;

                    for i = 1:length(stimes)
                        [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                        HHdum(idx1,pid) = 1;
                    end
                end
                pid = pid+1;
            end
        end
        
        if blnkflag

            blnkvec = zeros(length(HHdum(:,1)),1);
            blnkid = find(isnan(seqs{T}.oriseq));
            stimes = synctimes{cond,rep}(blnkid)*1000;
            for i = 1:length(stimes)
                [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                blnkvec(idx1) = 1;
            end
            HHdum = [HHdum blnkvec];
            
        end
        
        %Truncate to timepoints with stimulation
        sid = sum(HHdum,2);
        id = find(sid);
        HHdum = HHdum(id(1):id(end),:);
        y = y(id(1):id(end));
        
        %Create stimulus matrix
        HH = zeros(length(HHdum(:,1)),length(HHdum(1,:))*tauP); %preallocate
        HHdum = [.0*randn(tauP-1,length(HHdum(1,:))); HHdum];  %Pad with zeros
        for z = tauP:length(HHdum(:,1))
            chunk = squeeze(HHdum(z-tauP+1:z,:))';
            chunk = chunk(:)';
            HH(z-tauP+1,:) = chunk;
        end
        %HH = zscore(HH);
        
        if DC
            HH = [HH ones(length(HH(:,1)),1)];
        end        
            
        MA = HH*params{p}(1:end-ARord);            
        
        if ARord
            alpha = params{p}(end);
            yhat = 0;
            for i = 1:length(MA);  %Now the autoregression
                yhat(i+1) = yhat(i)*(1-alpha) + MA(i);
            end
            yhat = yhat(2:end);
        else
            yhat = MA;
        end
        
        figure
        tdom = (0:(length(y)-1))*acqPeriod;
        plot(tdom,(y))
        hold on
        plot(tdom,(yhat),'r')
 
        if ARord>0
            hold on
            plot(tdom,(MA)/2-std(MA)*4,'k')
            legend('observed Fluor','predicted Ca','predicted spikes')
        end
        title(num2str(tauConst(p)))
        drawnow
        
        figure,
        scatter(yhat,y,'.k')
        xlabel('predicted response'), ylabel('actual response')
        R = corrcoef(y(:),yhat(:));
        title(['r = ' num2str(R(1,2))])

    end

    
end

function y = processTcourse(y,hh,polyorder,sp)

%mu = mean(y);
%First subtract a low-order polynomial fit:
yfit = polyfitter(y,polyorder);
y = y-yfit';

%Linear Bandpass Filter from GUI
if ~isempty(hh)
    y = ifft(fft(y).*hh);
end

%y = zscore(y);

%figure,plot(abs(fft(y)))

%Get rid of any peaks around the breathing rate
yP = getResPeaks(y,sp,10000);
fdom = linspace(0,1/(sp/1000),length(y)+1);
fdom = fdom(1:end-1);
atten = [.4 .2 .15 .1 .15 .2 .4];
aN = length(atten);
idbreath = find(fdom>1 & fdom<2.5);
Np = 0;
Npeaks = 3;
%yf = getResPeaks(y,sp,10000);
for i = 1:Npeaks
    yf = fft(y);
    [ma id] = max(abs(yf(idbreath)));  %Use yP to find the peaks
    if ma > 3.5*std(abs(yf(idbreath))) + mean(abs(yf(idbreath)));
        idpeak = id+idbreath(1)-1;
        yf(idpeak-floor(aN/2):idpeak+floor(aN/2)) = yf(idpeak-floor(aN/2):idpeak+floor(aN/2)).*atten;
        yP(idpeak-floor(aN/2):idpeak+floor(aN/2)) = yP(idpeak-floor(aN/2):idpeak+floor(aN/2)).*atten;
        yf = fliplr(yf);
        [dum id] = max(abs(yf(idbreath)));
        idpeak = id+idbreath(1)-1;
        yf(idpeak-floor(aN/2):idpeak+floor(aN/2)) = yf(idpeak-floor(aN/2):idpeak+floor(aN/2)).*atten;
        yP(idpeak-floor(aN/2):idpeak+floor(aN/2)) = yP(idpeak-floor(aN/2):idpeak+floor(aN/2)).*atten;
        yf = fliplr(yf);        
        y = real(ifft(yf));
        Np = Np+1;
    else
        break
    end

end

%"Adaptive line enhancement" to remove heart rate
% ybp = LFPfilt4(y,0,1000/(sp),.375-.02,.375);
% ord = 2; del = 0; Win = round(2000/sp);
% ybp = LinPred(ybp,sp,ord,5,Win);
% y = y-ybp;


noise = fftshift(xcov(y,'unbiased'));
noise = noise(end) - noise(end-1); %estimate of noise variance; assuming its white, additive, and independent from everything else
%[y noise] = wiener2(y, [1 round(110/sp)],noise);

%hold on,
%plot(abs(fft(y)),'r')
%title(num2str(Np))

%Nonlinear highpass filter to subtract estimated baseline
h = fspecial('gaussian', [1 length(y)], 10); %Use stats from heavily smoothed version
ydum = ifft(fft(y-mean(y)).*abs(fft(h)));
y = y - ordfilt2(ydum, 20, ones(1,600));

y = zscore(y);

y = y - prctile(y,10);

%y = y-mean(y);

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

function x = invertNL(x,y,yhat)


Nbins = 5;
ptsbin = floor(length(y_all{p})/Nbins);
[yhatS id] = sort(yhat_all{p});
yS = y_all{p}(id);
for i = 1:Nbins
    idsamp = ((i-1)*ptsbin+1):i*ptsbin;
    if i == Nbins
        idsamp = ((i-1)*ptsbin+1):length(y_all{p});
    end
    sampyhat = yhatS(idsamp);
    sampy = yS(idsamp);
    
    hatdum = sampyhat-trimmean(sampyhat,10);
    dum = sampy - trimmean(sampy,10);
    slpe = regress(dum,hatdum);
    
    
    muyhat(i) = trimmean(sampyhat,10);
    muy(i) = trimmean(sampy,10);
    sigyhat(i) = nanstd(sampyhat)/sqrt(length(sampyhat));
    sigy(i) = nanstd(sampy)/sqrt(length(sampy));
end


function plotkerns(kern,oridom,sfdom,taudom,acqPeriod,blnkResp,alpha)

kernplot = fliplr(squeeze(mean(kern(:,1:2,:),2)));  %Mean over sfreq... Get ori/time kernel
smoother = zeros(size(kernplot));
s = [.5 1 .5]';
s = s/sum(s);
smoother(1:length(s),:) = s*ones(1,length(kernplot(1,:)));
smoother = abs(fft(smoother));
kernplot = ifft(fft(kernplot).*smoother);
kernplot = real(kernplot);

figure
subplot(3,3,1)
%kernplot = kernplot-ones(length(oridom),1)*blnkResp;
imagesc(0:length(oridom)-1, 0:length(taudom)-1,kernplot')
labs = str2num(get(gca,'XTickLabel'));
set(gca,'XTickLabel',labs*(oridom(2)-oridom(1)));
labs = str2num(get(gca,'YTickLabel'));
set(gca,'YTickLabel',round(labs*acqPeriod)/1000);
set(gca,'TickDir','out')

[idy idx] = find(kernplot == max(kernplot(:)));
subplot(3,3,2)
id1 = max([idx-1 1]);
id2 = min([idx+1 length(kernplot(1,:))]);
kdum = kernplot(:,[id1 idx id2]);
plot(oridom,kdum), xlim([oridom(1) oridom(end)]), %legend('-T','peak','+T')
for i = 1:length(oridom)
    orilab{i} = num2str(oridom(i));
end
set(gca,'XTick',oridom,'XTickLabel',orilab)
ylim([min(kdum(:))-.05 max(kdum(:))+.05])
subplot(3,3,3)
plot(taudom,mean(kernplot(idy-1:idy+1,:),1),'.-k'), ylim([-.05 .21])
xlim([taudom(1) taudom(end)])

kernplot = fliplr(squeeze(mean(kern(3:6,:,:),1)));  %Get sf/time kernel
%kernplot = fliplr(squeeze(mean(kern(6:12,:,:))));  %Get sf/time kernel
smoother = zeros(size(kernplot));
s = [.0 1 .0]';
s = s/sum(s);
smoother(1:length(s),:) = s*ones(1,length(kernplot(1,:)));
smoother = abs(fft(smoother));
kernplot = ifft(fft(kernplot).*smoother);

subplot(3,3,4)
%kernplot = kernplot-ones(length(sfdom),1)*blnkResp;
imagesc(0:length(sfdom)-1,0:length(taudom)-1,kernplot')
labs = str2num(get(gca,'XTickLabel'));
set(gca,'XTickLabel',sfdom(1)+labs*(sfdom(2)-sfdom(1)));
labs = str2num(get(gca,'YTickLabel'));
set(gca,'YTickLabel',round(labs*acqPeriod)/1000);
set(gca,'TickDir','out')

[idy idx] = find(kernplot == max(kernplot(:)));
subplot(3,3,5)
id1 = max([idx-1 1]);
id2 = min([idx+1 length(kernplot(1,:))]);
kdum = kernplot(:,[id1 idx id2]);
plot(sfdom,kdum), xlim([sfdom(1) sfdom(end)]) %legend('-T','peak','+T')
for i = 1:length(sfdom)
    sflab{i} = num2str(sfdom(i));
end
set(gca,'XTick',sfdom,'XTickLabel',sflab)
ylim([min(kdum(:))-.05 max(kdum(:))+.05])
subplot(3,3,6)
plot(taudom,mean(kernplot(idy:idy,:),1),'.-k'), ylim([-.05 .21])
xlim([taudom(1) taudom(end)])
drawnow
    
%Make the Ca transient
% ARord = length(alpha);
% clear y
% N = 50;
% x = [zeros(1,ARord) 1 zeros(1,N-1)];
% y = zeros(1,ARord);
% for i = ARord+1:length(x)
%     y(i) = x(i);
%     for j = 1:ARord
%         y(i) = y(i) + y(i-j)*alpha(j);
%     end
% end
% tdom = (0:(length(y)-1))*acqPeriod;
% tdom = tdom - acqPeriod*ARord;
% subplot(3,3,7),plot(tdom,y,'.-k')
% [dum id] = min(abs(y-exp(-1)));
% tau = tdom(id);
% title(['tau = ' num2str(tau)])
% xlim([tdom(1) tdom(end)])