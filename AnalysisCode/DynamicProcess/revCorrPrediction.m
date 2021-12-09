function kern = revCorrPrediction(cellMat,trialdom_build,trialdom_predict)

global ACQinfo maskS Analyzer G_RChandles idExamp

DC = 1;
ARflag = 1;
polyorder = 1;

nID = getNeuronMask;  %get the index values for the neurons
masklabel = bwlabel(maskS.neuronmask,4);
celldom = unique(masklabel);
Ncell = length(nID);

ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 

logfileroot = get(G_RChandles.logfilePath,'string');
expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];

load([logfileroot Analyzer.M.anim '\' expt],'frate')

if ~exist('frate')
    frate = 60.0044;
end

Tf = 1000/frate;  %Frame period in ms (frate obtained from log file) 

[domains seqs] = getSeqInfo(trialdom_build);

oridom = domains{trialdom_build(1)}.oridom;
sfdom = domains{trialdom_build(1)}.sfdom;
phasedom = domains{trialdom_build(1)}.phasedom;

paramP = length(oridom)*length(sfdom);
if ~ARflag
    tauP = round(1200/acqPeriod)+1;
else
    tauP = round(250/acqPeriod)+1;
end
%tauP = 40;
taudom = 0:acqPeriod:((tauP-1)*acqPeriod);
covMat = 0;
xCorr = cell(1,Ncell);
covMat = cell(1,Ncell);
for i = 1:Ncell
    xCorr{i} = 0;
    covMat{i} = 0;
end

hh = makeTemporalfilter;

%idExamp = [6 27]
%idExamp = [79 59 33]
%idExamp = [6 34 54];  %ab3 002_041 (trial 10), ab2  000_068
%idExamp = [9 54 78 83];  %aa9 001_027
%idExamp = [35 54];  %ab2 000_014

%idExamp = [6 7 9 10 13]; %zc1 025
%idExamp = [1 7]; %zc1 017

%First build the kernels
figure
for p = 1:length(idExamp)

    pID = nID(idExamp(p));
    [idcelly idcellx] = find(masklabel == celldom(idExamp(p)));
    CoM = [mean(idcelly) mean(idcellx)];  %center of mass
    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);
    tauConst(p) = NaN;

    for trialid = 1:length(trialdom_build)

        T = trialdom_build(trialid);
        [cond rep] = getcondrep(T);
        
        hper = gethper(cond);         
        Tupdate = Tf*hper;   

        y = squeeze(cellMat{cond}(pID,:,rep));     

        %y = squeeze(mean(cellMat{cond}(:,:,rep),1));
        y = processTcourse(y,hh,polyorder,acqPeriod);
        
        tdom = (0:length(y)-1)*acqPeriod;
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain
        %of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;

        HHdum = zeros(length(y),length(oridom)*length(sfdom));
        pid = 1;
        for ori = 1:length(oridom)
            for sf = 1:length(sfdom)

                id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq == sfdom(sf));
                %id = find(seqs{T}.oriseq == oridom(ori));

                if ~isempty(id)

                    stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times (ms)

                    for i = 1:length(stimes)
                        [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                        HHdum(idx1,pid) = 1;
                    end
                end
                pid = pid+1;
            end
        end
        %HHdum = zscore(HHdum);
        HH = zeros(length(HHdum(:,1)),length(HHdum(1,:))*tauP); %preallocate
        HHdum = [.0*randn(tauP-1,length(HHdum(1,:))); HHdum];  %Pad with zeros
        for z = tauP:length(HHdum(:,1))
            chunk = squeeze(HHdum(z-tauP+1:z,:))';
            chunk = chunk(:)';
            HH(z-tauP+1,:) = chunk;
        end
        
        if DC
            HH = [HH ones(length(HH(:,1)),1)];
        end
        if ARflag
            
            %dydt = [(y(2:end)-y(1:end-1)) 0];
            %dydt = [0 (y(3:end)-y(1:end-2))/2 0];                 
            %HH = [HH dydt'];           
            
            HH = [HH [0; y(1:end-1)']];
        end

        covMat{p} = covMat{p} + HH'*HH;

        xCorr{p} = HH'*y(:) + xCorr{p};

    end
    
    if ~ARflag
        covMat{p} = covMat{p}.*eye(size(covMat{p}));
    end
    
    params{p} = (covMat{p})\xCorr{p};    
    
    %params{p} = xCorr{p};
    
    paramsdum = params{p};
    if ARflag
        alpha = params{p}(end);
        tauConst(p) = acqPeriod/-log(alpha);
        paramsdum = params{p}(1:end-1);
        xCorr{p} = xCorr{p}(1:end-1);
    end
    if DC
        xCorr{p} = xCorr{p}(1:end-1);
        paramsdum = paramsdum(1:end-1);  %Get rid of DC shift
    end

    kerndum = reshape(paramsdum,paramP,tauP);

    for i = 1:length(kerndum(1,:))
        dum = kerndum(:,i);
        kern{p}(:,:,i) = reshape(dum,length(sfdom),length(oridom))';  %first dimension will be ori
    end
    
%This is the same thing as subtracting the response at time zero while
%building the kernel...
%     for i = 1:length(oridom)
%         for j = 1:length(sfdom)
%             kern{p}(i,j,:) = kern{p}(i,j,:) - kern{p}(i,j,1);
%         end
%     end
    

    %subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    plotkerns(kern{p},oridom,sfdom,taudom,acqPeriod,tauConst(p))
    
    %imagesc(fliplr(Xcorrkern))

    %plot(kern(6,:))
end


%% Now Prediction
figure
hh = makeTemporalfilter;
%pdom = 2:2; %Id of kernels that were extracted above... index in 'idExamp'
pdom = 1:length(idExamp);
for ploop = 1:length(pdom)
    p = pdom(ploop);
    
    pID = nID(idExamp(p));
    [idcelly idcellx] = find(masklabel == celldom(idExamp(p)));
    CoM = [mean(idcelly) mean(idcellx)];  %center of mass
    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    
    yall = [];
    yhatall = [];
    for trialid = 1:length(trialdom_predict)
        T = trialdom_predict(trialid);
        [cond rep] = getcondrep(T);
        
        hper = gethper(cond);         
        Tupdate = Tf*hper;   

        y = squeeze(cellMat{cond}(pID,:,rep));
        
        %y = processTcourse(y,hh,polyorder,acqPeriod,normdom(ploop));
        y = processTcourse(y,hh,polyorder,acqPeriod);

        tdom = (0:length(y)-1)*acqPeriod;
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain
        %of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;

        HHdum = zeros(length(y),length(oridom)*length(sfdom));
        pid = 1;
        for ori = 1:length(oridom)
            for sf = 1:length(sfdom)

                id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq == sfdom(sf));
                %id = find(seqs{T}.oriseq == oridom(ori));

                if ~isempty(id)

                    stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times (ms)

                    for i = 1:length(stimes)
                        [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                        HHdum(idx1,pid) = 1;
                    end
                end
                pid = pid + 1;
            end
        end

        HH = zeros(length(HHdum(:,1)),length(HHdum(1,:))*tauP); %preallocate
        HHdum = [zeros(tauP-1,length(HHdum(1,:))); HHdum];  %Pad with zeros
        for z = tauP:length(HHdum(:,1))
            chunk = squeeze(HHdum(z-tauP+1:z,:))';
            chunk = chunk(:)';
            HH(z-tauP+1,:) = chunk;
        end
        
        if DC
            HH = [HH ones(length(HH(:,1)),1)];
        end
        if ARflag
            MA = HH*params{p}(1:end-1);  %The moving average (i.e. spike rate)
            yhat = 0;
            for i = 1:length(MA);  %Now the autoregression
                yhat(i+1) = yhat(i)*params{p}(end) + MA(i);
            end
            yhat = yhat(2:end);
        else
            yhat = HH*params{p};
        end
        
        %MA = MA(50:end-100)/max(MA)*max(yhat)/2;
        strt = round(getparam('predelay')/(acqPeriod/1000));
        stp = round((getparam('predelay')+getparam('stim_time'))/(acqPeriod/1000));
        yhat = yhat(strt:stp);
        y = y(strt:stp);
        if ARflag
            MA = MA(strt:stp);
        end
        tdom = (0:length(y)-1)*acqPeriod/1000;
        
        figure
        %subplot(length(pdom),1,ploop)
        %subplot(2,2,ploop*2)
        
        id = find(y<prctile(y,50) & y>prctile(y,0));
        y = y-median(y); 
        y = y/std(y(id));
        id = find(yhat<prctile(yhat,50) & yhat>prctile(yhat,0));
        yhat = yhat-median(yhat(id)); 
        yhat = yhat/std(yhat(id));
        
        plot(tdom,(y),'k')
        hold on
        plot(tdom,(yhat),'r')
        xlabel('seconds'), ylabel('standard deviation')
        if ARflag
            hold on
            plot(tdom,zscore(MA)/2-4,'k')
            legend('Observed fluorescence','predicted Ca2+','predicted spikes')
        else
            legend('Observed fluorescence','Predicted Ca2+')
        end
        title(num2str(tauConst(p)))
        ylim auto
        drawnow
        hold off
        
%         figure(22)
%         subplot(length(pdom),1,ploop)
%         [hst hdom] = hist(y,[-3:.5:5]);
%         bar(hdom,hst/sum(hst)), set(gca,'TickDir','out')
%         xlim([hdom(1) hdom(end)])
%         
        [muy muyhat sigy] = makeEbars(y,yhat);
        
%         figure,
%         scatter(yhat,y,'.k')
%         xlabel('predicted response'), ylabel('actual response')
%         
%         [param] = NL(yhat(:),y(:),[1 muyhat(1) 1 muy(1)]);        
%         domu = unique(yhat);
%         ffit = param(1)*phi(domu-param(2)).^param(3) + param(4);
%         
%         %hold on
%         %errorbar(muyhat,muy,sigy,'r','linewidth',3)     
%         hold on
%         plot(domu,ffit,'b')
%         R = corrcoef(y(:),yhat(:));
%         title(['r = ' num2str(R(1,2))])
        
        yall = [yall; y(:)];
        yhatall = [yhatall; yhat(:)];
%         
%         yhat2 = param(1)*phi(yhat-param(2)).^param(3) + param(4);
        
        
%         figure
%         plot(tdom,y)
%         hold on
%         plot(tdom,yhat2,'r')
%         legend('actual','predicted')
%         xlabel('sec')
%         ylabel('Z')
%         title('after NL')
        
        %corrcoef(y(:),yhat2(:));
        
        
    end
    
    [muy muyhat sigy] = makeEbars(yall,yhatall);
    figure,
    scatter(yhatall,yall,'.k')
    xlabel('predicted response'), ylabel('actual response')
    
    [param] = NL(yhatall(:),yall(:),[1 muyhat(1) 1 muy(1)]);
    domu = unique(yhatall);
    ffit = param(1)*phi(domu-param(2)).^param(3) + param(4);
    
    %hold on
    %errorbar(muyhat,muy,sigy,'r','linewidth',3)
    hold on
    plot(domu,ffit,'b')
    R = corrcoef(y(:),yhat(:));
    title(['r = ' num2str(R(1,2))])
    
    
end


function xfit = polyfitter(x,order)

dom = (0:length(x)-1)';

H = ones(length(dom),order+1);  %last column for DC
for i = 1:order
    H(:,i) = dom.^i;
end

p = inv(H'*H)*H'*x';
xfit = H*p;

function plotkerns(kern,oridom,sfdom,taudom,acqPeriod,tauConst)

kernplot = fliplr(squeeze(mean(kern(:,1:end,:),2)));  %Mean over sfreq... Get ori/time kernel
smoother = zeros(size(kernplot));
s = [.5 1 .5]';
smoother(1:length(s),:) = s*ones(1,length(kernplot(1,:)));
smoother = abs(fft(smoother));
kernplot = ifft(fft(kernplot).*smoother);

figure
subplot(4,2,1)
imagesc(0:length(taudom)-1,0:length(oridom)-1,kernplot)
labs = (get(gca,'YTickLabel'));
for i = 1:length(labs)
   labsnum(i) = str2num(labs{i});
end
set(gca,'YTickLabel',labsnum*(oridom(2)-oridom(1)));
labs = (get(gca,'XTickLabel'));
clear labsnum
for i = 1:length(labs)
   labsnum(i) = str2num(labs{i});
end
set(gca,'XTickLabel',round(labsnum*acqPeriod)/1000);
title(num2str(tauConst))

[idy idx] = find(kernplot == max(kernplot(:)));
subplot(4,2,2)
plot(oridom,kernplot(:,idx-1:idx+1)), xlim([oridom(1) oridom(end)]), %legend('-T','peak','+T')
xlabel('ori')
subplot(4,2,3)
plot(taudom,mean(kernplot(idy:idy,:),1))
xlabel('ms')

kernplot = fliplr(squeeze(mean(kern(1:end,:,:))));  %Get sf/time kernel
smoother = zeros(size(kernplot));
s = [.5 1 .5]';
smoother(1:length(s),:) = s*ones(1,length(kernplot(1,:)));
smoother = abs(fft(smoother));
kernplot = ifft(fft(kernplot).*smoother);

subplot(4,2,5)
imagesc(0:length(taudom)-1,1:length(sfdom),kernplot)
labs = get(gca,'YTickLabel'); 
clear labsnum
for i = 1:length(labs)
   labsnum(i) = str2num(labs{i});
end
set(gca,'YTickLabel',round(sfdom(labsnum)*100)/100);
labs = get(gca,'XTickLabel');
clear labsnum
for i = 1:length(labs)
   labsnum(i) = str2num(labs{i});
end
set(gca,'XTickLabel',round(labsnum*acqPeriod)/1000);


[idy idx] = find(kernplot == max(kernplot(:)));
subplot(4,2,6)
plot(log2(sfdom),kernplot(:,idx:idx)), %legend('-T','peak','+T')
xlabel('log_2(sf)')
xlim([log2(sfdom(1)) log2(sfdom(end))])
subplot(4,2,7)
plot(taudom,mean(kernplot(idy:idy,:),1))
xlabel('ms')
drawnow
    
function [muy muyhat sigy] = makeEbars(y,yhat)

Nbins = 5;
ptsbin = floor(length(y(:))/Nbins);
[yhatS id] = sort(yhat(:));
yS = y(id);
for i = 1:Nbins
    idsamp = ((i-1)*ptsbin+1):i*ptsbin;
    if i == Nbins
        idsamp = ((i-1)*ptsbin+1):length(y(:));
    end
    sampyhat = yhatS(idsamp);
    sampy = yS(idsamp);

    muyhat(i) = trimmean(sampyhat,10);
    muy(i) = trimmean(sampy,10);
    sigyhat(i) = nanstd(sampyhat)/sqrt(length(sampyhat));
    sigy(i) = nanstd(sampy)/sqrt(length(sampy));
end
