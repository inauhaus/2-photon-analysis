function getRFfromRandPos4
%2 computes the fit for the 'on', 'off', and 'on' + 'off' kernel

%3 was used for the paper: it computes size from the 1d profile at the
%optimal orientation

global ACQinfo G_RChandles cellS DM TC idExamp

Monkbit = 1;
if Monkbit
    %For monkey
    kernTsig = 50;
    kernOrisig = 10;
    kernXsig = .05;
    kernTsigMa = 50;
    kernOrisigMa = 15;
    kernXsigMa = .3;
    % maDel = getparam('h_per')*10 + 150;
    % delayWin = [100 maDel];
    delayWin = [50 400];  %assume the peak response is within this time window
else
    %For mouse
    kernTsig = 50;
    kernOrisig = 20;
    kernXsig = 2;
    kernTsigMa = 50;
    kernOrisigMa = 15;
    kernXsigMa = 2;
    delayWin = [50 500];  %assume the peak response is within this time window
end
%%%%

%Downsample kernel and domains
DSx = round(getparam('barWidth')/(getparam('x_size')/getparam('Nx')));
if DSx>1
    'Downsampling xpos'
end
DSo = 0;
if getparam('n_ori') >= 16
    if ~rem(getparam('n_ori'),2);
        DSo = 2;
        'Downsampling ori'
    end
end

%if ~isempty(find(isnan(kernAll{1})))

trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
getCountMat(trialdom,cellS.cellMat) % loads kernCount

eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with acqPeriod spacing, and end at an estimate of kernDel(2)

kernCount = cellS.kernCount;
kernAll = cellS.kernAll;
kernSigAll = cellS.kernSigAll;


for c = 1:length(DM.colordom)
    for kElem = 1:3
        for p = 1:length(kernAll)
            id = find(taudom>100 & taudom<400);
            dumB = mean(kernAll{p}(:,:,1,:,id),5);
            dumW = mean(kernAll{p}(:,:,2,:,id),5);
            Bamp = nanmean(dumB(:).^2);
            Wamp = nanmean(dumW(:).^2);
            Bamp = max(dumB(:));
            Wamp = max(dumW(:));
            
            TC.BWdiff{c}{kElem}(p) = (Bamp-Wamp)/(Bamp+Wamp);  %its just convenient to compute this here
            %TC.BWdiff{c}(p) = log2(norm(dumB(:))./norm(dumW));  %its just convenient to compute this here
        end
    end
end

DSbw = 0; %Flag to combine black and white
[kernSigAll kernCountdum] = DSrandposkernel(kernSigAll,kernCount,DSx,DSo,DSbw);
[kernAll kernCount] = DSrandposkernel(kernAll,kernCount,DSx,DSo,DSbw);

if DSx>1
     DM.xdom = DM.xdom(2:2:end);
end

if DSo
%     dori = median(diff(DM.oridom));
%     dum = exp(1i*2*(DM.oridom-dori/2)*pi/180);
%     dum = angle(dum)/2*180/pi;
%     id = find(dum<0);
%     dum(id) = dum(id)+180;
%     DM.oridom = dum;

    DM.oridom = DM.oridom(1:2:end);
end

for c = 1:length(DM.colordom)
    for p = 1:length(kernAll)
        kernC{c}{p} = squeeze(kernAll{p}(:,:,:,c,:));
        kernSigC{c}{p} = squeeze(kernSigAll{p}(:,:,:,c,:));
        countC{c}{p} = squeeze(kernCount{p}(:,:,:,c,:));
    end
end


[dum delayWinID(1)] = min(abs(delayWin(1)-taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-taudom));

Ncell = length(kernAll);

%Get the receptive field

kernsmooth = getSmoother(kernTsigMa,kernOrisigMa,kernXsigMa,DM.taudom,DM.oridom,DM.xdom); %used to find maxima

param(5) = -36;

for c = 1:length(DM.colordom)

    for p = 1:Ncell
        
        %%%%Make OFF kernel%%%%%%%
        kernplot_b = squeeze(kernC{c}{p}(:,:,1,:));
        countplot_b = squeeze(countC{c}{p}(:,:,1,:));
        kernsigplot_b = squeeze(kernSigC{c}{p}(:,:,1,:))./sqrt(countplot_b);
        
        kernplot_b = kernplot_b.*countplot_b; %unnormalize
        kernplot_b = smoothkern(kernplot_b,kernTsig,kernOrisig,kernXsig,DM.taudom,DM.xdom); %smooth
        countplot_b = smoothkern(countplot_b,kernTsig,kernOrisig,kernXsig,DM.taudom,DM.xdom);
        kernplot_b = kernplot_b./countplot_b; %re-normalize
        
        %%%%Make ON kernel%%%%%%%
        kernplot_w = squeeze(kernC{c}{p}(:,:,2,:));
        countplot_w = squeeze(countC{c}{p}(:,:,2,:));
        kernsigplot_w = squeeze(kernSigC{c}{p}(:,:,2,:))./sqrt(countplot_w);
        
        kernplot_w = kernplot_w.*countplot_w; %unnormalize
        kernplot_w = smoothkern(kernplot_w,kernTsig,kernOrisig,kernXsig,DM.taudom,DM.xdom); %smooth
        countplot_w = smoothkern(countplot_w,kernTsig,kernOrisig,kernXsig,DM.taudom,DM.xdom);
        kernplot_w = kernplot_w./countplot_w; %re-normalize
        
        %%%%Make ON+OFF kernel%%%%%%%
        kernplot_bw = (kernplot_b.*countplot_b + kernplot_w.*countplot_w)./(countplot_w+countplot_b);
        kernsigplot_bw = (kernsigplot_b.*countplot_b + kernsigplot_w.*countplot_w)./(countplot_w+countplot_b);
        
        %%%%Get time slice from ON+OFF kernel%%%%%%%
        kerndum = ifftn(fftn(kernplot_bw).*abs(fftn(kernsmooth)));
        
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        
        [ma idma] = max(kerndum,[],3); %find maxima from smoothed version
        [bestoriid bestposid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestposid) + delayWinID(1) - 1;
        
        tcourseBestMu{c}{p} = squeeze(kernplot_bw(bestoriid,bestposid,:));
        tcourseBestSig{c}{p} = squeeze(kernsigplot_bw(bestoriid,bestposid,:));
        
        [mi idmi] = min(kerndum,[],3); %find minima from smoothed version
        [worstoriid worstposid] = find(mi == min(mi(:)));
        tcourseWorstMu{c}{p} = squeeze(kernplot_bw(worstoriid,worstposid,:));
        tcourseWorstSig{c}{p} = squeeze(kernsigplot_bw(worstoriid,worstposid,:));
        
        kernplotIm{1}{c}{p} = kernplot_bw(:,:,tau);
        kernplotIm{2}{c}{p} = kernplot_b(:,:,tau);
        kernplotIm{3}{c}{p} = kernplot_w(:,:,tau);
        
    end
end

%%

for c = 1:length(DM.colordom)

    idex = 1;
    figure
    for p = 1:Ncell
        
        subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
        
        for kernElement = 1:3 %[b+w b w]
                
            im = kernplotIm{kernElement}{c}{p};
            degperpix = getparam('x_size')/(length(im(1,:))-1);  %can't use xdom because of the padding
            
            [bestoriid bestposid] = find(im == max(im(:)));
            
            rawProfile = im(bestoriid,:);
            [param profileFit] = Gaussfit(1:length(rawProfile),rawProfile,0);
            param(2) = abs(param(2));  %This is the same since it gets squared
            
            profileSize = 2*param(2)*degperpix; %2sigma (degrees)
            
            %     for i = 1:length(DM.oridom)
            %         tc = im_b(i,:) + im_b(i,:);
            %         [param ffit varaccount] = Gaussfit(1:length(tc),tc,0);
            %         imfit(:,i) = ffit;
            %     end
            
            
            imfit = im;
            
            amp = norm(imfit(:));  %do this before any normalization
            
            %if q == 1
            thresh1 = prctile(imfit(:),2);
            %end
            imfit = phi(imfit-thresh1);
            padder = zeros(length(imfit(:,1)),0);
            imfit = [padder imfit padder]';
            %imfit(:,end/2:end) = 0;
            
            %%%%%%%Control%%%%%%%%%
            %         imfit = [imfit flipud(imfit)];
            %         shiftori = round(rand(1)*(size(imfit,2)-1));
            %         imfit = circshift(imfit,[0 shiftori]);
            %         imfit = imfit(:,1:size(imfit,2)/2);
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            RF = iradon(imfit,DM.oridom,'linear','none',1,length(imfit(:,1))); %if you don't do "filtered" back-projection, it smears the image
            
            %         hh = zeros(size(RF));
            %         hh(1:5,1:5) = [.1 .3 1 .3 .1]'*[.1 .3 1 .3 .1];
            %         RF = ifft2(fft2(RF).*abs(fft2(hh)));
            
            %if q == 1
            thresh2 = prctile(RF(:),20) + (max(RF(:))-prctile(RF(:),20))*.3;
            %thresh2 = prctile(RF(:),50);
            %thresh = median(RF(:));
            %end
            RF = phi(RF-thresh2);
            
            %     [TC.xpos{c}(p) TC.ypos{c}(p)] = getmaxLoc(RF,30);  %returns 'degrees'
            
            otc = max(imfit);
            
            [dum ohat] = orifind(otc',DM.oridom);
            
            [param ffit varaccount] = Gaussfit2Drot(RF,ohat); %This does not fit an amplitude or a baseline
            ffit = ffit*max(RF(:));
            [cc pp] = corrcoef(ffit(:),RF(:));
            
%             if param(1) > param(2) %if ysize > xsize.
%                 param(5) = param(5)+90;
%             end
            if param(5) < 0
                param(5) = param(5)+180;
            end
            if param(5) > 180
                param(5) = param(5)-180;
            end
            
            dataYielder = varaccount;
            dataYieldThresh = 0.6;
            
            if dataYielder > dataYieldThresh
                xpos = param(2)*degperpix;
                ypos = param(1)*degperpix;
                xsize = param(3)*degperpix*2;
                ysize = param(4)*degperpix*2;
                princAxis = param(5);  %Not necessarily the ori preference, but the axis of elongation.
                
                %pkamp = param(6);
                if xpos > 10 || ypos > 10 || xpos < -10 || ypos < -10
                    xpos = NaN;
                    ypos = NaN;
                    xsize = NaN;
                    ysize = NaN;
                end
            else
                xpos = NaN;
                ypos = NaN;
                xsize = NaN;
                ysize = NaN;
                profileSize = NaN;
                princAxis = NaN;
                %pkamp = NaN;
            end

            %             Opref{c}(p) = ohat;
            %             Opreffit{c}(p) = param(5);
            %             Omagfit{c}(p) = max([TC.xsize{c}(p) TC.ysize{c}(p)])/min([TC.xsize{c}(p) TC.ysize{c}(p)]);  %aspect ratio
            
            xdom = (0:length(RF(1,:))-1)*degperpix; %degrees
            ydom = (0:length(RF(:,1))-1)*degperpix;   %don't use 'y_size'
            %RF = RF/max(RF(:));
            if p == 1
                [yma xma] = find(RF == max(RF(:)));
   
                yran = 1:length(ydom);
                xran = 1:length(xdom);
            end
            dim = size(RF);
            %RF = RF(yran,xran); ffit = ffit(yran,xran); ydom = ydom(yran); xdom = xdom(xran);
            
            %%%Store all the important variables
            TC.xpos{c}{kernElement}(p) = xpos;
            TC.ypos{c}{kernElement}(p) = ypos;
            TC.xsize{c}{kernElement}(p) = xsize;
            TC.ysize{c}{kernElement}(p) = ysize;
            TC.princAxis{c}{kernElement}(p) = princAxis; 
            TC.amp{c}{kernElement}(p) = amp;
            %TC.amp{c}{kernElement}(p) = pkamp + thresh1;
            TC.RF{c}{kernElement}{p} = RF;
            TC.RFfit{c}{kernElement}{p} = ffit;
            TC.RFraw{c}{kernElement}{p} = im;
            TC.rawProfile{c}{kernElement}{p} = rawProfile;
            TC.profileSize{c}{kernElement}(p) = profileSize;
            TC.princAxis{c}{kernElement}(p) = princAxis;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ffit = ffit/max(RF(:));
            RF = RF/max(RF(:));
            cT = .6*max(RF(:));
            switch kernElement
                
                case 1   %the on+off
                    imagesc(xdom,ydom,RF,[0 1]), colormap gray %title(num2str(param(5)))
                    hold on
%                     contour(xdom,ydom,ffit,[cT cT],'k')
%                     axis image

                    drawnow
                case 2  %OFF
                    hold on
                    contour(xdom,ydom,ffit,[cT cT],'b','LevelStepMode','manual')
                    axis image
                    drawnow
                case 3   %ON
                    hold on
                    contour(xdom,ydom,ffit,[cT cT],'r','LevelStepMode','manual')
                    axis image
                    drawnow
                    
            end
            axis off
%             
%             imagesc(xdom(xran),ydom(yran),RF(yran,xran),[0 max(RF(:))]), colormap jet %title(num2str(param(5)))
%             hold on
%             contour(xdom(xran),ydom(yran),ffit(yran,xran),[.5*max(ffit(:)) .5*max(ffit(:))],'LineColor',[0 0 0],'LevelStepMode','manual')
%             axis image
%             axis off
%             drawnow
            
            
            if p == 1
                title([num2str(ydom(end)-ydom(1)) 'degrees'])
            end
            
           
            %plot([mean(kernplot_b(4:5,:))' mean(kernplot_w(4:5,:))'])
            
            
            %Get stuff to plot the examples
            if ~isempty(find(p == idExamp))
                
                fit_examp{c}{idex}{kernElement} = ffit;
                RF_examp{c}{idex}{kernElement} = RF;
                rawKern{c}{idex}{kernElement} = im;
                rawProf{c}{idex}{kernElement} = rawProfile;
                profFit{c}{idex}{kernElement} = profileFit;
                
                cc_examp{c}{idex}{kernElement} = cc(1,2);
                
                timecourseBestMu{c}{idex} =  tcourseBestMu{c}{p};
                timecourseBestSig{c}{idex} =  tcourseBestSig{c}{p};
                timecourseWorstMu{c}{idex} =  tcourseWorstMu{c}{p};
                timecourseWorstSig{c}{idex} =  tcourseWorstSig{c}{p};
                
                
                
            end
           
            
        end
        if ~isempty(find(p == idExamp))
            idex = idex+1;
        end
    end
end
axis image
%%

for kernElement = 1:3
    id = find(isnan(TC.xpos{1}{kernElement}));
    TC.amp{c}{kernElement}(id) = NaN; %amp is not based on a fit, so its not necessary to exclude as usual
end

%TC.BWdiff{c} = (TC.ampB{c}-TC.ampW{c})./(TC.ampB{c}+TC.ampW{c});  %its just convenient to compute this here 

%TC.BWdiff{c} = log2(norm(TC.ampB{c})./norm(TC.ampW{c}));  %its just convenient to compute this here 



%% Plot examples
colorstr = ['k' 'b' 'r']

if ~isempty(idExamp)
    for c = 1:length(DM.colordom)
        figure
        Nex = length(fit_examp{1});
        for i = 1:Nex
            
            if length(find(DM.oridom>180))
                nshift = floor(90/(DM.oridom(2)-DM.oridom(1)));
            else
                nshift = 0;
            end
            odum = circshift(DM.oridom,[1 nshift]);
            q = find(odum>180);
            odum(q) = odum(q)-180;
            odum = odum-odum(1);
            
            for kernElement = 1:3
                
                rawdum = [rawKern{c}{i}{kernElement}' flipud(rawKern{c}{i}{kernElement}')];
                rawdum = circshift(rawdum,[1 -nshift]);
                rawdum = rawdum(:,1:(length(DM.oridom)));
                
                subplot(Nex,5,(i-1)*5+(5-4))
                if kernElement == 1
                    imagesc([odum 180],xdom,[rawdum flipud(rawdum(:,1))]), colorbar, colormap gray
                    set(gca,'Xtick',[0 90 180])
                    axis square
                    title(num2str(cc_examp{c}{i}{kernElement}))
                    xlabel('ori')
                end
                
                subplot(Nex,5,(i-1)*5+(5-3))
                if kernElement == 1
                    imagesc(xdom,ydom,RF_examp{c}{i}{kernElement}), colorbar, colormap gray
                    hold on,
%                     contour(xdom,ydom,fit_examp{c}{i}{kernElement},[.61 .61]*max(fit_examp{c}{i}{kernElement}(:)),'LineColor',[0 0 0])
%                     caxis([min(RF_examp{c}{i}{kernElement}(:)) max(RF_examp{c}{i}{kernElement}(:))])                    
                else                    
                    hold on
                    contour(xdom,ydom,fit_examp{c}{i}{kernElement},[.61 .61]*max(fit_examp{c}{i}{kernElement}(:)),colorstr(kernElement),'LevelStepMode','manual')
                    axis image
                    drawnow                   
                end
                

                axis image
                
                subplot(Nex,5,(i-1)*5+(5-2))
                tcdum = mean(rawdum(:,:));
                
                [param] = Gaussfit(odum,tcdum,1);
                fitdom = 0:180;
                ffit = exp(-(fitdom-90).^2/(2*(param(2)).^2));
                ffit = circshift(ffit,[1 round(param(1))-90]);
                ffit = param(3)*ffit + param(4);
                %ffit = circshift(ffit,[0 round(param(1))-90]);
                
                if kernElement ~= 1
                    plot([odum odum(1)],[tcdum tcdum(1)],['.' colorstr(kernElement)])
                    hold on
                    plot(0:180,ffit,[colorstr(kernElement)])
                    xlim([0 180])
                    set(gca,'Xtick',[0 90 180])
                    xlabel('ori')
                end
                %ylim([-.1 .25])
                
                %             subplot(Nex+1,2,Nex*2+2)
                %             plot(mean(RF_examp{c}{i}),'k') hold on
                
                if kernElement ~= 1
                    xdomI = xdom(1):.01:xdom(end);
                    pFitI = interp1(xdom,profFit{c}{i}{kernElement},xdomI,'spline');
                    
                    subplot(Nex,5,(i-1)*5+(5-1))
                    plot(xdomI,pFitI, colorstr(kernElement))
                    hold on
                    plot(xdom,rawProf{c}{i}{kernElement},['.' colorstr(kernElement)])
                    hold on
                    
                   xlabel('pos (deg)')
                end
                
                 
                %ylim([-.5 1.5])
                
            end
            
            subplot(Nex,5,(i-1)*5+(5-0))
            fill([DM.taudom fliplr(DM.taudom)],[timecourseBestMu{c}{i}-timecourseBestSig{c}{i}; flipud(timecourseBestMu{c}{i}+timecourseBestSig{c}{i})]',[.0 .0 1])
            hold on
            plot(DM.taudom,timecourseBestMu{c}{i},'k'),xlabel('ms')
            
           hold on
            fill([DM.taudom fliplr(DM.taudom)],[timecourseWorstMu{c}{i}-timecourseWorstSig{c}{i}; flipud(timecourseWorstMu{c}{i}+timecourseWorstSig{c}{i})]',[1 .0 0])
            hold on
            plot(DM.taudom,timecourseWorstMu{c}{i},'k'),xlabel('ms')
            
            xlim([DM.taudom(1) DM.taudom(end)])        
            
            
            
        end
    end
end


%% Get ori tuning

% thet = 40;
% xposp = TC.xpos*cos(thet) + TC.ypos*sin(thet);
% yposp = -TC.xpos*sin(thet) + TC.ypos*cos(thet);

%Get and plot ori tuning 
colorstr = ['k' 'b' 'r'];
figure
for c = 1:length(DM.colordom)
    
    for p = 1:Ncell
        for kElem = 1:3
            if kElem == 1
                kernplot = squeeze(kernC{c}{p}(:,:,1,:)) + squeeze(kernC{c}{p}(:,:,2,:));
                squeeze(kernC{c}{p}(:,:,1,:));
            elseif kElem == 2
                kernplot = squeeze(kernC{c}{p}(:,:,1,:)); %OFF
            else
                kernplot = squeeze(kernC{c}{p}(:,:,2,:)); %ON
            end
            
            kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
            
            kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
            [ma idma] = max(kerndum,[],3); %find maxima from smoothed version
            [bestoriid bestposid] = find(ma == max(ma(:)));
            tau = idma(bestoriid,bestposid) + delayWinID(1) - 1;
            
            kernplot = smoothkern(kernplot,50,5,1,DM.taudom,DM.xdom);  %smooth heavily in spatial domain to get ori curve
            tcori = max(kernplot(:,:,tau),[],2);
            
            xp = TC.xpos{c}{kElem}(p);
            yp = TC.ypos{c}{kElem}(p);
            Amp = sqrt(xp^2 + yp^2);
            phase = atan2(TC.ypos{c}{kElem}(p),TC.xpos{c}{kElem}(p));
            OptPosID = Amp*cos(DM.oridom*pi/180 - phase);
            
            TC.tcoriall{c}{kElem}(p,:) = tcori';
            
            [TC.OMag{c}{kElem}(p) TC.OAng{c}{kElem}(p)] = orifind(tcori,DM.oridom);
            
            [param ffit varacc] = Gaussfit(DM.oridom,tcori',1);
            tcoriFit{c}{p} = exp(-((0:179)-param(1)).^2/(2*param(2)^2));
            
            param(1) = param(1)-DM.oridom(1);
            if param(1)<0
                param(1) = param(1)+180;
            end
            
            if varacc < .6
                param = param*NaN;
            end
            
            TC.OAng{c}{kElem}(p) = param(1);
            TC.OSig{c}{kElem}(p) = param(2);
            
            subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
            
            plot(tcori,colorstr(kElem))
            hold on
        end
        
    end
    
end

id = find(isnan(TC.xpos{1}{1}));
for kElem = 1:3
    TC.BWdiff{1}{kElem}(id) = NaN;
end

%% Plot ON and OFF scatter

xOn = TC.xpos{1}{3} - (nanmedian(TC.xpos{1}{3}) + nanmean(TC.xpos{1}{2}))/2;
xOff = TC.xpos{1}{2} - (nanmedian(TC.xpos{1}{3}) + nanmean(TC.xpos{1}{2}))/2;
yOn = TC.ypos{1}{3} - (nanmedian(TC.ypos{1}{3}) + nanmean(TC.ypos{1}{2}))/2;
yOff = TC.ypos{1}{2} - (nanmedian(TC.ypos{1}{3}) + nanmean(TC.ypos{1}{2}))/2;

xOn_mu = nanmean(xOn);
xOff_mu = nanmean(xOff);
yOn_mu = nanmean(yOn);
yOff_mu = nanmean(yOff);


figure, 
subplot(3,2,1)
plot(xOff,yOff,'ob')
hold on
plot(xOn,yOn,'or')
hold on

plot(xOff_mu,yOff_mu,'.b','MarkerSize',30)
hold on
plot(xOn_mu,yOn_mu,'.r','MarkerSize',30)
hold on

slope = (yOn_mu-yOff_mu)/(xOn_mu-xOff_mu);
yint = yOn_mu-slope*xOn_mu;

xran = [min([xOn xOff]) max([xOn xOff])];
yran = [min([yOn yOff]) max([yOn yOff])];

hold on
plot(xran,xran*slope + yint,'k')

xlim([-.8 .8]),ylim([-.8 .8])
axis square
legend('Off center','On center')

rot = atan(slope);

OnProj = xOn*cos(rot) + yOn*sin(rot);
OffProj = xOff*cos(rot) + yOff*sin(rot);

%OnProj = sqrt((xOn-xOff_mu).^2 + (yOn-yOff_mu).^2);
%OffProj = sqrt((xOff-xOff_mu).^2 + (yOff-yOff_mu).^2);


[dum bins] = hist([OnProj(:); OffProj(:)],15);
id = find(~isnan(OffProj));
hOff = hist(OffProj(id),bins);
id = find(~isnan(OnProj));
hOn = hist(OnProj(id),bins);

[H p] =ttest2(OnProj,OffProj)

subplot(3,2,3), bar(bins,hOn,'FaceColor',[1 0 0])
title(['p = ' num2str(p)])
xlabel('ON position (deg)')
subplot(3,2,5),bar(bins,hOff,'FaceColor',[0 0 1])
xlabel('OFF position (deg)')

dX = xOn-xOff;
dY = yOn-yOff;
subplot(3,2,2), 
xran = [min(dX) max(dX)];
yran = [min(dY) max(dY)]
plot([xran(1) xran(2)],[0 0],'k'), hold on, plot([0 0],[yran(1) yran(2)],'k')
plot(dX,dY,'ok')
xlabel('dipole dX (deg)')
ylabel('dipole dY (deg)')
grid on
axis image
%% Dipole axis vs. orientation preference vs. principal aspect ratio

dX = xOn-xOff;
dY = yOn-yOff;
dipoleSep = sqrt(dX.^2+dY.^2);
dipoleAxis = atan2(dY,dX)*180/pi-90;
id = find(dipoleAxis<0);
dipoleAxis(id) = dipoleAxis(id) + 180;
id = find(dipoleAxis>180);
dipoleAxis(id) = dipoleAxis(id) - 180;

id = find(dipoleSep>prctile(dipoleSep,0));
id2 = find(~isnan(dipoleSep));
id = unique([id id2]);

%id = [3 4 6]
figure, 
subplot(3,1,1),plot(TC.OAng{1}{1}(id),dipoleAxis(id),'o'), xlim([0 180]), ylim([0 180])
xlabel('ori preference'), ylabel('dipole axis + 90'), axis square
[r p] = circCorr(TC.OAng{1}{1}(id)*pi/180*2,dipoleAxis(id)*pi/180*2,2*pi,1);
title(['r = ' num2str(r) '; p = ' num2str(p)])

subplot(3,1,2), plot(TC.princAxis{1}{1}(id),dipoleAxis(id),'o'), xlim([0 180]), ylim([0 180])
xlabel('RF elongation axis'), ylabel('dipole axis + 90'), axis square
[r p] = circCorr(TC.princAxis{1}{1}(id)*pi/180*2,dipoleAxis(id)*pi/180*2,2*pi,2);
title(['r = ' num2str(r) '; p = ' num2str(p)])

subplot(3,1,3), plot(TC.princAxis{1}{1}(id),TC.OAng{1}{1}(id),'o'), xlim([0 180]), ylim([0 180])
xlabel('RF elongation axis'), ylabel('ori preference'), axis square
[r p] = circCorr(TC.princAxis{1}{1}(id)*pi/180*2,TC.OAng{1}{1}(id)*pi/180*2,2*pi,1);
title(['r = ' num2str(r) '; p = ' num2str(p)])

%%
kernElement = 1;
for c = 1:length(DM.colordom)
    figure
    for p = 1:Ncell

        subplot(1,3,1)
        %plot(mean(rawKern{c}{p}'),'k')
        plot(tcoriFit{c}{p},'k');
        xlabel('orientation')
        hold on

        subplot(1,3,2)
        dum = mean(TC.RFfit{c}{kernElement}{p});
        plot(xdom,dum/max(dum),'k')
        xlabel('x position')
        hold on

        subplot(1,3,3)
        dum = mean(TC.RFfit{c}{kernElement}{p}');
        plot(ydom,dum/max(dum),'k')
        xlabel('y position')
        hold on
        
    end
    
end


function smoother = getSmoother(tausig,orisig,possig,taudom,oridom,posdom)

ktau = exp(-(taudom-mean(taudom)).^2/(2*tausig^2)); %tausig is in ms

dom = posdom-mean(posdom);
kpos = exp(-dom.^2/(2*possig^2)); %possig is in degrees

oridomdum = linspace(0,180,length(oridom)+1);
oridomdum = oridomdum(1:end-1);  %I use this one in case oridom wraps around
kori = exp(-(oridomdum-oridomdum(ceil(end/2))).^2/(2*orisig^2)); %orisig is in degrees

kdum = kori'*kpos;
smoother = zeros(length(kori),length(kpos),length(ktau));
for i = 1:length(ktau)
    smoother(:,:,i) = kdum*ktau(i);
end
smoother = smoother/sum(smoother(:));


function kern = smoothkern(kern,tausig,orisig,possig,taudom,posdom)

Po = round(2*possig/(posdom(2)-posdom(1)));  %pad position by 2sig

dim = [size(kern,1)*2 size(kern,2)+2*Po size(kern,3)];

kpad = zeros(dim);
for i = 1:size(kern,3)
    kdum = kern(:,:,i);
    kdum = [kdum; fliplr(kdum)];
    kdum = [kdum(:,1)*ones(1,Po) kdum kdum(:,end)*ones(1,Po)];
    kpad(:,:,i) = kdum;    
end

%Make new domains to account for padded kernel
dim = size(kpad);
posdom = (1:dim(2))*(posdom(2)-posdom(1)); %remake posdom to account for padding

oridom = linspace(0,360,dim(1)+1);  %needs to go to 360
oridom = oridom(1:end-1);  %I use this one in case oridom wraps around

%Make Gaussian smoothers
ktau = exp(-(taudom-mean(taudom)).^2/(2*tausig^2)); %tausig is in ms
kpos = exp(-(posdom-mean(posdom)).^2/(2*possig^2)); %possig is in degrees
kori = exp(-(oridom-mean(oridom)).^2/(2*orisig^2)); %orisig is in degrees

kdum = kori'*kpos;
smoother = zeros(length(kori),length(kpos),length(ktau));
for i = 1:length(ktau)
    smoother(:,:,i) = kdum*ktau(i);
end
smoother = smoother/sum(smoother(:));

kern = ifftn(fftn(kpad).*abs(fftn(smoother)));

kern = kern(1:dim(1)/2,Po+1:dim(2)-Po,:);  %Truncate back to original size

function [OMag OAng] = orifind(G,oridomain)

R = sum(G'.*exp(1i*oridomain*pi/90));
OAng = angle(R);                 %-pi to pi
OAng = OAng + pi*(1-sign(OAng+eps));  %0 to 2pi
OAng = OAng*90/pi;               %0 to 180
OMag = abs(R)/(sum(G));


function [xpos ypos] = getmaxLoc(RF,D)

xdom = linspace(0,1,length(RF(1,:)))*getparam('x_size'); %degrees
ydom = linspace(0,1,length(RF(:,1)))*getparam('x_size');   %don't use 'y_size'
[xmat ymat] = meshgrid(xdom,ydom);

xdomI = linspace(xdom(1),xdom(end),D*length(xdom));
ydomI = linspace(ydom(1),ydom(end),D*length(ydom));
[xmatI ymatI] = meshgrid(xdomI,ydomI);
RFI = interp2(xmat,ymat,RF,xmatI,ymatI);

[yma xma] = find(RFI == max(RFI(:)));
xpos = xdomI(xma);
ypos = ydomI(yma);


function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;