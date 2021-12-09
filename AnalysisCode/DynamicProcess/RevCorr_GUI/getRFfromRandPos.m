function getRFfromRandPos

global ACQinfo G_RChandles cellS DM TC

%%%%

%Downsample kernel and domains
if ~isempty(find(isnan(cellS.kernAll{1})))
    for p = 1:length(cellS.kernAll)
        
        dumA = cellS.kernAll{p}(:,1:2:end,:,:,:,:); dumB = cellS.kernAll{p}(:,2:2:end,:,:,:,:);
        id = find(isnan(dumA)); dumA(id) = dumB(id);
        id = find(isnan(dumB)); dumB(id) = dumA(id);

        cellS.kernAll{p} = (dumA + dumB)/2; %oridomain        

        dumA = cellS.kernAll{p}(1:2:end,:,:,:,:,:); dumB = cellS.kernAll{p}(2:2:end,:,:,:,:,:);
        id = find(isnan(dumA)); dumA(id) = dumB(id);
        id = find(isnan(dumB)); dumB(id) = dumA(id);

        cellS.kernAll{p} = (dumA + dumB)/2; %spatial domain
        
    end
    
    DM.oridom = (DM.oridom(1:2:end) + DM.oridom(2:2:end))/2;
    DM.xdom = (DM.xdom(1:2:end) + DM.xdom(2:2:end))/2;

    'Downsampling because not enough presentations'
end



for c = 1:length(DM.colordom)
    for p = 1:length(cellS.kernAll)
        kernC{c}{p} = squeeze(cellS.kernAll{p}(:,:,:,c,:));
    end    
end


ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with acqPeriod spacing, and end at an estimate of kernDel(2)

maDel = getparam('h_per')*10 + 150;
delayWin = [100 maDel];
delayWin = [200 500];  %assume the peak response is within this time window
[dum delayWinID(1)] = min(abs(delayWin(1)-taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-taudom));

Ncell = length(cellS.kernAll);

%Get the receptive field

kernsmooth = getSmoother(70,15,.3,DM.taudom,DM.oridom,DM.xdom); %used to find maxima

param(5) = -36;

for c = 1:length(DM.colordom)
    TC.xpos{c} = NaN*ones(1,Ncell);
    TC.ypos{c} = NaN*ones(1,Ncell);
    
    figure
    for p = 1:Ncell

        subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

        kernplot_b = squeeze(kernC{c}{p}(:,:,1,:));
        kernplot_w = squeeze(kernC{c}{p}(:,:,2,:));

        kerndum_b = ifftn(fftn(kernplot_b).*abs(fftn(kernsmooth)));
        kerndum_b = kerndum_b(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        kerndum_w = ifftn(fftn(kernplot_w).*abs(fftn(kernsmooth)));
        kerndum_w = kerndum_w(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window

        kernplot_b = smoothkern(kernplot_b,50,5,.1,DM.taudom,DM.xdom);  %this will append to account for wrapping, then smooth
        kernplot_w = smoothkern(kernplot_w,50,5,.1,DM.taudom,DM.xdom);

        [ma idma] = max(kerndum_b+kerndum_w,[],3); %find maxima from smoothed version
        [bestoriid bestposid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestposid) + delayWinID(1) - 1;

        im_b = kernplot_b(:,:,tau);
        im_w = kernplot_w(:,:,tau);
        
%         im_b = phi(im_b-prctile(im_b(:),50));
%         im_w = phi(im_w-prctile(im_w(:),50));

        degperpix = getparam('x_size')/(length(im_b(1,:))-1);  %can't use xdom because of the padding

        %     for i = 1:length(DM.oridom)
        %         tc = im_b(i,:) + im_b(i,:);
        %         [param ffit varaccount] = Gaussfit(1:length(tc),tc,0);
        %         imfit(:,i) = ffit;
        %     end

        imfit = im_w+im_b;
        imfit = phi(imfit-prctile(imfit(:),20));
        padder = zeros(length(imfit(:,1)),0);
        imfit = [padder imfit padder]';
        
        %imfit(:,end/2:end) = 0;
        
        %%%%%%%Control%%%%%%%%%
%         imfit = [imfit flipud(imfit)];
%         shiftori = round(rand(1)*(size(imfit,2)-1));
%         imfit = circshift(imfit,[0 shiftori]);
%         imfit = imfit(:,1:size(imfit,2)/2);
        %%%%%%%%%%%%%%%%%%%%%%%%

        RF = iradon(imfit,DM.oridom,'none','linear',1,length(imfit(:,1)));
        %RF = backproj(imfit,DM.oridom);

        hh = zeros(size(RF));
        hh(1:5,1:5) = [.1 .3 1 .3 .1]'*[.1 .3 1 .3 .1];
        %RF = ifft2(fft2(RF).*abs(fft2(hh)));
        
        thresh = prctile(RF(:),20) + (max(RF(:))-prctile(RF(:),20))*.3;
        RF = phi(RF-thresh);

        %     [TC.xpos{c}(p) TC.ypos{c}(p)] = getmaxLoc(RF,30);  %returns 'degrees'

        otc = max(imfit);

        [dum ohat] = orifind(otc',DM.oridom);

        [param ffit varaccount] = Gaussfit2Drot(RF,ohat);
        
        if varaccount > 0.6
            TC.xpos{c}(p) = param(2)*degperpix;
            TC.ypos{c}(p) = param(1)*degperpix;
            TC.xsize{c}(p) = param(3)*degperpix*2;
            TC.ysize{c}(p) = param(4)*degperpix*2;
        else
            TC.xpos{c}(p) = NaN;
            TC.ypos{c}(p) = NaN;
            TC.xsize{c}(p) = NaN; 
            TC.ysize{c}(p) = NaN;
        end

        if TC.ysize{c}(p)>TC.xsize{c}(p)
            param(5) = param(5)+90;
        end
        if param(5) < 0
            param(5) = param(5)+180;
        end
        if param(5) > 180
            param(5) = param(5)-180;
        end
        Opref{c}(p) = ohat;
        Opreffit{c}(p) = param(5);
        Omagfit{c}(p) = max([TC.xsize{c}(p) TC.ysize{c}(p)])/min([TC.xsize{c}(p) TC.ysize{c}(p)]);  %aspect ratio
        
 

        %     [yma xma] = find(RF == max(RF(:)));
        %     xdom = linspace(0,1,length(RF(1,:)))*getparam('x_size');
        %     ydom = linspace(0,1,length(RF(:,1)))*getparam('x_size');   %don't use 'y_size'
        %     xmarg = mean(RF);
        %     ymarg = mean(RF,2)';
        %     xmarg = RF(yma,:);
        %     ymarg = RF(:,xma)';
        %     xmarg = xmarg/sum(xmarg);
        %     ymarg = ymarg/sum(ymarg);
        %     TC.xpos{c}(p) = sum(xmarg.*xdom);
        %     TC.ypos{c}(p) = sum(ymarg.*ydom);


        xdom = (0:length(RF(1,:))-1)*degperpix; %degrees
        ydom = (0:length(RF(:,1))-1)*degperpix;   %don't use 'y_size'
        RF = RF/max(RF(:));
        if p == 1
            [yma xma] = find(RF == max(RF(:)));
            yran = yma-5:yma+5;
            xran = xma-5:xma+5;
            yran = 1:length(ydom);
            xran = 1:length(xdom);
        end
        dim = size(RF);
        %RF = RF(yran,xran); ffit = ffit(yran,xran); ydom = ydom(yran); xdom = xdom(xran);
        
        TC.RF{c}{p} = RF;

        imagesc(xdom,ydom,RF), %title(num2str(param(5)))
        hold on
        contour(xdom,ydom,ffit,.6,'k')
        axis image
        drawnow

        [idy idx] = find(imfit' == max(imfit(:)));

        %plot([im_b(idy,:)' im_w(idy,:)'])
        
%        imagesc([imfit flipud(imfit)])
%         drawnow
%         axis off

        if p == 1
            title([num2str(ydom(end)-ydom(1)) 'degrees'])
        end

        %     contour(xdom,ydom,RF,.5,'k'), axis ij
        %     hold on
        % %     ylim([1.5 3.8])
        % %     xlim([.2 2.3])
        %     drawnow

        %plot([mean(kernplot_b(4:5,:))' mean(kernplot_w(4:5,:))'])
        
        

    end
end

% thet = 40;
% xposp = TC.xpos*cos(thet) + TC.ypos*sin(thet);
% yposp = -TC.xpos*sin(thet) + TC.ypos*cos(thet);

%Get and plot ori tuning 
    
figure
for c = 1:length(DM.colordom)
    for p = 1:Ncell
    
        kernplot_b = squeeze(kernC{c}{p}(:,:,1,:));
        kernplot_w = squeeze(kernC{c}{p}(:,:,2,:));

        kerndum_b = ifftn(fftn(kernplot_b).*abs(fftn(kernsmooth)));
        kerndum_b = kerndum_b(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        kerndum_w = ifftn(fftn(kernplot_w).*abs(fftn(kernsmooth)));
        kerndum_w = kerndum_w(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        
        kernplot_b = smoothkern(kernplot_b,50,10,.5,DM.taudom,DM.xdom);  %this will append to accound for wrapping, then smooth
        kernplot_w = smoothkern(kernplot_w,50,10,.3,DM.taudom,DM.xdom);

        [ma idma] = max(kerndum_b+kerndum_w,[],3); %find maxima from smoothed version
        [bestoriid bestposid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestposid) + delayWinID(1) - 1;

        %     tcori_b = kernplot_b(:,bestposid,tau);  %Can't do this dumb assssss!!!
        %     tcori_w = kernplot_w(:,bestposid,tau);

        %     tcori_b = mean(kernplot_b(:,:,tau),2);  %Shouldn't do this either. If it is actually a projection then they will all be equal
        %     tcori_w = mean(kernplot_w(:,:,tau),2);
        %     tcori = tcori_b+tcori_w;

        tcori = max(kernplot_b(:,:,tau)+kernplot_w(:,:,tau),[],2);
        
        TC.tcoriall{c}(p,:) = tcori;

        [TC.OMag{c}(p) TC.OAng{c}(p)] = orifind(tcori,DM.oridom);

        [param ffit varacc] = Gaussfit(DM.oridom,tcori',1);
        param(1) = param(1)-DM.oridom(1);
        if param(1)<0
            param(1) = param(1)+180;
        end

        if varacc < .6
            param = param*NaN;
        end

        TC.OAng{c}(p) = param(1);
        TC.OSig{c}(p) = param(2);

        subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

        plot(tcori)

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