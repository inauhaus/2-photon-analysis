function cellMat = getGeoTrxTimeCourse2(CH,Px,Py)

%Get cell time courses
%2 uses the original mask, and smooths it with a shifted (by Px/Py)
%Gaussian to create the new mask

global ACQinfo maskS

global CoMx2 CoMy2

M = ACQinfo.linesPerFrame;
N = ACQinfo.pixelsPerLine;

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);

Ncell = length(celldom);
Nframe = length(CH(1,1,:));

[xmicperpix ymicperpix] = getImResolution;

micperpix = (xmicperpix + ymicperpix)/2;  %used to define the size of the mask

sp = ACQinfo.linesPerFrame*ACQinfo.msPerLine;

%Make a smoothing kernel to better determine local maxima
smoother = fspecial('gaussian',[M N],1);
smoother = abs(fft2(smoother));

tempfilt = ifft2(fft2(maskS.im{1}).*smoother); %filtered template

m = 1:M;
n = 1:N;

CoMy = zeros(1,Ncell);
CoMx = zeros(1,Ncell);
for p = 1:Ncell
    [idy0{p} idx0{p}] = find(masklabel == celldom(p));
    CoMy(p) = mean(idy0{p});
    CoMx(p) = mean(idx0{p});
    
    idCvec{p} = (round(idx0{p})-1)*M + round(idy0{p}); 
    cbright(p) = mean(maskS.im{1}(idCvec{p}));
end
cbright = cbright/sum(cbright);
%%

%Preallocate
dxm = zeros(Nframe,Ncell); dym = zeros(Nframe,Ncell);
dLocmaxX = zeros(Nframe,Ncell); dLocmaxY = zeros(Nframe,Ncell);
CoMx2 = cell(1,Nframe); CoMy2 = cell(1,Nframe);
cellMat = zeros(Ncell,Nframe);

ma = max(max(max(CH(10:end-10,10:end-10,2:end-1),[],1),[],2),[],3);
mi = min(min(min(CH(10:end-10,10:end-10,2:end-1),[],1),[],2),[],3);
nI = 1:.1:N;

%First Wcell is only to find the local maxima, so it should be smaller than
%the one set for the next loop
Wcell = 7; %must be 'odd'

%This first loop gets the cell location and difference from the local max,
%and stores it for the next loop
for f = 1:Nframe
    
    %Smoothed frame to better determine local maxima
    %tensdum = ifft2(fft2(CH(:,:,f)).*smoother); %needed if I compute dLocmaxX
    
    CoMx2{f} = zeros(1,Ncell);
    CoMy2{f} = zeros(1,Ncell);
    idframe = ((f-1)*M+1):f*M;      
    yI = 1:.1:M;
    xI = 1:.1:M;
    PyI = interp1(1:M,Py(idframe),yI);
    PxI = interp1(1:M,Px(idframe),xI);             

    for p = 1:Ncell               
        
        %Get the row and column of the cell's new location
        %-Py tells me the new position of the entire brain, at each column, in pixels.
        %So, -Py+yI gives the new position of each column.  Therefore, we
        %are finding the index that most closely matches CoMy in order to
        %find its new location.
        [dum idy] = min(abs(-PyI+yI - CoMy(p))); %get y index of cell location within the interpolated domain, after shifting by PyI
        [dum idx] = min(abs(-PxI(idy)+nI - CoMx(p))); %get x index of cell location within the interpolated domain, after shifting by PxI
        CoMy2{f}(p) = yI(idy); %new y location for center of mass... I save this only so that I can plot the dots later
        CoMx2{f}(p) = nI(idx); %new x location for center of mass
        
        %dy = CoMy2{f}(p)-CoMy(p);  %CoM Difference from template location
        %dx = CoMx2{f}(p)-CoMx(p);  
        
        
        if CoMy2{f}(p)>=1 & CoMy2{f}(p)<=M & CoMx2{f}(p)>=1 & CoMx2{f}(p)<=N     %if it is in the ROI still
            
            %residual location 
            relid = 1;  %This can be set to 1 or f.  Setting it to 1 makes it easier to visualize how well things are working in the movie... 
            %But if slow movement is really big in the trial, things might get a bit messed up.                      
            
            dxm(f,p) = CoMx2{f}(p) - CoMx(p);  %the shift for the Gaussian in the next loop
            dym(f,p) = CoMy2{f}(p) - CoMy(p);                        
   
            
        else
            dxm(f,p) = NaN;
            dym(f,p) = NaN;
            CoMx2{f}(p) = NaN;
            CoMy2{f}(p) = NaN;
            
            if f == 1  %if I don't do this it will screw up the rest of the trial for this cell
                CoMx2{f}(p) = CoMx(p);
                CoMy2{f}(p) = CoMy(p);
            end
        end
    end
    
            
end

%% 


sig = 3/micperpix;  %2.505 is a heuristically determined constant.
Wcell = round(sig*10);  %width of window to cell
Wcell = Wcell + (1-rem(Wcell,2)); %must be 'odd', 
[mdomx1 mdomy1] = meshgrid(ceil(-Wcell/2):floor(Wcell/2),ceil(-Wcell/2):floor(Wcell/2));  %domain before shift

global masksig
masksig = 2;

%This loop takes the cell location from the previous loops to create the
%actual time course.
clear F
maskdum = zeros(size(maskS.bwCell{1}));
for p = 1:Ncell
    
    maskdum = maskdum*0;
    maskdum(idCvec{p}) = 1;
    
    %maskdum(round(CoMy(p)),round(CoMx(p))) = 1; %if we just want a Gaussian
    
    for f = 1:length(CH(1,1,:))
        tensdum = CH(:,:,f);   
        
        if ~isnan(dxm(f,p))
            %make the mask for this cell and image:
            mdomx = mdomx1-dxm(f,p); mdomy = mdomy1 - dym(f,p); 
            rdom = sqrt(mdomx.^2 + mdomy.^2);
            Gs = exp(-rdom.^2/(2*masksig^2));
%             id = find(rdom(:) > 2*masksig^2); %Truncate at sig^2
%             Gs(id) = 0; Gs = Gs/sum(Gs(:));           
            
            %Cut out a piece around the cell, approximately centered on the cell, the dxm shift will take care of the rest:
            yran = (round(CoMy(p))-floor(Wcell/2)):(round(CoMy(p))+floor(Wcell/2));
            xran = (round(CoMx(p))-floor(Wcell/2)):(round(CoMx(p))+floor(Wcell/2));
            
            id = find(xran<1 | xran > N);
            xran(id) = []; Gs(:,id) = [];
            id = find(yran<1 | yran > M);
            yran(id) = []; Gs(id,:) = [];
            
            if ~isnan(CoMy2{relid}(p))
                impiece = tensdum(yran,xran);
                
                maskpiece = maskdum(yran,xran);
                mask = conv2(maskpiece,Gs,'same');
                mask = mask/sum(mask(:));
                
                %Dot product between mask and image:
                cellMat(p,f) = mean(impiece(:).*mask(:));                
            else
                cellMat(p,f) = NaN;
            end
% 
% 
%             figure(9)
%               clf
%             subplot(1,2,1), 
%             imagesc(impiece,[0 5000]), colormap gray
%             hold on
%            contour(mask,exp(-masksig^2/2)/max(mask(:)),'r')
%            contour(mask,1,'r')
%            axis image
%             subplot(1,2,2), imagesc(mask), colormap gray
%             axis image
%             pause(.05)
%             drawnow
%             F(f) = getframe;
        else
            cellMat(p,f) = NaN;
        end        
        
                
    end

 
    % Filter out the residual motion component
 
%     signal = cellMat(p,:)';
%     signalNaN = find(isnan(signal));
%     signal(signalNaN) = nanmedian(signal);
%     
%     refx = dxm(:,p);   
%     refy = dym(:,p); 
%     
%     %First get power spectra
%     W = 20;
%     W = round(1000*W/sp)-1;    
%     [refxpow fdom] = getPowerSpect(refx,sp,W,round(W/2));
%     [refypow fdom] = getPowerSpect(refy,sp,W,round(W/2));
%     [sigpow fdom] = getPowerSpect(signal,sp,W,round(W/2));
%     Pow = refypow+refxpow;
%     
%     %Adaptive removal of the peaks
%     W = 20;
%     W = round(1000*W/sp)-1;
%     idbreath = find(fdom<2.5 & fdom>.5);  %range of breathing and heart rate
%     for q = 1:2
% 
%         [dum id] = max(abs(Pow(idbreath)));
%         idpeak = id+idbreath(1)-1;        
% 
%         signal = adaptiveLineFilt_ref(signal,refx,refy,fdom(idpeak),sp,W,round(W/2));
%         
%         %Pow(idpeak-2:idpeak+2) = 0;
% 
%     end
%     
%     signal(signalNaN) = NaN;
%     cellMat(p,:) = signal;
% 
end

cellMat = cellMat(:,1:end-1); %Get rid of the last frame

%%     
clear F
figure(21)
for f = 1:length(CH(1,1,:))
    tensdum = CH(:,:,f);
    clf
    imagesc(tensdum,[mi ma]), colormap gray
    hold on
    plot(CoMx2{f},CoMy2{f},'.r')
    
    axis image
%     xlim([32 75])
%     ylim([1 60])    
    axis off
    drawnow
    % F(f) = getframe;
    
    %pause(.001)
end

% figure(88)
% fid = [10 round(length(CH(1,1,:))/2) length(CH(1,1,:))-10];
% for i = 1:length(fid)
%     subplot(1,length(fid),i)
%     imagesc(CH(:,:,fid(i)),[mi ma]), colormap gray
%     hold on
%     plot(CoMx2{fid(i)},CoMy2{fid(i)},'.r')
% end

%%  plot traces

% tdom = (0:(length(dym)-1))*sp;
% id = 20:(length(tdom)-20);
% figure,
% subplot(2,1,1)
% plot(tdom(1:length(id))/1000,dxm(id,55)*xmicperpix)
% hold on
% plot(tdom(1:length(id))/1000,dym(id,55)*ymicperpix,'r')
% legend('xmotion','ymotion')
% ylabel('microns')
% xlabel('sec')
% xlim([0 60])
% 
% subplot(2,1,2)
% plot(fdom,refxpow)
% hold on
% plot(fdom,refypow,'r')
% xlabel('frequency (Hz)')
% xlim([0 fdom(end)])
