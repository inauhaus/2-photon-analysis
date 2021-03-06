function cellMat = getGeoTrxTimeCourse2(CH,Px,Py)

%This older one is the same, but it was overcomplicated.

%Get cell time courses
%2 uses the original mask, and smooths it with a shifted (by Px/Py)
%Gaussian to create the new mask

global ACQinfo maskS

M = ACQinfo.linesPerFrame;
N = ACQinfo.pixelsPerLine;

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);

Ncell = length(celldom);
Nframe = length(CH(1,1,:));

xVolt = ACQinfo.scanAmplitudeX/ACQinfo.zoomFactor;
yVolt = ACQinfo.scanAmplitudeY/ACQinfo.zoomFactor;
xmic = 94*xVolt-2;  %Kristina fit these lines
ymic = 135*yVolt+0.5;
xmicperpix = xmic/ACQinfo.pixelsPerLine;
ymicperpix = ymic/ACQinfo.linesPerFrame;
micperpix = (xmicperpix + ymicperpix)/2;  %used to define the size of the mask

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
        [dum idy] = min(abs(PyI+yI - CoMy(p))); %get y index of cell location within the interpolated domain, after shifting by PyI
        [dum idx] = min(abs(PxI(idy)+nI - CoMx(p))); %get x index of cell location within the interpolated domain, after shifting by PxI
        idy = yI(idy); %new y location for center of mass
        idx = nI(idx); %new x location for center of mass
        
        dy = idy-CoMy(p);  %CoM Difference from template location
        dx = idx-CoMx(p);
         
        idcelly = idy0{p} - dy;  %New cell location (vector of pixels from mask)
        idcellx = idx0{p} - dx;        
        idbad = find(idcelly > M | idcelly < 1 | idcellx > N | idcellx < 1);
        idcelly(idbad) = [];
        idcellx(idbad) = [];
        
        idvec = (round(idcellx)-1)*M + round(idcelly);  %New cell location (vector of pixels from mask) vectorized location         
        idvec = round(idvec);
        
        if ~isempty(idvec)    
            
           CoMy2{f}(p) = mean(idcelly);  %New cell location, center of mass  (could have also used idy/idx)
           CoMx2{f}(p) = mean(idcellx);
            
            CoMy3{f}(p) = idy;
            CoMx3{f}(p) = idx;
            
            %residual location 
            relid = 1;  %This can be set to 1 or f.  Setting it to 1 makes it easier to visualize how well things are working in the movie... 
            %But if slow movement is really big in the trial, things might get a bit messed up.
            
            %CoMx2 is used to get a window centered on the cell.  dxm/dxy
            %will be used to refine the mask on the subpixel scale
            dxm(f,p) = CoMx2{f}(p) - round(CoMx2{relid}(p));
            dym(f,p) = CoMy2{f}(p) - round(CoMy2{relid}(p));
            
            %Get difference between local maxima and estimate of cell
            %location            
%             yran = (round(CoMy2{relid}(p))-floor(Wcell/2)):(round(CoMy2{relid}(p))+floor(Wcell/2));
%             xran = (round(CoMx2{relid}(p))-floor(Wcell/2)):(round(CoMx2{relid}(p))+floor(Wcell/2));
            
%             id = find(xran<1 | xran > N);
%             xran(id) = []; 
%             id = find(yran<1 | yran > M);
%             yran(id) = [];            

%             impiece = tensdum(yran,xran);
%             [lY lX] = find(impiece == max(impiece(:)));
%             dLocmaxX(f,p) = lX(1)-dxm(f,p)-ceil(Wcell/2);  %Why do I need this?
%             dLocmaxY(f,p) = lY(1)-dym(f,p)-ceil(Wcell/2);            
            
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
    
   %CoMx2{f} = CoMx2{f} + sum(dLocmaxX(f,:).*cbright);
   %CoMy2{f} = CoMy2{f} + sum(dLocmaxY(f,:).*cbright);
            
end

% dumx = dLocmaxX.*(ones(Nframe,1)*cbright)/Nframe; dumx = sum(dumx(:));
% dumy = dLocmaxY.*(ones(Nframe,1)*cbright)/Nframe; dumy = sum(dumy(:));
% for f = 1:Nframe
%     CoMx2{f} = CoMx2{f} + dumx;
%     CoMy2{f} = CoMy2{f} + dumy;
%     
%     dxm(f,:) = dxm(f,:) + dumx;
%     dym(f,:) = dym(f,:) + dumy;
% end
%%

%dxm = dLocmaxX;
%dym = dLocmaxY;
 
sig = 3/micperpix;  %2.505 is a heuristically determined constant.
Wcell = round(sig*8);  %width of window to cell
Wcell = Wcell + (1-rem(Wcell,2)); %must be 'odd', 
[mdomx1 mdomy1] = meshgrid(ceil(-Wcell/2):floor(Wcell/2),ceil(-Wcell/2):floor(Wcell/2));  %domain before shift

sig = 1;

%This loop takes the cell location from the previous loops to create the
%actual time course.
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
            Gs = exp(-rdom.^2/(2*sig^2));
            %id = find(rdom(:) > 2*sig^2); %Truncate at sig^2
            %Gs(id) = 0; Gs = Gs/sum(Gs(:));           
            
            %Cut out a piece around the cell, approximately centered on the cell, the dxm shift will take care of the rest:
            yran = (round(CoMy2{relid}(p))-floor(Wcell/2)):(round(CoMy2{relid}(p))+floor(Wcell/2));
            xran = (round(CoMx2{relid}(p))-floor(Wcell/2)):(round(CoMx2{relid}(p))+floor(Wcell/2));
            
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
%             subplot(1,2,1), imagesc(impiece,[0 5000]), colormap gray
%             hold on
%            % contour(mask,exp(-sig^2/2)/max(mask(:)),'r')
%            contour(mask,1,'r')
%             subplot(1,2,2), imagesc(mask), colormap gray
%             pause(.05)
%             drawnow
        else
            cellMat(p,f) = NaN;
        end        
        
                
    end


end

cellMat = cellMat(:,1:end-1); %Get rid of the last frame

%%     
figure(21)
for f = 1:length(CH(1,1,:))
    tensdum = CH(:,:,f);
    clf
    imagesc(tensdum,[mi ma]), colormap gray
    hold on
    plot(CoMx2{f},CoMy2{f},'.r')
        hold on
    plot(CoMx3{f},CoMy3{f},'.g')
    drawnow
    pause(.08)
end

figure(88)
fid = [10 round(length(CH(1,1,:))/2) length(CH(1,1,:))-10];
for i = 1:length(fid)
    subplot(1,length(fid),i)
    imagesc(CH(:,:,fid(i)),[mi ma]), colormap gray
    hold on
    plot(CoMx2{fid(i)},CoMy2{fid(i)},'.r')
end
