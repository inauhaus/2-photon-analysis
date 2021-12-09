function [kmap_hor kmap_vert] = processkret_cellMask(f1)

global ACQinfo maskS

%f1 is a cell containing the result from 'f1meanimage'.  varargin is the optional
%filter kernel.  Each of the images are smoothed with a Gaussian with a std 
%dev of 'stdev' and a width of 'width'.

ang{2} = f1{2}; %for one axis
ang{4} = f1{4}; 
ang{1} = f1{1};
ang{3} = f1{3};


% ang0 = ang0*exp(-j*180*pi/180);
% ang1 = ang1*exp(-j*180*pi/180);
% ang2 = ang2*exp(-j*180*pi/180);
% ang3 = ang3*exp(-j*180*pi/180);

%The negative is to show where it peaks in the range of -180 to 180.
%i.e. -180 is the left most side of the stimulus.  Without the negative,
%an angle of -180 would have been the middle of the stimulus.
%angle(FourierTX(cos(wt-0))) == 0

ang{2} = angle(-ang{2});
ang{4} = angle(-ang{4});
ang{1} = angle(-ang{1});
ang{3} = angle(-ang{3});

%% Plot

sPerx = 360;
sPery = 360;

[xmicperpix ymicperpix] = getImResolution;

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

maskdum = maskS.bwCell{1};

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
cdom = jet;
for i = 1:length(ang)    
    angmap{i} = zeros(size(maskS.bwCell{1}));
    for p = 1:length(celldom)        
        cellidx = find(masklabel == celldom(p));
        angmap{i}(cellidx) = ang{i}(p);
        if isnan(ang{i}(p))
           maskdum(cellidx) = 0;
        end
    end
end

figure,
subplot(2,2,1),imagesc(angmap{2}*180/pi,'AlphaData',maskdum,[-180 180]), colorbar, colormap hsv, title('90')
axis image
%hold on, contour(ang1,'k')
subplot(2,2,2),imagesc(angmap{4}*180/pi,'AlphaData',maskdum,[-180 180]), colorbar, colormap hsv, title('270')
axis image
%hold on,contour(ang3,'k')

subplot(2,2,3),imagesc(angmap{3}*180/pi,'AlphaData',maskdum,[-180 180]), colorbar, colormap hsv, title('180')
axis image
%hold on, contour(ang2,'k')
subplot(2,2,4),imagesc(angmap{1}*180/pi,'AlphaData',maskdum,[-180 180]), colorbar, colormap hsv, title('0')
axis image



%%

unwrapflag = 0;
if unwrapflag
    %%Find delay as the sum of the 2 vectors
    delay_hor = angle(exp(1i*ang{1}) + exp(1i*ang{3})); %Will be 0 or 180 if delay is 0
    delay_vert = angle(exp(1i*ang{2}) + exp(1i*ang{4}));
    
    %Make delay go from 0 to pi and 0 to pi, instead of 0 to pi and 0 to -pi.
    %The delay can't be negative.  If the delay vector is in the bottom two
    %quadrants, it is assumed that the it started at 180.  The delay always
    %pushes the vectors counter clockwise.
    delay_hor = delay_hor + pi/2*(1-sign(delay_hor));
    delay_vert = delay_vert + pi/2*(1-sign(delay_vert));
    
    %Use delay vector to calculate retinotopy.
    kmap_hor = .5*(angle(exp(1i*(ang{1}-delay_hor))) - angle(exp(1i*(ang{3}-delay_hor))));
    kmap_vert = .5*(angle(exp(1i*(ang{2}-delay_vert))) - angle(exp(1i*(ang{4}-delay_vert))));
    
else
    
%  
    delay_hor = angle(exp(1i*ang{1}) + exp(1i*ang{3})); %Will be 0 or 180 if delay is 0
    delay_vert = angle(exp(1i*ang{2}) + exp(1i*ang{4}));
    
    delay_hor = angle(exp(1i*delay_hor*2))/2;
    delay_vert = angle(exp(1i*delay_vert*2))/2;

%    delay_hor = delay_hor + pi*(1-sign(delay_hor)); %add 360 to negative values
%    delay_vert = delay_vert + pi*(1-sign(delay_vert));
% 
%    
% %     delay_hor = angle(exp(1i*delay_hor) + exp(1i*delay_vert));
% %     delay_vert = delay_hor;
%     
%     kmap_hor = .5*(exp(1i*ang{1}) + exp(-1i*ang{3}));       
%     id = find(delay_hor < pi/2);   
%     kmap_hor(id) = angle(kmap_hor(id));
%     id = find(delay_hor > pi/2);
%     kmap_hor(id) = angle(kmap_hor(id)*exp(1i*pi));
%     
%     %This is theoretically correct, but screws up sometimes
%     kmap_vert = .5*(exp(1i*ang{2}) + exp(-1i*ang{4}));
%     id = find(delay_vert < pi/2);   
%     kmap_vert(id) = angle(kmap_vert(id));
%     id = find(delay_vert > pi/2);
%     kmap_vert(id) = angle(kmap_vert(id)*exp(1i*pi));
%    
    %%The code below is much simpler, but doesn't account for delays longer than pi/2:
    kmap_hor = .5*(exp(1i*ang{1}) + exp(-1i*ang{3}));
    kmap_vert = .5*(exp(1i*ang{2}) + exp(-1i*ang{4}));
    kmap_hor = angle(kmap_hor);
    kmap_vert = angle(kmap_vert);

    %Get rid of cells with a delay longer than pi/2.
    idbad = find((delay_vert)<-pi/8  |  (delay_vert)>pi/4 );
    kmap_vert(idbad) = NaN;
    idbad = find((delay_hor)<-pi/8 | (delay_hor)>pi/4);
    kmap_hor(idbad) = NaN;
    
    %d = angle(median(exp(1i*kmap_hor))-exp(1i*kmap_hor));

    
end

kmap_vert = real(kmap_vert);
kmap_hor = real(kmap_hor);

%radians to degrees
%delay_hor = delay_hor*180/pi.*bw;
kmap_hor = kmap_hor*180/pi;
%delay_vert = delay_vert*180/pi.*bw;
kmap_vert = kmap_vert*180/pi;


%Create shadow of ROI coverage.
x = floor(100/360*(kmap_hor+180))+1;
y = floor(100/360*(-kmap_vert+180))+1;
%sh = shadow(x,y,100,100);
sh = [];


kmap_hor = kmap_hor/360*getparam('x_size');
kmap_vert = kmap_vert/360*getparam('y_size');

figure,hist(kmap_hor)
