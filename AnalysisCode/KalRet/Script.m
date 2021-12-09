%Script for analyzing kalatsky retinotopy with 2photon data

%f1 = f1meanimage;  %Build F1 images (takes the longest)
global Tens Analyzer
%First set dirs and hit 'Compute F0 images' within pF0
f1 = CondF1;
L = fspecial('gaussian',30,1);  %make spatial filter
bw = ones(size(f1{1}));

%%
h = fspecial('gaussian',size(f1{1}),3);

[kmap_hor kmap_vert] = processkret(f1,bw,h);  %Make maps to plot, delete L if no smoothing

projectorAdjustment = 1;

if projectorAdjustment
    
    kmap_horx = kmap_hor;
    kmap_vertx = kmap_vert;
    kmap_hor = kmap_vertx;
    kmap_vert = kmap_horx; %This should be negative if run with LCD rotated 90 clockwise
    
end

[kmap_vert kmap_hor] = unwrapKmap(kmap_vert,kmap_hor);

%%

figure
sPerx = getparam('x_size');
sPery = getparam('y_size');

imagesc(kmap_hor,[-sPery/2 sPery/2])
title('Horizontal Retinotopy')
%colormap jet
colorbar
axis image

figure
imagesc(kmap_vert,[-sPerx/2 sPerx/2])
title('Vertical Retinotopy')
colorbar
%colormap jet
axis image

%bw = roipoly;
%% Convert to screen coordinates

%Correction to make zero x_pos y_pos, not the perpendicular bisector,  the
%zero point
dxperp_cm = getparam('dx_perpbis');
dyperp_cm = getparam('dy_perpbis');
dxperp_deg = atan(dxperp_cm/Analyzer.M.screenDist)*180/pi; %Shift in deg of visual field
dyperp_deg = atan(dyperp_cm/Analyzer.M.screenDist)*180/pi; %Shift in deg of visual field
kmap_hor_re = kmap_hor+dyperp_deg; 
kmap_vert_re = kmap_vert+dxperp_deg;

kmap_vert_pix = tan(kmap_vert_re*pi/180)*Analyzer.M.screenDist*Analyzer.M.xpixels/Analyzer.M.screenXcm;
kmap_hor_pix = tan(kmap_hor_re*pi/180)*Analyzer.M.screenDist*Analyzer.M.ypixels/Analyzer.M.screenYcm;

xPer_pix = tan(getparam('x_size')*pi/180/2)*Analyzer.M.screenDist*Analyzer.M.xpixels/Analyzer.M.screenXcm*2;
yPer_pix = tan(getparam('y_size')*pi/180/2)*Analyzer.M.screenDist*Analyzer.M.ypixels/Analyzer.M.screenYcm*2;

ymid = getparam('y_pos')*0;
xmid = getparam('x_pos')*0;
kmap_hor_pix = kmap_hor_pix+ymid;
kmap_vert_pix = kmap_vert_pix+xmid;

yrange = [-yPer_pix/2 yPer_pix/2]+ymid; %beginning and ending of stimulus, in pixels
xrange = [-xPer_pix/2 xPer_pix/2]+xmid;

figure

subplot(2,2,1)
imagesc(kmap_hor_pix.*bw,yrange)
title('Horizontal Retinotopy (pixels)')
colormap jet
colorbar
axis image

subplot(2,2,2)
imagesc(kmap_vert_pix.*bw,xrange)
title('Vertical Retinotopy (pixels)')
colorbar
colormap jet
axis image

horsub = kmap_hor_pix(find(bw(:)));
id = find(horsub<prctile(horsub,2) | horsub>prctile(horsub,98));
horsub(id) = []; %Remove outliers

subplot(2,2,3)
hist(horsub,40);
title(['median = ' num2str(median(horsub)) ' pixels'])
xlabel('RF position - ypos (pixels)')

vertsub = kmap_vert_pix(find(bw(:)));
id = find(vertsub<prctile(vertsub,2) | vertsub>prctile(vertsub,98));
vertsub(id) = []; %Remove outliers

subplot(2,2,4)
hist(vertsub,40);
title(['median = ' num2str(median(vertsub)) ' pixels'])
xlabel('RF position - xpos (pixels)')



