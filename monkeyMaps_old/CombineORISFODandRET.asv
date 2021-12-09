global bw f0m f0m_var funcmap ACQinfo symbolInfo Analyzer G_handles idExamp kernelsIm 

set(G_handles.datadir,'string','e:\2p_data\')
set(G_handles.analyzedir,'string','e:\2p_data\AnalyzerFiles\')

%set(G_handles.datadir,'string','C:\2p_data\')
%set(G_handles.analyzedir,'string','C:\2p_data\AnalyzerFiles\')

%%%%%%%First, left hemisphere%%%%%%%%%%%

idExamp = [];  %left hemi

%make sure to hit "Set Directory" in GUI

%% ac1; u000_115 and u000_138  

ex = 3

switch ex
    case 1
        anim = 'ac1';
        expt = 'u000_115';
    case 2
        anim = 'ac1';
        expt = 'u001_011';
    case 3
        anim = 'ac1';
        expt = 'u002_019';
end

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

Temp1 = GetTrialData([0 1 0 0],1); 
Temp1 = mean(Temp1{1}(:,:,3:end-2),3);

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 1.5;
[oriang orimag sfpref sfmag dpmask] = OriSf_WideField(dthresh);
[dum ocdommap_raw dprimemask dum2 ocdommap] = OriOcdom_WideField(dthresh);
[dum dum dum] = SfOcdom_WideField(dthresh);

%NOW the randpos experiment
%load expt

switch ex
    case 1
        expt = 'u000_138';
    case 2
        expt = 'u001_006';
    case 3
        expt = 'u002_025';
end

set(G_handles.loadexp,'string',expt)
Gsetdirectories

N = ACQinfo.linesPerFrame; NI = size(oriang,1);
[domx domy] = meshgrid(1:N,1:N); 
[domxI domyI] = meshgrid(linspace(1,N,NI),linspace(1,N,NI));

Temp2 = GetTrialData([0 1 0 0],1); 
Temp2 = mean(Temp2{1}(:,:,3:end-2),3);
Temp2I = interp2(domx,domy,Temp2,domxI,domyI); 

set(G_handles.searchRange,'String','30');
[mbest nbest] = getShiftVals(Temp2I.^2,Temp1.^2,[0 0]);  %squaring seems to really help sometimes

%Temp2Ishift = circshift(Temp2I,[-round(mbest) -round(nbest) 0]);

load(['C:\2ph_code\Beta\randposTensors\' anim '_' expt(2:end)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blockSize = 16; vaThresh = .90;
[xpos ypos omap BWmap] = plotRandPosMap_xy2(blockSize,vaThresh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(xpos,1); NI = size(oriang,1); D = round(NI/N);
[domxI domyI] = meshgrid(1:NI,1:NI);
[domx domy] = meshgrid(D/2:D:(NI-D/2),D/2:D:(NI-D/2)); 

xposI = interp2(domx,domy,xpos,domxI,domyI); 
xposIshift = circshift(xposI,[-round(mbest) -round(nbest) 0]); 
xposIshift(:,1:phi(round(-nbest))) = NaN;
xposIshift(:,end-phi(round(nbest))-1:end) = NaN;
xposIshift(1:phi(round(-mbest)),:) = NaN;
xposIshift(end-phi(round(mbest))-1:end,:) = NaN;

yposI = interp2(domx,domy,ypos,domxI,domyI);
yposIshift = circshift(yposI,[-round(mbest) -round(nbest) 0]);
yposIshift(:,1:phi(round(-nbest))) = NaN;
yposIshift(:,end-phi(round(nbest))-1:end) = NaN;
yposIshift(1:phi(round(-mbest)),:) = NaN;
yposIshift(end-phi(round(mbest))-1:end,:) = NaN;

oriI = interp2(domx,domy,omap,domxI,domyI); 
oriIshift = circshift(oriI,[-round(mbest) -round(nbest) 0]); 
oriIshift(:,1:phi(round(-nbest))) = NaN;
oriIshift(:,end-phi(round(nbest))-1:end) = NaN;
oriIshift(1:phi(round(-mbest)),:) = NaN;
oriIshift(end-phi(round(mbest))-1:end,:) = NaN;

BWI = interp2(domx,domy,BWmap,domxI,domyI); 
BWIshift = circshift(BWI,[-round(mbest) -round(nbest) 0]); 
BWIshift(:,1:phi(round(-nbest))) = NaN;
BWIshift(:,end-phi(round(nbest))-1:end) = NaN;
BWIshift(1:phi(round(-mbest)),:) = NaN;
BWIshift(end-phi(round(mbest))-1:end,:) = NaN;

subplot(1,2,1), imagesc(xposIshift), axis image
subplot(1,2,2), imagesc(yposIshift), axis image

%%%%%%%%%%%%%%%%%%%
hh = fspecial('gaussian',size(oriang),8); hh = hh/sum(hh(:));
oriangH = angle(ifft2(abs(fft2(hh)).*fft2(exp(1i*oriang*2*pi/180))))*180/pi/2;
id = find(oriangH<0); oriangH(id) = oriangH(id)+180;

sfprefH = ifft2(abs(fft2(hh)).*fft2(sfpref));
%%%%%%%%%%%%%%%%%%%%%%
%% Retinotopy vs. Ocdom

D = round(size(xposIshift,1)/size(ypos,1))/1;
[JacvecOD{ex} ocdomvec{ex}] = RetvsOcdom(xposIshift,yposIshift,ocdommap,D);

xposIshift = xposIshift-min(xposIshift(:));
yposIshift = yposIshift-min(yposIshift(:));
ma = max([xposIshift(:); yposIshift(:)]); cdom = linspace(0,ma,10);
figure,imagesc(abs(ocdommap-mean(ocdommap(:)))*100)
hold on,contour(xposIshift,cdom,'k'), hold on, contour(yposIshift,cdom,'w'), axis ij

%% Retinotopy vs. sfreq

D = round(size(xposIshift,1)/size(ypos,1))/1;
[JacvecSF{ex} SFvec{ex}] = RetvsSfpref(xposIshift,yposIshift,sfprefH,D);
 
xposIshift = xposIshift-min(xposIshift(:));
yposIshift = yposIshift-min(yposIshift(:));
ma = max([xposIshift(:); yposIshift(:)]); cdom = linspace(0,ma,20);
figure,imagesc(sfprefH)
hold on,contour(xposIshift,cdom,'k'), hold on, contour(yposIshift,cdom,'w'), axis ij
%% Retinotopy vs. ori

% D = round(size(xposIshift,1)/size(ypos,1))/2;
% dorivsdpos(xposIshift,yposIshift,oriangH,D)

[oridum{ex} prefAxisMF{ex}] = dorivsdpos(xpos,ypos,omap,1);


figure,imagesc(oriang), colormap hsv
hold on,contour(xposIshift,cdom,'k'), 
hold on, contour(yposIshift,cdom,'k'), axis ij
axis square

%%
plotRetOriStuff(oridum,prefAxisMF,JacvecSF,SFvec,JacvecOD,ocdomvec);

%% Global plane

globalPlane(xposIshift,yposIshift)

%% dBW vs. dPos

% D = round(size(xposIshift,1)/size(ypos,1))/2;
% dBWvsdpos(xposIshift,yposIshift,BWIshift,D)

dBWvsdpos(xpos,ypos,BWmap,1)

%% dBW vs. dORI

% D = round(size(xposIshift,1)/size(ypos,1))/2;
% dorivsdBW(oriangH,BWIshift,D)

id = find(isnan(xpos.*ypos));
omapdum = omap; omapdum(id) = NaN;
dorivsdBW(omap,BWmap,1)

%% BW vs. Ocdom

D = round(size(xposIshift,1)/size(ypos,1));

id = find(isnan(xposIshift.*yposIshift));
BWmapdum = BWIshift; BWmapdum(id) = NaN;
ocdommapdum = ocdommap;  ocdommapdum(id) = NaN;
BWvsOcDom(BWmapdum,ocdommapdum,D)

figure,imagesc(abs(BWIshift))
hold on, contour(ocdommap/10,'k')
%% BW vs. SF

D = round(size(xposIshift,1)/size(ypos,1));

id = find(isnan(xposIshift.*yposIshift));
BWmapdum = BWIshift; BWmapdum(id) = NaN;
BWvsSF(BWmapdum,sfprefH,D)

figure,imagesc(abs(BWIshift))
hold on, contour(log2(sfprefH)/10,'k')