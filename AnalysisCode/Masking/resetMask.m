function resetMask

global maskS G_handles ACQinfo

[xmicperpix ymicperpix] = getImResolution;

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

bw = ZThresh(maskS.imZ,str2num(get(G_handles.maskThresh,'string')),str2num(get(G_handles.maskMorph,'string')));
bw = cellMorph(bw,str2num(get(G_handles.maskMorph,'string'))); 

%Apply ROI
bw{1} = bw{1}.*maskS.bw;
%bw{2} = bw{2}.*maskS.bw;

%minimum cell size (should be done after applying ROI)
eval(['mima = ' get(G_handles.minCellArea,'string') ' '])
maskS.bwCell = cellMinMaxSize(bw,mima);

figure(40), 
imagesc(xdom,ydom,maskS.im{1}), colormap gray
hold on
contour(xdom,ydom,maskS.bwCell{1},[.5 .5],'r')
axis image
hold off