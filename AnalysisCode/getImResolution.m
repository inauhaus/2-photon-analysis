function [xmicperpix ymicperpix] = getImResolution
global ACQinfo

%% Measured 1/28/16 with scanbox version 2

imWidth_um = 980; %microns (N.B.  this is the "unblanked" region)
imHeight_um = 716; %microns

imWidth_pix = 676; %pixels (N.B.  this is the "unblanked" region)
imHeight_pix = 512; %pixels

ymicperpix = imHeight_um/imHeight_pix; %Default
xmicperpix = imWidth_um/imWidth_pix; %Default

%% Scale by the magnification setting

magvec = [1 1.2 1.4 1.7 2 2.4 2.8 3.4 4 4.8 5.7 6.7 8];

%old scanbox, ????
%magvec = [1 2 3];


ymicperpix = ymicperpix/magvec(ACQinfo.SBInfo.config.magnification);
xmicperpix = xmicperpix/magvec(ACQinfo.SBInfo.config.magnification);




% 
% switch ACQinfo.SBInfo.config.magnification
%     
%     case 1
%         
%         ymicperpix = 723/ACQinfo.linesPerFrame;
%         xmicperpix = 1166/ACQinfo.pixelsPerLine;
%         
%     case 2       
%         
%         ymicperpix = 723/ACQinfo.linesPerFrame/;
%         xmicperpix = 511/ACQinfo.pixelsPerLine;
%         
%     case 3       
%         
%         ymicperpix = 364/ACQinfo.linesPerFrame;
%         xmicperpix = 511/ACQinfo.pixelsPerLine;
%         
%     case 4
%         
%         ymicperpix = 185/ACQinfo.linesPerFrame;
%         xmicperpix = 260/ACQinfo.pixelsPerLine;
%         
% end
