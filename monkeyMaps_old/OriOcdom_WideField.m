function [oriang_raw ocdommap_raw dprimemask oriang ocdommap dphaseROI dphaseROISEL] = OriOcdom_WideField(dpthresh,varargin)

global bw f0m f0m_var funcmap ACQinfo symbolInfo Analyzer G_handles idExamp flipeyebit SelectivityThresh

[xmicperpix ymicperpix] = getImResolution(1);

set(G_handles.HPflag,'Value',0);
set(G_handles.LPflag,'Value',1);

set(G_handles.Lwidth,'string',num2str(9/xmicperpix*3.823));
hh = makeMapFilter;
set(G_handles.Lwidth,'string',num2str(.8/xmicperpix*3.823));
hh2 = makeMapFilter;


%% get orimap
for i = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp(Analyzer.loops.conds{1}.symbol{i},'ori');
        idsym = i;
        break
    end
end

symbolInfo.ID(1) = idsym;
set(G_handles.primSymbol,'value',idsym); 
if idsym == 1
    set(G_handles.secSymbol,'value',2);
else
    set(G_handles.secSymbol,'value',1);
end

setsymbolstruct

funcmap = GprocessAxis(f0m,hh);  %output is a vector image

oriang = angle(funcmap);
orimag = abs(funcmap);
oriang = (oriang+pi*(1-sign(oriang)))/2*180/pi;  %Put domain as [0 180].
orisel = abs(funcmap);

funcmap = GprocessAxis(f0m,hh2);  %output is a vector image

oriang_raw = angle(funcmap);
orimag_raw = abs(funcmap);
oriang_raw = (oriang_raw+pi*(1-sign(oriang_raw)))/2*180/pi;  %Put domain as [0 180].


%% get ocdom map

for i = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp(Analyzer.loops.conds{1}.symbol{i},'Leye_bit');
        idsym = i;
        break
    end
end

symbolInfo.ID(1) = idsym;
set(G_handles.primSymbol,'value',idsym); 
if idsym == 1
    set(G_handles.secSymbol,'value',2);
else
    set(G_handles.secSymbol,'value',1);
end
    
setsymbolstruct

ocdommap = GprocessBinary2(f0m,hh,hh);
ocdommap_raw = GprocessBinary2(f0m,hh2,hh);

%This is to make it so that positive is for the contralateral eye
if flipeyebit
    ocdommap = -ocdommap;
    ocdommap_raw = -ocdommap_raw;
end

edgeTrunc = 8;
dim = size(ocdommap);
trY = (edgeTrunc+1):(dim(1)-edgeTrunc);
trX = (edgeTrunc+1):(dim(2)-edgeTrunc);

%%
dim = size(oriang)
dim = size(oriang);
dsx = 1:10:dim(2);
dsy = 1:10:dim(1);
figure,

subplot(1,2,1)
imagesc(oriang), colormap hsv

subplot(1,2,2)
imagesc(ocdommap_raw)

hyp = 4;

for i = 1:length(dsy)
    for j = 1:length(dsx)
        
        hold on
        adj = hyp*cos(oriang(dsy(i),dsx(j))*pi/180 + pi/2);
        opp = hyp*sin(oriang(dsy(i),dsx(j))*pi/180 +pi/2);
        plot([dsx(j)-adj dsx(j)+adj],[dsy(i)-opp dsy(i)+opp] ,'k')
        
    end
end



%% Get mask for data selection

if isempty(varargin)
    [dprime2 dprimemask] = getMapMask(dpthresh,hh2,hh,edgeTrunc);
else
    dprimemask = varargin{1};
    set(G_handles.Lwidth,'string','1');
    h = makeMapFilter;
    dprime2 = ifft2(fft2(double(dprimemask)).*abs(fft2(h)));
    dprime2(:,1) = -.3;
end


%% plot

N_ODlines = 10;

mid = nanmean(ocdommap(:));
ocdommap_raw = ocdommap_raw - mid;
ocdommap = ocdommap - mid;

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

if ~isempty(idExamp)
    idExamp2(:,1) = idExamp(:,1)/ACQinfo.pixelsPerLine*xdom(end);
    idExamp2(:,2) = idExamp(:,2)/ACQinfo.pixelsPerLine*ydom(end);
end

dim = size(oriang);
oridom = 0:10:170;
ocdomdom = linspace(min(ocdommap(:)),max(ocdommap(:)),N_ODlines);

if isempty(bw)
    bw = ones(size(dprimemask));
end

%plot orimap w/ ocdom contour
figure
subplot(1,3,1)

leg = {'ow','sw','dw','pw'};
plotaxismap(dprime2(trY,trX),oriang_raw(trY,trX),0,xdom(trX),ydom(trY))
hold on
contour(xdom(trX),ydom(trY),ocdommap(trY,trX),ocdomdom,'k')
hold on
if ~isempty(idExamp)
    for i = 1:length(idExamp2(:,1))
        plot(idExamp2(i,1),idExamp2(i,2),leg{i},'markersize',10,'linewidth',2)
        hold on
    end
end

%[idy idx] = find(1-bw);  %cover up outside the ROI
%plot(idx*xmicperpix,idy*ymicperpix,'.w','MarkerSize',15)  
%hold off
title('Ori map w/ OcDom contour')

axis([0 1000 0 1000])  %This is to make all experiments scaled appropriately, relative to each other
axis off
plot([0 500],[995 995],'k')  %Plot scale bar
text(0,1050,'500um')  
 
%plot ocdommap w/ ori contour

subplot(1,3,2)
[x y] = meshgrid(-128:127,-128:127);

mi = prctile(ocdommap_raw(:),.5);
ma = prctile(ocdommap_raw(:),99.5);
mi = max([-1 mi]);
ma = min([1 ma]);
mima = max([-mi ma]);
plotODmap(dprime2(trY,trX),ocdommap_raw(trY,trX),0,xdom(trX),ydom(trY),[-mima mima])

%imagesc(xdom,ydom,sign(ocdommap))
hold on
%contour(xdom,ydom,(oriang-90)/1800,(oridom-90)/1800,'k')
contour(xdom(trX),ydom(trY),oriang(trY,trX),oridom,'k')
hold on
if ~isempty(idExamp)
    for i = 1:length(idExamp2(:,1))
        plot(idExamp2(i,1),idExamp2(i,2),leg{i},'markersize',10,'linewidth',2)
        hold on
    end
end
%[idy idx] = find(1-bw);  %cover up outside the ROI
%plot(idx*xmicperpix,idy*ymicperpix,'.w','MarkerSize',15)  
hold off
title('Ocdom map w/ Ori contour')
%Make colorbar
colormap jet
ocdombardom = round(linspace(-mima,mima,5)*1000)/1000;
clear domcell
for i = 1:length(ocdombardom)
    domcell{i} = ocdombardom(i);
end
iddom = linspace(1,max(oridom(:)),length(ocdombardom));
colorbar('YTick',iddom,'YTickLabel',domcell)

axis image
axis([0 1000 0 1000])  %This is to make all experiments scaled appropriately, relative to each other
axis off

%plot both ori/ocdom contours together
subplot(1,3,3)
contour(xdom(trX),ydom(trY),oriang(trY,trX),oridom,'k')
hold on
contour(xdom(trX),ydom(trY),ocdommap(trY,trX),ocdomdom,'r')
axis ij
hold on
[idy idx] = find(round(1-dprimemask(trY,trX)));  %cover up outside the ROI
plot(trX(idx)*xmicperpix,trY(idy)*ymicperpix,'.w','MarkerSize',15)  
axis image
title('Sfreq contour & Ori contour')

axis image
axis([0 1000 0 1000])  %This is to make all experiments scaled appropriately, relative to each other
axis off

%% plot anatomy
CH = GetTrialData([1 0 0 0],1);
imanat = (mean(CH{1},3)).^1.8;
figure,
%subplot(2,2,4)
imagesc(xdom(trX),ydom(trY),imanat(trY,trX),[prctile(imanat(:),2) prctile(imanat(:),99.5)]), colormap gray
axis image
axis([0 1000 0 1000])  %This is to make all experiments scaled appropriately, relative to each other


hold on
if ~isempty(idExamp)
    for i = 1:length(idExamp2(:,1))
        plot(idExamp2(i,1),idExamp2(i,2),leg{i},'markersize',10,'linewidth',2)
        hold on
    end
end

hold off

%% Compute bimodality

% %%
idROI = find(dprimemask>0);

dim = size(oriang);
%dorix = oriang(:,3:end) - oriang(:,1:end-2);
%doriy = oriang(3:end,:) - oriang(1:end-2,:);

dorix = oridiff(oriang(:,3:end)*pi/180,oriang(:,1:end-2)*pi/180);
doriy = oridiff(oriang(3:end,:)*pi/180,oriang(1:end-2,:)*pi/180);

dorix = [zeros(dim(1),1) dorix zeros(dim(1),1)];
doriy = [zeros(1,dim(2)); doriy; zeros(1,dim(2))];

magorigrad = sqrt(dorix.^2 + doriy.^2);
angorigrad = atan2(doriy,dorix);

dodx = ocdommap(:,3:end) - ocdommap(:,1:end-2);
dody = ocdommap(3:end,:) - ocdommap(1:end-2,:);

dodx = [zeros(dim(1),1) dodx zeros(dim(1),1)];
dody = [zeros(1,dim(2)); dody; zeros(1,dim(2))];

magodgrad = sqrt(dodx.^2 + dody.^2);
angodgrad = atan2(dody,dodx);

%Find regions with reliable gradients
set(G_handles.Lwidth,'string',num2str(20/xmicperpix*3.823));
h = makeMapFilter;
odGradSel = gradientReliability(ocdommap,h);
oriGradSel = gradientReliability(oriang,h);

idGrad = find(odGradSel> SelectivityThresh & oriGradSel > SelectivityThresh);
idHist = intersect(idGrad,idROI);  %Use idHist instead of idROI to limit to the regions with reliable gradients within the ROI


%% Compute intersection relation betweeen maps 

dphase = angle(exp(1i*angorigrad).*exp(-1i*angodgrad));  
dphaseROI = dphase(idROI)*180/pi;
dphaseROISEL = dphase(idHist)*180/pi;
plotGradIntersection(dphaseROI,dphaseROISEL);



function plotaxismap(mag,ang,anatflag,xdom,ydom,varargin)

global fh
%mag = log(mag)

%This is because of the funny stuff with bidirectional scanning
%mag = mag(3:end-2,3:end-2); ang = ang(3:end-2,3:end-2); 

if ~isempty(varargin)
    transparentMask = varargin{1};
    %transparentMask = transparentMask(3:end-2,3:end-2);
    transparentID = find(transparentMask == 1);
end


mag = phi(mag-prctile(mag(:),.05));
mag = mag/prctile(mag(:),98);
mag(find(mag>1)) = 1;

dim = size(ang);
set(gcf,'Color',[1 1 1]);

if anatflag
    
    CH = GetTrialData([1 0 0 0],1);
%     if get(G_handles.fastMotionFlag,'Value')
%         [Px_fast Py_fast] = getTrialMotion3(CH{1});
%         CH{1} = makeGeoTrx(CH{1},Px_fast,Py_fast);
%     end
    imanat = mean(CH{1}(:,:,2:end-1),3);
    
    %imanat = imanat(3:end-2,3:end-2);
    
    mi = prctile(imanat(:),0);    
    imanat = phi(imanat-mi);
    ma = prctile(imanat(:),100);
    imanat = imanat/ma;
    
    %%%
    mag = sqrt(imanat.*mag);
    %%%
    
    imfunc = ang;
    imfunc = imfunc/180;
    imfunc = round(imfunc*63+1);
    %imanat = round(imanat*63+1);

    hsvid = hsv;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)            
            imout(i,j,:) = mag(i,j)*hsvid(imfunc(i,j),:);
        end
    end
    
    imanat(:,:,2) = imanat;
    imanat(:,:,3) = imanat(:,:,1);
    
    %imout = 3*imout.^3+.3*(imanat).^.3;
    
    imout = imout + imanat;

    imout = imout/max(imout(:));
    
    %imout = imout(1:end-4,8:end,:);
    
    x = image(xdom,ydom,imout,'CDataMapping','direct','AlphaDataMapping','none');

else
%     imout = ang;
%     imout = imout/180;
%     imout = round(imout*63+1);
%     x = image(1:length(ang(1,:)),1:length(ang(:,1)),imout,'CDataMapping','direct','AlphaData',mag,'AlphaDataMapping','none');
    
    imfunc = ang;
    imfunc = imfunc/180;
    imfunc = round(imfunc*63+1);
    %imanat = round(imanat*63+1);

    hsvid = hsv;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)            
            imout(i,j,:) = mag(i,j)*hsvid(imfunc(i,j),:);
        end
    end
    
    if ~isempty(varargin)
        for i = 1:3
           imdum = imout(:,:,i);
           imdum(transparentID) = max(imout(:))*ones(size(transparentID)); 
           imout(:,:,i) = imdum;
        end
    end
    

    imout = imout/max(imout(:));
    
    %imout = imout(1:end-4,8:end,:);
    
    x = image(xdom,ydom,imout,'CDataMapping','direct');
    
end

axis image;

fh = gcf;

colormap hsv
%colorbar('YTick',[1 16:16:64],'YTickLabel',{'0','45','90','135','180'})

%Create the orientation legend%%%%%%%%%%%%%%%%%%
legdom = 0:30:180;
hsvdom = hsv;
id = round(linspace(1,64,length(legdom)));
hsvdom = hsvdom(id,:);
R = 20;
rid = linspace(1,ydom(end),length(legdom));
cid = xdom(end)+ 50;
xpts_o = [0 0];
ypts_o = [1-R 1+R];

for i = 1:length(legdom)
   
    xpts = xpts_o*cos(legdom(i)*pi/180) + ypts_o*sin(legdom(i)*pi/180);
    ypts = xpts_o*sin(legdom(i)*pi/180) - ypts_o*cos(legdom(i)*pi/180);
    ypts = ypts + rid(i);
    xpts = xpts + cid;
    hold on
    line(xpts,ypts,'Color',hsvdom(i,:),'Clipping','off','LineWidth',3);
    
end
hold off


