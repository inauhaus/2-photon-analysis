function [oriang_raw orimag_raw sfpref_raw sfmag_raw dprimemask dphaseROI dphaseROISEL] = OriSf_WideField(dpthresh,varargin)

global bw f0m f0m_var funcmap ACQinfo symbolInfo Analyzer G_handles idExamp SelectivityThresh

[xmicperpix ymicperpix] = getImResolution(1);

set(G_handles.HPflag,'Value',0);
set(G_handles.LPflag,'Value',1);

set(G_handles.Lwidth,'string',num2str(9/xmicperpix*3.823));
hh = makeMapFilter;
set(G_handles.Lwidth,'string',num2str(.8/xmicperpix*3.823));
hh2 = makeMapFilter;

anatomyflag = 0;
bwCellPlot = ones(size(funcmap));

sfdom = getdomain('s_freq');
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

%% get sfmap

for i = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp(Analyzer.loops.conds{1}.symbol{i},'s_freq');
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

funcmap = GprocessLog(f0m,bwCellPlot,hh);   %output is complex
sfmag = real(funcmap);
sfpref = imag(funcmap);

funcmap = GprocessLog(f0m,bwCellPlot,hh2);   %output is complex
sfmag_raw = real(funcmap);
sfpref_raw = imag(funcmap);

% sfpref_raw = medfilt2(sfpref_raw,[3 3]);
% sfpref = medfilt2(sfpref,[3 3]);
% 
% id = find(sfpref_raw<0.5);
% sfpref_raw(id) = 0.5;
% id = find(sfpref_raw>8);
% sfpref_raw(id) = 8;
% 
% id = find(sfpref<0.5);
% sfpref(id) = 0.5;
% id = find(sfpref>8);
% sfpref(id) = 8;


edgeTrunc = 8;
dim = size(sfpref);
trY = (edgeTrunc+1):(dim(1)-edgeTrunc);
trX = (edgeTrunc+1):(dim(2)-edgeTrunc);

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

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

if ~isempty(idExamp)
    idExamp2(:,1) = idExamp(:,1)/ACQinfo.pixelsPerLine*xdom(end);
    idExamp2(:,2) = idExamp(:,2)/ACQinfo.pixelsPerLine*ydom(end);
end

dim = size(oriang);
logdom = logspace(log10(.5),log10(8),25);
oridom = 0:10:170;

if isempty(bw)
    bw = ones(size(dprimemask));
end

%plot orimap w/ sf contour
figure
subplot(1,3,1)

leg = {'ow','sw','dw','pw'};
plotaxismap(dprime2(trY,trX),oriang_raw(trY,trX),0,xdom(trX),ydom(trY))
hold on
contour(xdom(trX),ydom(trY),sfpref(trY,trX),logdom,'k')
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
title('Ori map w/ Sfreq contour')

axis([0 1000 0 1000])  %This is to make all experiments scaled appropriately, relative to each other
axis off
hold on, plot([0 500],[995 995],'k')  %Plot scale bar
text(0,1050,'500um') 
 
%plot sfmap w/ ori contour
sfpref_raw(1,1) = sfdom(1); sfpref_raw(1,end) = sfdom(end); 
subplot(1,3,2)
plotlogmap(dprime2(trY,trX),sfpref_raw(trY,trX),0,xdom(trX),ydom(trY))
hold on
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
title('Sfreq map w/ Ori contour')
%Make colorbar
colormap jet
sfbardom = [.5 1 2 4 8];
clear domcell
for i = 1:length(sfbardom)
    domcell{i} = sfbardom(i);
end
iddom = linspace(1,max(oridom(:)),length(sfbardom));
colorbar('YTick',iddom,'YTickLabel',domcell)

axis([0 1000 0 1000])  %This is to make all experiments scaled appropriately, relative to each other
axis off

%plot both ori/sf contours together
subplot(1,3,3)
contour(xdom(trX),ydom(trY),oriang(trY,trX),oridom,'k')
hold on
contour(xdom(trX),ydom(trY),sfpref(trY,trX),logdom,'r')
axis ij
hold on
[idy idx] = find(1-dprimemask(trY,trX));  %cover up outside the ROI
plot(trX(idx)*xmicperpix,trY(idy)*ymicperpix,'.w','MarkerSize',15)  
axis image
title('Sfreq contour & Ori contour')


axis([0 1000 0 1000])  %This is to make all experiments scaled appropriately, relative to each other
axis off

CH = GetTrialData([1 0 0 0],1);
imanat = (mean(CH{1},3)).^1.8;
figure,
%subplot(2,2,4)
imagesc(xdom(trX),ydom(trY),imanat(trY,trX),[prctile(imanat(:),2) prctile(imanat(:),99.5)]), colormap gray
axis off
axis image
hold on,  line([20 20],[0 100])

hold on
if ~isempty(idExamp)
    for i = 1:length(idExamp2(:,1))
        plot(idExamp2(i,1),idExamp2(i,2),leg{i},'markersize',10,'linewidth',2)
        hold on
    end
end




%%
% global sfdist oridist SFatPinwheel
% 
% aue = [Analyzer.M.anim Analyzer.M.unit Analyzer.M.expt];
% 
% switch aue
% 
%     case 'ab8000093'
%         id = 1;
%         pinwheelLocs = [68 93; 116 108; 148 61; 169 71; 204 13]; %x/y ab8 LH
%     case 'ab8001007'
%         id = 2;
%         pinwheelLocs = [63 109; 193 48; 192 209; 24 223; 34 156; 38 137]; %ab8 RH
%     case 'ab9000084'
%         id = 3;
%         pinwheelLocs = [86 69; 71 113; 63 222]; %ab9
%     case 'ab9000059'
%         id = 4;
%         pinwheelLocs = [70 83; 26 51]; %ab9
%     case 'ac0000047'
%         id = 5;
%         pinwheelLocs = [47 195; 73 206; 155 144; 85 94; 138 56; 165 56]; %ab8 RH
% end
% 
% % figure
% % imagesc(oriang), colormap hsv
% 
% SFatPinwheel{id} = [];
% for i = 1:length(pinwheelLocs(:,1))
%     %hold on, plot(pinwheelLocs(i,1),pinwheelLocs(i,2),'o')
%     
%     if dprimemask(pinwheelLocs(i,2),pinwheelLocs(i,1))
%         SFatPinwheel{id} = [SFatPinwheel{id} sfpref(pinwheelLocs(i,2),pinwheelLocs(i,1))];
%     end
%     
%     sfdist{id} = sfpref_raw(find(dprimemask));
%     oridist{id} = oriang_raw(find(dprimemask));
% end
% 


%% Compute bimodality
figure
id = find(dprimemask == 1);

hdom = linspace(log2(.5),log2(8),80);
pdf = hist(log2(sfpref_raw(id)),hdom);
%[pdf] = histc((sfpref_raw(id)),[0 .8 1.5 3 6 10]); pdf = pdf(1:end-1); hdom = log2([.5 1 2 4 8]);    

pdf = pdf/sum(pdf);

dhdom = hdom(2)-hdom(1);
stairs(hdom,pdf,'k')

set(gca,'Xtick',log2([.5 1 2 4 8]))
xlim(log2([.5 8]))

%set(gca,'Xtick',[-1:.2:1])

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[1 1 1],'EdgeColor',[0 0 0])
set(gca,'TickDir','out')

xlabel('octaves')

%[dip, p_value, xlow,xup]=HartigansDipSignifTest(pdf,1000);
%p_value

%%
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

dsfx = log2(sfpref(:,3:end)./sfpref(:,1:end-2));
dsfy = log2(sfpref(3:end,:)./sfpref(1:end-2,:));

dsfx = [zeros(dim(1),1) dsfx zeros(dim(1),1)];
dsfy = [zeros(1,dim(2)); dsfy; zeros(1,dim(2))];

magsfgrad = sqrt(dsfx.^2 + dsfy.^2);
angsfgrad = atan2(dsfy,dsfx);

%Find regions with reliable gradients
set(G_handles.Lwidth,'string',num2str(20/xmicperpix*3.823));
h = makeMapFilter;
sfGradSel = gradientReliability(log(sfpref),h);
oriGradSel = gradientReliability(oriang,h);

idGrad = find(sfGradSel> SelectivityThresh & oriGradSel > SelectivityThresh);
idHist = intersect(idGrad,idROI);

%% Plot intersection

dphase = angle(exp(1i*angorigrad).*exp(-1i*angsfgrad));  
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


%%
function plotlogmap(mag,pref,anatflag,xdom,ydom,varargin)

global fh symbolInfo


if ~isempty(varargin)
    transparentMask = varargin{1};
    transparentID = find(transparentMask == 1);
end

mag = phi(mag-prctile(mag(:),0));
mag = mag/prctile(mag(:),98);
mag(find(mag>1)) = 1;

pref = log2(pref);
dim = size(mag);
set(gcf,'Color',[1 1 1]);

id = find(isnan(pref));
mag(id) = 0;
pref(id) = min(pref(:));

if anatflag
    
    [imanat] = getExptMean([1 0 0 0],2);
    imanat = imanat{1};
    
    mi = prctile(imanat(:),1);    
    imanat = phi(imanat-mi);
    ma = prctile(imanat(:),99);
    imanat = imanat/ma;

    imfunc = pref;
    imfunc = imfunc-min(imfunc(:));
    imfunc = imfunc/max(imfunc(:));
    imfunc = round(imfunc*63+1);

    jetid = jet;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)
            imout(i,j,:) = mag(i,j)*jetid(imfunc(i,j),:);
        end
    end
    
    imanat(:,:,2) = imanat;
    imanat(:,:,3) = imanat(:,:,1);
    
    imout = imout+sqrt(imanat);

    imout = imout/max(imout(:));
    
    x = image(xdom,ydom,imout,'CDataMapping','direct','AlphaDataMapping','none');

else      

    imfunc = pref;
    imfunc = imfunc-min(imfunc(:));
    imfunc = imfunc/max(imfunc(:));
    imfunc = round(imfunc*63+1);

    jetid = jet;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)
            imout(i,j,:) = mag(i,j)*jetid(imfunc(i,j),:);
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
    
    image(xdom,ydom,imout,'CDataMapping','direct'); 
    
       
end
axis image


