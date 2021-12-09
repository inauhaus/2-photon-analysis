function plotRandPosMaps

%Ian Nauhaus

global TC DM MK PIs ACQinfo idExamp Analyzer


global posprincax oriprincax


[xmicperpix ymicperpix] = getImResolution;

oripref = TC.OAng{1}{1};

%Now plot color image of tuning
orimagIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
oriprefIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
posmagIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
xposIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
yposIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sizeIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
BWIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);

xposIm_hat = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
yposIm_hat = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
oriprefIm_hat = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

xysize = [.33*xdom(end)/200 .34*ydom(end)/250];  %Size of the panels, normalized by the ROI size
%%
kElem = 1;
for c = 1:length(DM.colordom)
    
    %sizevals = min([TC.xsize{c}{kElem}; TC.ysize{c}{kElem}]);
    %sizevals = sqrt((TC.xsize{c}{kElem}.*TC.ysize{c}{kElem}));  %us
    sizevals = TC.profileSize{c}{kElem};
    
    H = [MK.CoM ones(length(MK.CoM(:,1)),1)];
    
    id = find(~isnan((TC.xpos{c}{kElem})) & ~isinf((TC.xpos{c}{kElem})));
    PIs.posplanex = inv(H(id,:)'*H(id,:))*H(id,:)'*(TC.xpos{c}{kElem}(id))'; %[dxdv dxdu]  u are columns, v are rows 
    [dum1,dum2,dum3,dum4,stats] = regress(TC.xpos{c}{kElem}(id)',H(id,:)); 
    r2_x = stats(1); pval_x = stats(3);
    xpos_hat = H*PIs.posplanex;    
    
    id = find(~isnan((TC.ypos{c}{kElem})) & ~isinf((TC.ypos{c}{kElem})));
    PIs.posplaney = inv(H(id,:)'*H(id,:))*H(id,:)'*(TC.ypos{c}{kElem}(id))';
    [dum1,dum2,dum3,dum4,stats] = regress(TC.ypos{c}{kElem}(id)',H(id,:)); 
    r2_y = stats(1); pval_y = stats(3);
    ypos_hat = H*PIs.posplaney;   
    
    xerr = TC.xpos{c}{kElem}'-xpos_hat;  %used later to compute the scatter
    yerr = TC.ypos{c}{kElem}'-ypos_hat;
    TC.xerr{c}{kElem} = xerr;
    TC.yerr{c}{kElem} = yerr;
    
    [xg yg] = meshgrid(1:length(xdom),1:length(ydom)); %keep it in pixel units
    ypos_hat2 = yg*PIs.posplaney(1) + xg*PIs.posplaney(2) + PIs.posplaney(3);
    xpos_hat2 = yg*PIs.posplanex(1) + xg*PIs.posplanex(2) + PIs.posplanex(3);
    
    %%%%Normalize the ranges based on the fits
    yposnorm = TC.ypos{c}{kElem} - min(ypos_hat2(:));
    xposnorm = TC.xpos{c}{kElem} - min(xpos_hat2(:));
    
    ypos_hat2 = ypos_hat2 - min(ypos_hat2(:));
    xpos_hat2 = xpos_hat2 - min(xpos_hat2(:));
    
    posrangeX = [0 max([ypos_hat2(:); xpos_hat2(:)])];
    posrangeY = posrangeX;
    
    xdomDum = linspace(posrangeX(1),posrangeX(2),5);
    ydomDum = linspace(posrangeY(1),posrangeY(2),5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %id = find(~isnan(TC.OAng{c}{kElem}));
    %H = H(id,:); opref = TC.OAng{c}{kElem}(id);  
    [oriplane oripref_hat] = Planefit(H(:,2),H(:,1),TC.OAng{c}{kElem});
    
    ori_ang = atan(oriplane(2)/oriplane(1))*180/pi;
    pos_ang = atan(PIs.posplanex(1)/PIs.posplanex(2))*180/pi;
    maporidiff = abs(oridiff(ori_ang*pi/180,pos_ang*pi/180)*180/pi);
    
    mapmag = TC.OMag{c}{kElem};
    mapmag = mapmag-min(mapmag);
    mapmag = mapmag/max(mapmag);  

    
   
    for p = 1:MK.Ncell

        idcell = find(MK.masklabel(:) == MK.celldom(MK.nID(p)));

        orimagIm(idcell) = mapmag(p);
        oriprefIm(idcell) = TC.OAng{c}{kElem}(p);
        
        oriprefIm_hat(idcell) = oripref_hat(p);
        
        
        if isnan(oriprefIm(idcell))
            oriprefIm(idcell) = 0;
        end

        posmagIm(idcell) = 1;
        
        xposIm(idcell) = xposnorm(p);        
        xposIm_hat(idcell) = xpos_hat(p);
        
        yposIm(idcell) = yposnorm(p);        
        yposIm_hat(idcell) = ypos_hat(p);
        
        sizeIm(idcell) = sizevals(p);
        
        BWIm(idcell) = TC.BWdiff{c}{kElem}(p);
        
    end

    orimagIm = log10(orimagIm+.01);
    orimagIm = orimagIm-min(orimagIm(:));
    orimagIm = orimagIm/max(orimagIm(:));
    
    %%%Orientation
    figure, subplot(3,2,1)
    IMtens = getImTens(oriprefIm,sign(orimagIm),[0 180],'hsv',1);
    image(xdom,ydom,IMtens)
    title(['color ' num2str(c)]),  axis image
    makeOriLegend(IMtens,xdom,ydom)
    plotEXcirc(idExamp,xdom,ydom)
    hold on
    hyp = 10; orig = -5;
    plot([orig hyp*cos(oriprincax*pi/180)+orig],[orig hyp*sin(oriprincax*pi/180)+orig],'k','Clipping','off')
    %set(gca,'Position',[.13 .58 xysize])
    
    subplot(3,2,2)
    IMtens = getImTens(oriprefIm_hat,sign(orimagIm),[0 180],'hsv',1);
    image(xdom,ydom,IMtens)    
    axis image
    plotEXcirc(idExamp,xdom,ydom)    
    %set(gca,'Position',[.53 .58 xysize])
    
    %%%Retinotopy X
    subplot(3,2,3)
    IMtens = getImTens((xposIm),sign(posmagIm),posrangeX,'jet',1);
    image(xdom,ydom,IMtens), colorbar, colormap jet
    for i = 1:length(xdomDum)
        domcell{i} = round(xdomDum(i)*100)/100;
    end    

    posvec = round((xdomDum)*100)/100;
    posvec(end) = floor((xdomDum(end))*100)/100;
    posvec(1) = ceil((xdomDum(1))*100)/100;
    posvec = (posvec - (posrangeX(1)))/((posrangeX(2))-(posrangeX(1)));    
    colorbar('YTick',posvec,'YTickLabel',domcell)
    title(['color ' num2str(c)]), axis image
    plotEXcirc(idExamp,xdom,ydom)
    hold on
    hyp = 10; orig = -5;
    plot([orig hyp*cos(posprincax*pi/180)+orig],[orig hyp*sin(posprincax*pi/180)+orig],'k','Clipping','off')
    %set(gca,'Position',[.13 .1 xysize])
    dxdu = PIs.posplanex(1)/xmicperpix*1000; %degrees per mm
    dxdv = PIs.posplanex(2)/xmicperpix*1000;
    slp = 1/sqrt(dxdu^2 + dxdv^2);
    title(['r^2 ' num2str(round(r2_x*100)/100) '; p ' num2str(pval_x) '; mm/deg ' num2str(round(slp*100)/100)  ])
    
    subplot(3,2,4)
    %IMtens = getImTens((xposIm_hat),sign(posmagIm),posrangeX,'jet',1);
    IMtens = getImTens(xpos_hat2,ones(size(xpos_hat2)),posrangeX,'jet',1); %plane
    image(xdom,ydom,IMtens)
    hold on
    contour(xdom,ydom,sign(MK.masklabel),[.5 .5],'k')
    
    %colorbar('YTick',posvec,'YTickLabel',domcell)
    axis image
    plotEXcirc(idExamp,xdom,ydom)    
    
    %set(gca,'Position',[.53 .1 xysize])
    
    
    %%%Retinotopy Y
    subplot(3,2,5)
    IMtens = getImTens(yposIm,sign(posmagIm),posrangeY,'jet',1);
    image(xdom,ydom,IMtens), colorbar, colormap jet
    for i = 1:length(ydomDum)
        domcell{i} = round(ydomDum(i)*100)/100;
    end

    posvec = round((ydomDum)*100)/100;
    posvec(end) = floor((ydomDum(end))*100)/100;
    posvec(1) = ceil((ydomDum(1))*100)/100;
    posvec = (posvec - (posrangeY(1)))/((posrangeY(2))-(posrangeY(1)));  
 
    colorbar('YTick',posvec,'YTickLabel',domcell)
    title(['color ' num2str(c)]), axis image
    plotEXcirc(idExamp,xdom,ydom)
    hold on
    hyp = 10; orig = -5;
    plot([orig hyp*cos(posprincax*pi/180)+orig],[orig hyp*sin(posprincax*pi/180)+orig],'k','Clipping','off')
    %set(gca,'Position',[.13 .1 xysize])
    dydu = PIs.posplaney(1)/ymicperpix*1000; %degrees per mm
    dydv = PIs.posplaney(2)/ymicperpix*1000;
    slp = 1/sqrt(dydu^2 + dydv^2);
    title(['r^2 ' num2str(round(r2_y*100)/100) '; p ' num2str(pval_y) '; mm/deg ' num2str(round(slp*100)/100)  ])
    
    subplot(3,2,6)
    %IMtens = getImTens((yposIm_hat),sign(posmagIm),posrangeY,'jet',1); plane of cells
    IMtens = getImTens(ypos_hat2,ones(size(ypos_hat2)),posrangeY,'jet',1); %plane
    image(xdom,ydom,IMtens)
    %colorbar('YTick',posvec,'YTickLabel',domcell)
    axis image
    plotEXcirc(idExamp,xdom,ydom)        
    hold on
    contour(xdom,ydom,sign(MK.masklabel),[.5 .5],'k')
    %set(gca,'Position',[.53 .1 xysize])
    
    %Compute magnification factor
    dxdu = PIs.posplanex(1)/xmicperpix*1000; %degrees per mm
    dxdv = PIs.posplanex(2)/xmicperpix*1000;
    dydu = PIs.posplaney(1)/ymicperpix*1000;
    dydv = PIs.posplaney(2)/ymicperpix*1000;
    PIs.MF = 1./abs(dxdu*dydv - dxdv*dydu);    
    subplot(3,2,6)
    title(['MagFac = '  num2str(PIs.MF) ' mm^2/deg^2'])
    
    %%%RF size
    %%%RF size
    figure,    
    %sizeImLog = log2(sizeIm); 
    sizeImLog = log2(1./sizeIm); 
    sizeImLog(find(isinf(-sizeImLog))) = 0;
    id = find(sizeImLog~=0);
    %posrange = [min(sizeImLog(id)) max(sizeImLog(id))];
    posrange = log2([.5 8]);
    IMtens = getImTens(sizeImLog,sign(abs(sizeImLog)),posrange,'jet',1);
    image(xdom,ydom,IMtens)
    title(['color ' num2str(c)]),  axis image
    plotEXcirc(idExamp,xdom,ydom)    
    posvec = linspace(0,1,5); 
    for i = 1:length(posvec)
        domcell{i} = posvec(i) * (posrange(2)-posrange(1)) + posrange(1);
        domcell{i} = round(2.^domcell{i}*100)/100;
    end
    colorbar('YTick',posvec,'YTickLabel',domcell), colormap jet
    title('RF size')
    
    %%%Off - On response
    figure,
    id = find(BWIm~=0);
    posrange = [min(BWIm(id)) max(BWIm(id))];
    posrange = [-1 1];
    IMtens = getImTens(BWIm,abs(sign(BWIm)),posrange,'jet',1);
    image(xdom,ydom,IMtens)
    title(['color ' num2str(c)]),  axis image
    plotEXcirc(idExamp,xdom,ydom)    
    posvec = linspace(0,1,5); 
    for i = 1:length(posvec)
        domcell{i} = posvec(i) * (posrange(2)-posrange(1)) + posrange(1);
        domcell{i} = round(domcell{i}*100)/100;
    end
    colorbar('YTick',posvec,'YTickLabel',domcell), colormap jet
    title('Off - On response')
    %hold on
    %set(gca,'Position',[.13 .58 xysize])
    
    %Plot Anatomy with example circles
    plotFilteredAnatomy(xdom,ydom)
    if ~isempty(idExamp)
        plotEXcirc(idExamp,xdom,ydom)
    else
        id = find(~isnan(TC.xpos{c}{kElem}));
        plotAllcirc(id,xdom,ydom)
    end
    
end

%% Point Image stuff

figure

% subplot(1,3,1), hist(sizevals,[.5:.1:2])
% title(['mean RFsize (2sigma) = ' num2str(mean(sizevals))])
% xlim([.5 2])

sizesig = sizevals/2;

%If we wanted to correct for bar width:
% barwidth = getparam('barWidth'); 
% barvar = (barwidth^2)/12;
% sizesig = sqrt(sizesig.^2 - barvar);

hdom = [-.2:.04:.2];
for i = 1:length(hdom)
   hdomstr{i} = num2str(hdom(i)); 
end
subplot(2,3,1), h = hist(xerr,hdom);
bar(hdom,h)
title(['mean((X - Xhat)^2) = ' num2str(nanmean(xerr.^2))])
xlim([hdom(1) hdom(end)])
xlabel('X - Xhat')
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.5 .5 .5])

hdom = [-.2:.04:.2];
for i = 1:length(hdom)
   hdomstr{i} = num2str(hdom(i)); 
end
subplot(2,3,2), h = hist(yerr,hdom);
bar(hdom,h)
title(['mean((Y - Yhat)^2) = ' num2str(nanmean(yerr.^2))])
xlim([hdom(1) hdom(end)])
xlabel('Y - Yhat')
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.5 .5 .5])

possig = sqrt(xerr.^2 + yerr.^2); %positional sigma
%PIs.posvarnorm = sqrt(posvar)/mean(sizevals/2);  %relative stdev
PIs.posvarnorm = (possig'./sizesig);  %we want relative variance
hdom = [0:.1:1.6];
for i = 1:length(hdom)
   hdomstr{i} = num2str(hdom(i)); 
end
subplot(2,3,3), h = histc(PIs.posvarnorm,hdom);
bar(h)
set(gca,'XTick',(1:length(hdom))-.5);
set(gca,'XTickLabel',hdomstr);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.5 .5 .5])
title(['mean(error/RFsig) = ' num2str(nanmean(PIs.posvarnorm))])
xlim([-.1 length(hdom)])
xlabel('dPos/RFsig from linear fit')

hdom = [0:.04:2];
for i = 1:length(hdom)
   hdomstr{i} = num2str(hdom(i)); 
end
subplot(2,3,4), h = hist(sizevals,hdom);
bar(hdom,h)
title(['mean RFsize (2sigma) = ' num2str(nanmean(sizevals))])
xlim([hdom(1) hdom(end)])
xlabel('RF size in degrees (2sigma)')
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.5 .5 .5])


subplot(2,3,5), 
h = hist(sizevals,hdom);
PIs.PI2sigma = 2*sqrt(nanmean(possig.^2) + nanmean(sizesig.^2))*sqrt(PIs.MF); %sqrt(sigscatter^2 + sigsize^2)*magfac
title(['Point-Image (2sigma) = ' num2str(PIs.PI2sigma) ' mm'])
axis off

%PI2sigma = sqrt(mean(sizesig.^2))*sqrt(PIs.MF);  %this is an approximation


PIs.oripref = oripref;



function dist = oridiff(angle1,angle2)

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;

function plotEXcirc(idExamp,xdom,ydom)

global MK

if ~isempty(idExamp)
    for q = 1:length(idExamp)
        hold on
        plot(xdom(round(MK.CoM(idExamp(q),2))),ydom(round(MK.CoM(idExamp(q),1))),'ok','MarkerSize',10)
    end
end

function plotAllcirc(idAll,xdom,ydom)

global MK

if ~isempty(idAll)
    for q = 1:length(idAll)
        hold on
        plot(xdom(round(MK.CoM(idAll(q),2))),ydom(round(MK.CoM(idAll(q),1))),'or','MarkerSize',20)
    end
end

function IMtens = getImTens(pref,mag,mima,maptype,BackG)

id = find(mag>1);
mag(id) = 1;

id = find(pref<mima(1));
pref(id) = mima(1);
id = find(pref>mima(2));
pref(id) = mima(2);


prefid = (pref-mima(1))/(mima(2)-mima(1));
prefid = round(prefid*63+1);  %normalize to be colormap index

dim = size(pref);
mapvals = eval(maptype);
IMtens = zeros(dim(1),dim(2),3);
for i = 1:dim(1)
    for j = 1:dim(2)
        if mag(i,j) == 0 | isnan(mag(i,j)) | isnan(prefid(i,j));
            IMtens(i,j,:) = [BackG BackG BackG];
        else
            IMtens(i,j,:) = mag(i,j)*mapvals(prefid(i,j),:);
        end
    end
end


function makeOriLegend(imTens,xdom,ydom)
%%
hold on
xlim manual
ylim manual
%Create the orientation legend%%%%%%%%%%%%%%%%%%
legdom = 0:30:180;
hsvdom = hsv;
id = round(linspace(1,64,length(legdom)));
hsvdom = hsvdom(id,:);
R = 20;
rid = linspace(ydom(1),ydom(end),length(legdom));
cid = xdom(end) + 30;
xpts_o = [0 0];
ypts_o = [1-R 1+R];

for i = 1:length(legdom)
   
    xpts = xpts_o*cos(legdom(i)*pi/180) + ypts_o*sin(legdom(i)*pi/180);
    ypts = -xpts_o*sin(legdom(i)*pi/180) + ypts_o*cos(legdom(i)*pi/180);
    ypts = ypts + rid(i);
    xpts = xpts + cid;
    hold on
    line(xpts,ypts,'Color',hsvdom(i,:),'Clipping','off','LineWidth',2);
    
end
hold off