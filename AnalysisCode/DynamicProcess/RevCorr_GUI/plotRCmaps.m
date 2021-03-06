function [maporidiff mapslopemag] = plotRCmaps

%Ian Nauhaus

global TC DM MK ACQinfo maskS idExamp

global sfprincax oriprincax

[xmicperpix ymicperpix] = getImResolution;

%Now plot color image of tuning
orimagIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
oriprefIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfmagIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfprefIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);

sfprefIm_hat = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
oriprefIm_hat = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);

F1F0Im = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
phaseIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;
%%
sfrange = [DM.sfdom(1) DM.sfdom(end)];
%sfdom = logspace(log10(sfrange(1)),log10(sfrange(end)),5);
sfdom = DM.sfdom;
xysize = [.33*xdom(end)/200 .34*ydom(end)/250];  %Size of the panels, normalized by the ROI size

for c = 1:length(DM.colordom)
    c
    H = [MK.CoM ones(length(MK.CoM(:,1)),1)];
    H = [MK.CoM MK.CoM.^2 MK.CoM(:,1).*MK.CoM(:,2) MK.CoM.^3 ones(length(MK.CoM(:,1)),1)];
    
    id = find(~isnan(log2(TC.sfpref{c})) & ~isinf(log2(TC.sfpref{c})));
    sfplane = inv(H(id,:)'*H(id,:))*H(id,:)'*log2(TC.sfpref{c}(id))';
    sfpref_hat = 2.^(H*sfplane);    

    %id = find(~isnan(TC.OAng{c}));
    %H = H(id,:); opref = TC.OAng{c}(id);  
    [oriplane oripref_hat] = Planefit(H(:,2),H(:,1),TC.OAng{c});
    
    ori_ang = atan(oriplane(2)/oriplane(1))*180/pi;
    sf_ang = atan(sfplane(1)/sfplane(2))*180/pi;
    maporidiff = abs(oridiff(ori_ang*pi/180,sf_ang*pi/180)*180/pi);
    
    mapslopemag = sqrt(oriplane(1)^2 + oriplane(2)^2) * sqrt(sfplane(1)^2 + sfplane(2)^2);
    
%     mapmag = mapmag-min(mapmag);
%     mapmag = mapmag/max(mapmag);

    mapmag = TC.OMag{c};
    id = find(~isnan(mapmag));
    mapmag(id) = 1;
    id = find(isnan(mapmag));
    mapmag(id) = 0;
    
    sfmag = TC.sfmag{c};
    id = find(~isnan(sfmag));
    sfmag(id) = 1;
    id = find(isnan(sfmag));
    sfmag(id) = 0;
    
    for p = 1:MK.Ncell

        idcell = find(MK.masklabel(:) == MK.celldom(MK.nID(p)));

        orimagIm(idcell) = mapmag(p);
        oriprefIm(idcell) = TC.OAng{c}(p);
        
        oriprefIm_hat(idcell) = oripref_hat(p);
        
        F1F0Im(idcell) = TC.F1F0{c}(p);
        phaseIm(idcell) = TC.phase{c}(p);
        
        if isnan(oriprefIm(idcell))
            oriprefIm(idcell) = 0;
        end

        sfmagIm(idcell) = sfmag(p);
        sfprefIm(idcell) = TC.sfpref{c}(p);
        
        sfprefIm_hat(idcell) = sfpref_hat(p);

    end

    orimagIm = log10(orimagIm+.01);
    orimagIm = orimagIm-min(orimagIm(:));
    orimagIm = orimagIm/max(orimagIm(:));
    
    %%%Orientation
    figure, subplot(2,2,1)
    IMtens = getImTens(oriprefIm,sign(orimagIm),[0 180],'hsv',1);
    image(xdom,ydom,IMtens)
    title(['color ' num2str(c)]),  axis image
    makeOriLegend(xdom,ydom)
    plotEXcirc(idExamp,xdom,ydom)
    %hold on
    %hyp = 10; orig = -5;
    %plot([orig hyp*cos(oriprincax*pi/180)+orig],[orig hyp*sin(oriprincax*pi/180)+orig],'k','Clipping','off')
    %set(gca,'Position',[.13 .58 xysize])
    
    subplot(2,2,2)
    IMtens = getImTens(oriprefIm_hat,sign(orimagIm),[0 180],'hsv',1);
    image(xdom,ydom,IMtens)    
    axis image
    plotEXcirc(idExamp,xdom,ydom)    
    %set(gca,'Position',[.53 .58 xysize])
    
    %%%Spatial frequency
    subplot(2,2,3)
    IMtens = getImTens(log2(sfprefIm),sign(sfmagIm),log2(sfrange),'jet',1);
    image(xdom,ydom,IMtens), colorbar, colormap jet
    for i = 1:length(sfdom)
        domcell{i} = round(sfdom(i)*100)/100;
    end    

    sfvec = round(log2(sfdom)*100)/100;
    sfvec(end) = floor(log2(sfdom(end))*100)/100;
    sfvec(1) = ceil(log2(sfdom(1))*100)/100;
    sfvec = (sfvec - log2(sfrange(1)))/(log2(sfrange(2))-log2(sfrange(1)));
    %sfvec = round(sfvec*63+1);    
    colorbar('YTick',sfvec,'YTickLabel',domcell)
    title(['color ' num2str(c)]), axis image
    plotEXcirc(idExamp,xdom,ydom)
    hold on
    hyp = 10; orig = -5;
    plot([orig hyp*cos(sfprincax*pi/180)+orig],[orig hyp*sin(sfprincax*pi/180)+orig],'k','Clipping','off')
    %set(gca,'Position',[.13 .1 xysize])
    
    subplot(2,2,4)
    IMtens = getImTens(log2(sfprefIm_hat),sign(sfmagIm),log2([.5 8]),'jet',1);
    image(xdom,ydom,IMtens)
    %colorbar('YTick',sfvec,'YTickLabel',domcell)
    axis image
    plotEXcirc(idExamp,xdom,ydom)    
    
    %set(gca,'Position',[.53 .1 xysize])
    
    %Now modulation ratio F1/F0 and phase
    IMtens = getImTens(F1F0Im,sign(F1F0Im),[0 2],'jet',1);
    figure, s1 = subplot(1,2,1); image(xdom,ydom,phi(IMtens)), title('F1/F0'), axis image 
    plotEXcirc(idExamp,xdom,ydom)
    F1F0dom = {0 ,.5 ,1 ,1.5, 2};
    iddom = linspace(0,1,length(F1F0dom));
    colorbar('YTick',iddom,'YTickLabel',F1F0dom)
    colormap(s1,jet)
    
    IMtens = getImTens(phaseIm,sign(F1F0Im),[-180 180],'hsv',1);
    s2 = subplot(1,2,2); image(xdom,ydom,phi(IMtens)), title('spatial phase'), axis image 
    plotEXcirc(idExamp,xdom,ydom)
    phasedom = {-180 ,-90 ,0 ,90, 180};
    iddom = linspace(0,1,length(phasedom));
    colorbar('YTick',iddom,'YTickLabel',phasedom)
    colormap(s2,hsv)
    
    figure,scatter((TC.sfpref{c}),TC.F1F0{c},'.k'), xlabel('sfpref'),ylabel('F1/F0')

    
    %% Fit Gaussian
    [sfpref id] = sort(TC.sfpref{c});
    F1F0 = TC.F1F0{c}(id);
    id = find(~isnan(F1F0.*sfpref));
    F1F0 = F1F0(id);
    sfpref = sfpref(id);
    
    sigdom = linspace(.1,3,100);
    for i = 1:length(sigdom)
        ffit = 2*exp(-(sfpref.^2)/(2*sigdom(i)^2));
        E(i) = norm(ffit-F1F0);
    end
    [dum id] = min(E);
    sig = sigdom(id);
    ffit = 2*exp(-(sfpref.^2)/(2*sig^2));
    figure,
    subplot(1,2,1)
    plot(sfpref,ffit), hold on, plot(sfpref,F1F0,'.k')
    xlabel('SF preference'), ylabel('F1/F0')
    title(['sig = ' num2str(sig) 'cyc/deg'])
    
%     [sfpref id] = sort(TC.sfpref{c});
%     F1F0 = TC.F1F0{c}(id);
%     id = find(~isnan(F1F0.*sfpref));
%     F1F0 = F1F0(id);
%     sfpref = sfpref(id);    
%     sfpref = [fliplr(-sfpref) sfpref];
%     F1F0 = [fliplr(F1F0) F1F0];    
%     [param ffit varacc sigma] = Gaussfit2(sfpref,F1F0);
%     figure,scatter(sfpref,ffit), hold on, plot(sfpref,F1F0,'.')


subplot(1,2,2), hist(F1F0,20), xlabel('F1/F0')
    
    %%
    
    %Now timetopeak
%     [dum ttp] = max(TC.tcourse{1},[],2);
%     IMtens = getImTens(ttp,sign(F1F0Im),[0 20],'jet',1);
%     figure, subplot(1,2,1), image(xdom,ydom,phi(IMtens)), title('time-to-peak'), axis image 
%     plotEXcirc(idExamp,xdom,ydom)
%     F1F0dom = {0 ,.5 ,1 ,1.5, 2};
%     iddom = linspace(1,64,length(F1F0dom));
%     colorbar('YTick',iddom,'YTickLabel',F1F0dom)

end

%% Anatomy plot
CH = GetTrialData([1 0],5);
im{1} = mean(CH{1}(:,:,2:end-1),3);
dum = im{1};
dum = LocalZ(dum,150);
mi = prctile(dum(:),2);
ma = prctile(dum(:),99.9);
dum(find(dum>ma)) = ma;
dum(find(dum<mi)) = mi;

figure,imagesc(xdom,ydom,dum), colormap gray, axis image
%plotEXcirc(idExamp,xdom,ydom)
%Need circles to be red, not black
if ~isempty(idExamp)
    for q = 1:length(idExamp)
        hold on
        plot(xdom(round(MK.CoM(idExamp(q),2))),ydom(round(MK.CoM(idExamp(q),1))),'or','MarkerSize',15)
    end
end
%%

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


function makeOriLegend(xdom,ydom)
%%
hold on

%Create the orientation legend%%%%%%%%%%%%%%%%%%
legdom = 0:30:180;
hsvdom = hsv;
id = round(linspace(1,64,length(legdom)));
hsvdom = hsvdom(id,:);
R = 14;
rid = linspace(ydom(1),ydom(end),length(legdom));
cid = xdom(end)+15;
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

xlim([xdom(1) xdom(end)])
ylim([ydom(1) ydom(end)])
hold off