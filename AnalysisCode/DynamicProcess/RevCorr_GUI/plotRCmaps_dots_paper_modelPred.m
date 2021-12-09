function plotRCmaps_dots_paper_modelPred

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

xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

[xmesh ymesh] = meshgrid(xdom,ydom);
%%
sfrange = [DM.sfdom(1) DM.sfdom(end)];
%sfdom = logspace(log10(sfrange(1)),log10(sfrange(end)),5);
sfdom = DM.sfdom;
%sfdom = [0.25 .5 1 2 4 8];
xysize = [.33*xdom(end)/200 .34*ydom(end)/250];  %Size of the panels, normalized by the ROI size

for c = 1:length(DM.colordom)
    c
    H = [MK.CoM ones(length(MK.CoM(:,1)),1)];
    H = [MK.CoM MK.CoM.^2 MK.CoM(:,1).*MK.CoM(:,2) ones(length(MK.CoM(:,1)),1)];
    

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
        
       
        F1F0Im(idcell) = TC.F1F0{c}(p);

        phaseIm(idcell) = TC.phase{c}(p);
        
        if isnan(oriprefIm(idcell))
            oriprefIm(idcell) = 0;
        end

        sfmagIm(idcell) = sfmag(p);
        sfprefIm(idcell) = TC.sfpref{c}(p);
        
    end

    orimagIm = log10(orimagIm+.01);
    orimagIm = orimagIm-min(orimagIm(:));
    orimagIm = orimagIm/max(orimagIm(:));
    
    %%%Orientation
    figure
    plotDotMap(TC.OAng{c},'hsv',MK.masklabel,MK.celldom(MK.nID),[0 180]);

    %IMtens = getImTens(oriprefIm,sign(orimagIm),[0 180],'hsv',1);
    %image(xdom,ydom,IMtens)
    title('orientation gratings'), 
    makeOriLegend(xdom,ydom)
    %plotEXcirc(idExamp,xdom,ydom)
    

    set(gca,'Position',[.4 .5 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
   
    
    %% SF map
    figure,     
    %IMtens = getImTens(log2(sfprefIm),sign(sfmagIm),log2(sfrange),'jet',1);
    %image(xdom,ydom,IMtens), colorbar, colormap jet
    plotDotMap(log2(TC.sfpref{c}),'jet',MK.masklabel,MK.celldom(MK.nID),log2(sfrange));
    for i = 1:length(sfdom)
        domcell{i} = round(sfdom(i)*100)/100;
    end    

    sfvec = round(log2(sfdom)*100)/100;
    sfvec(end) = floor(log2(sfdom(end))*100)/100;
    sfvec(1) = ceil(log2(sfdom(1))*100)/100;
    sfvec = (sfvec - log2(sfrange(1)))/(log2(sfrange(2))-log2(sfrange(1)));
    %sfvec = round(sfvec*63+1);    
    colorbar('YTick',sfvec,'YTickLabel',domcell), colormap jet
    title('sf')
    %plotEXcirc(idExamp,xdom,ydom)
    
    set(gca,'Position',[.4 .5 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    
    %% Size SI
    
    alpha = pi;
    
    figure
    subplot(2,2,1)
    %sizedom = sort(1./(2*sfdom));
    sizedom = [1/32 1/16 1/8 1/4 1/2 1];
    
    sizeEst = 1./(alpha*TC.sfpref{c}); %1/4 the period 1/(4sf) is about sigma RF size
    plotDotMap(log2(sizeEst),'jet',MK.masklabel,MK.celldom(MK.nID),log2([sizedom(1) sizedom(end)]));
    clear domcell
    for i = 1:length(sizedom)
        domcell{i} = round(sizedom(i)*100)/100;
    end    

    sfvec = round(log2(sizedom)*100)/100;
    sfvec(end) = floor(log2(sizedom(end))*100)/100;
    sfvec(1) = ceil(log2(sizedom(1))*100)/100;
    sfvec = (sfvec - log2(sizedom(1)))/(log2(sizedom(end))-log2(sizedom(1)));
    %sfvec = round(sfvec*63+1);    
    colorbar('YTick',sfvec,'YTickLabel',domcell), colormap jet
    title('Size estimate (sig) = 1/(4*preferredSF)'),
    %plotEXcirc(idExamp,xdom,ydom)
    
    set(gca,'Position',[.1 .6 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('size SI')
    %%%%% size pooled SI
    sigx = 0.23;
    
    subplot(2,2,2)
    %sizedom = sort(1./(2*sfdom));
    sizedom = [1/32 1/16 1/8 1/4 1/2 1];
    sizeEst = 1./(alpha*TC.sfpref{c}); %1/4 the period 1/(4sf) is about sigma RF size
    sizeEst = sqrt(sizeEst.^2 + sigx^2);
    plotDotMap(log2(sizeEst),'jet',MK.masklabel,MK.celldom(MK.nID),log2([sizedom(1) sizedom(end)]));
    clear domcell
    for i = 1:length(sizedom)
        domcell{i} = round(sizedom(i)*100)/100;
    end    

    sfvec = round(log2(sizedom)*100)/100;
    sfvec(end) = floor(log2(sizedom(end))*100)/100;
    sfvec(1) = ceil(log2(sizedom(1))*100)/100;
    sfvec = (sfvec - log2(sizedom(1)))/(log2(sizedom(end))-log2(sizedom(1)));
    %sfvec = round(sfvec*63+1);    
    colorbar('YTick',sfvec,'YTickLabel',domcell), colormap jet
    title('Size estimate (sig) = 1/(4*preferredSF)'),
    %plotEXcirc(idExamp,xdom,ydom)
        
    set(gca,'Position',[.6 .6 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('size pooled SI')
    %% sf BW map
    figure
    subplot(2,2,1)
    %IMtens = getImTens(log2(sfprefIm),sign(sfmagIm),log2(sfrange),'jet',1);
    %image(xdom,ydom,IMtens), colorbar, colormap jet
    plotDotMap(log2(TC.sfBWLin{c}/2),'jet',MK.masklabel,MK.celldom(MK.nID),log2(sfrange));
    clear domcell
    for i = 1:length(sfdom)
        domcell{i} = round(sfdom(i)*100)/100;
    end    

    sfvec = round(log2(sfdom)*100)/100;
    sfvec(end) = floor(log2(sfdom(end))*100)/100;
    sfvec(1) = ceil(log2(sfdom(1))*100)/100;
    sfvec = (sfvec - log2(sfrange(1)))/(log2(sfrange(2))-log2(sfrange(1)));
    %sfvec = round(sfvec*63+1);    
    colorbar('YTick',sfvec,'YTickLabel',domcell), colormap jet
    title('sf bandwidth; sigma'),
    %plotEXcirc(idExamp,xdom,ydom)
    
    set(gca,'Position',[.1 .6 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('BW')
    
    %%%%% sf BW SI map
    subplot(2,2,2)
    
    bwEst = TC.sfpref{c}/2; %1/4 the period 1/(4sf) is about sigma RF size
    %bwEst = sqrt(sizeEst.^2 + sigx^2);
    %IMtens = getImTens(log2(sfprefIm),sign(sfmagIm),log2(sfrange),'jet',1);
    %image(xdom,ydom,IMtens), colorbar, colormap jet
    plotDotMap(log2(bwEst),'jet',MK.masklabel,MK.celldom(MK.nID),log2(sfrange));
    clear domcell
    for i = 1:length(sfdom)
        domcell{i} = round(sfdom(i)*100)/100;
    end    

    sfvec = round(log2(sfdom)*100)/100;
    sfvec(end) = floor(log2(sfdom(end))*100)/100;
    sfvec(1) = ceil(log2(sfdom(1))*100)/100;
    sfvec = (sfvec - log2(sfrange(1)))/(log2(sfrange(2))-log2(sfrange(1)));
    %sfvec = round(sfvec*63+1);    
    colorbar('YTick',sfvec,'YTickLabel',domcell), colormap jet
    title('sf bandwidth; sigma'),
    %plotEXcirc(idExamp,xdom,ydom)
    
    set(gca,'Position',[.6 .6 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('BW SI')
    
    % sf BW pooled SI map
    
    sigf = 0.84;
    
    subplot(2,2,3)
    
    pbwEst = TC.sfpref{c}/2; %1/4 the period 1/(4sf) is about sigma RF size
    pbwEst = sqrt(pbwEst.^2 + sigf^2);
    %IMtens = getImTens(log2(sfprefIm),sign(sfmagIm),log2(sfrange),'jet',1);
    %image(xdom,ydom,IMtens), colorbar, colormap jet
    plotDotMap(log2(pbwEst),'jet',MK.masklabel,MK.celldom(MK.nID),log2(sfrange));
    clear domcell
    for i = 1:length(sfdom)
        domcell{i} = round(sfdom(i)*100)/100;
    end    

    sfvec = round(log2(sfdom)*100)/100;
    sfvec(end) = floor(log2(sfdom(end))*100)/100;
    sfvec(1) = ceil(log2(sfdom(1))*100)/100;
    sfvec = (sfvec - log2(sfrange(1)))/(log2(sfrange(2))-log2(sfrange(1)));
    %sfvec = round(sfvec*63+1);    
    colorbar('YTick',sfvec,'YTickLabel',domcell), colormap jet
    title('sf bandwidth; sigma'),
    %plotEXcirc(idExamp,xdom,ydom)
    
    set(gca,'Position',[.1 .1 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('BW pooled SI')

    %% ori BW map 
    figure
    subplot(2,2,1)
    orisigrange = [4 40];
    %IMtens = getImTens(log2(sfprefIm),sign(sfmagIm),log2(sfrange),'jet',1);
    %image(xdom,ydom,IMtens), colorbar, colormap jet

    %plotDotMap((TC.orisig{c}),'jet',MK.masklabel,MK.celldom(MK.nID),orisigrange);
    plotDotMap(log2(TC.orisig{c}),'jet',MK.masklabel,MK.celldom(MK.nID),log2(orisigrange));
    clear domcell
    orisigdom = [0:15:45];
    for i = 1:length(orisigdom)
        domcell{i} = round(orisigdom(i)*100)/100;
    end    

    orisigvec = round((orisigdom)*100)/100;
    orisigvec(end) = floor((orisigdom(end))*100)/100;
    orisigvec(1) = ceil((orisigdom(1))*100)/100;
    orisigvec = (orisigvec - (orisigrange(1)))/((orisigrange(2))-(orisigrange(1)));
    %sfvec = round(sfvec*63+1);    
    colorbar('YTick',orisigvec,'YTickLabel',domcell), colormap jet
    title('ori bandwidth (scale invariance); sigma'),
    %plotEXcirc(idExamp,xdom,ydom)
    
    set(gca,'Position',[.1 .6 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('oriBW')
    
    %%%%% ori BW map SI
    subplot(2,2,2)
    alpha = pi;
    A = 2;
    oriSI = atan(alpha/(2*pi*A))*180/pi;
    %IMtens = getImTens(log2(sfprefIm),sign(sfmagIm),log2(sfrange),'jet',1);
    %image(xdom,ydom,IMtens), colorbar, colormap jet

    %plotDotMap((TC.orisig{c}),'jet',MK.masklabel,MK.celldom(MK.nID),orisigrange);
    plotDotMap(log2(oriSI*ones(size((TC.orisig{c})))),'jet',MK.masklabel,MK.celldom(MK.nID),log2(orisigrange));
    clear domcell
    orisigdom = [0:15:45];
    for i = 1:length(orisigdom)
        domcell{i} = round(orisigdom(i)*100)/100;
    end    

    orisigvec = round((orisigdom)*100)/100;
    orisigvec(end) = floor((orisigdom(end))*100)/100;
    orisigvec(1) = ceil((orisigdom(1))*100)/100;
    orisigvec = (orisigvec - (orisigrange(1)))/((orisigrange(2))-(orisigrange(1)));
    %sfvec = round(sfvec*63+1);    
    colorbar('YTick',orisigvec,'YTickLabel',domcell), colormap jet
    title('ori bandwidth (scale invariance); sigma'),
    %plotEXcirc(idExamp,xdom,ydom)
    
    set(gca,'Position',[.6 .6 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('oriBW SI')
    
    %%%% ori BW model fit to complex cells
    subplot(2,2,3)
    %IMtens = getImTens(log2(sfprefIm),sign(sfmagIm),log2(sfrange),'jet',1);
    %image(xdom,ydom,IMtens), colorbar, colormap jet
    
    %plotDotMap((TC.orisig{c}),'jet',MK.masklabel,MK.celldom(MK.nID),orisigrange);

    pbwEst = TC.sfpref{c}/2; %1/4 the period 1/(4sf) is about sigma RF size
    pbwEst = sqrt(pbwEst.^2 + sigf^2);
    
    OriSig_model = atan(pbwEst./TC.sfpref{c}/A)*180/pi;
    plotDotMap(log2(OriSig_model),'jet',MK.masklabel,MK.celldom(MK.nID),log2(orisigrange));
    clear domcell
    orisigdom = [0:15:45];
    for i = 1:length(orisigdom)
        domcell{i} = round(orisigdom(i)*100)/100;
    end    

    orisigvec = round((orisigdom)*100)/100;
    orisigvec(end) = floor((orisigdom(end))*100)/100;
    orisigvec(1) = ceil((orisigdom(1))*100)/100;
    orisigvec = (orisigvec - (orisigrange(1)))/((orisigrange(2))-(orisigrange(1)));
    %sfvec = round(sfvec*63+1);    
    colorbar('YTick',orisigvec,'YTickLabel',domcell), colormap jet
    title('ori bandwidth (Cmplx Cell model); sigma'),
    %plotEXcirc(idExamp,xdom,ydom)
        
     set(gca,'Position',[.1 .1 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('oriBW pooled SI')
    
    
    %% F1F0

          
    %Now modulation ratio F1/F0 and phase
    figure
    subplot(2,2,1)
    plotDotMap(TC.F1F0{c},'jet',MK.masklabel,MK.celldom(MK.nID),[0 pi]);
    %plotEXcirc(idExamp,xdom,ydom)
    F1F0dom = {0 ,pi/4 ,pi/2 ,3*pi/4, pi};
    iddom = linspace(0,1,length(F1F0dom));
    colorbar('YTick',iddom,'YTickLabel',F1F0dom)
    colormap(jet)
    title('F1/F0')
    set(gca,'Position',[.1 .6 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('F1/F0')
   
    % F1F0 SI
 
    %Now modulation ratio F1/F0 and phase
    subplot(2,2,2)
    F1F0SI = ones(size(TC.F1F0{c}))*pi;
    F1F0SI(find(isnan((TC.F1F0{c})))) = NaN;
    plotDotMap(F1F0SI,'jet',MK.masklabel,MK.celldom(MK.nID),[0 pi]);
    %plotEXcirc(idExamp,xdom,ydom)
    F1F0dom = {0 ,pi/4 ,pi/2 ,3*pi/4, pi};
    iddom = linspace(0,1,length(F1F0dom));
    colorbar('YTick',iddom,'YTickLabel',F1F0dom)
    colormap(jet)
    title('F1/F0')
    set(gca,'Position',[.6 .6 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('F1/F0 SI')
    
    % F1F0 pooled SI
 
    D = 0.54;
    sigx = 0.24;
    
    %Now modulation ratio F1/F0 and phase
    subplot(2,2,3)
    F1F0model = pi*exp(-(sigx*2*pi*TC.sfpref{c}*D).^2/2);
    plotDotMap(F1F0model,'jet',MK.masklabel,MK.celldom(MK.nID),[0 pi]);
    %plotEXcirc(idExamp,xdom,ydom)
    F1F0dom = {0 ,pi/4 ,pi/2 ,3*pi/4, pi};
    iddom = linspace(0,1,length(F1F0dom));
    colorbar('YTick',iddom,'YTickLabel',F1F0dom)
    colormap(jet)
    title('F1/F0')
    set(gca,'Position',[.1 .1 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('F1/F0 pooled SI')
    
    %% log sf BW map
    figure
    %IMtens = getImTens(log2(sfprefIm),sign(sfmagIm),log2(sfrange),'jet',1);
    %image(xdom,ydom,IMtens), colorbar, colormap jet
    BWoct = (TC.sfBWLin{c}/2+TC.sfpref{c})./TC.sfpref{c};
    plotDotMap(log2(BWoct),'jet',MK.masklabel,MK.celldom(MK.nID),log2(sfrange));
    clear domcell
    for i = 1:length(sfdom)
        domcell{i} = round(sfdom(i)*100)/100;
    end    

    sfvec = round(log2(sfdom)*100)/100;
    sfvec(end) = floor(log2(sfdom(end))*100)/100;
    sfvec(1) = ceil(log2(sfdom(1))*100)/100;
    sfvec = (sfvec - log2(sfrange(1)))/(log2(sfrange(2))-log2(sfrange(1)));
    %sfvec = round(sfvec*63+1);    
    colorbar('YTick',sfvec,'YTickLabel',domcell), colormap jet
    title('logarithmic sf bandwidth; sigma'),
    %plotEXcirc(idExamp,xdom,ydom)
    
    set(gca,'Position',[.4 .5 .25 .25],'Xtick',[],'Ytick',[])
    colorbar off
    title('logBW')

end

% %% Anatomy plot
% CH = GetTrialData([1 0],5);
% im{1} = mean(CH{1}(:,:,2:end-1),3);
% anatIm = im{1};
% anatIm = LocalZ(anatIm,150);
% mi = prctile(anatIm(:),2);
% ma = prctile(anatIm(:),99.9);
% anatIm(find(anatIm>ma)) = ma;
% anatIm(find(anatIm<mi)) = mi;
% 
% LXcorr = maskS.im{1};
% mi = prctile(LXcorr(:),1);
% ma = prctile(LXcorr(:),99.5);
% LXcorr(find(LXcorr>ma)) = ma;
% LXcorr(find(LXcorr<mi)) = mi;
% LXcorr = LocalZ(LXcorr,50,1);
% mi = prctile(LXcorr(:),1);
% ma = prctile(LXcorr(:),99.5);
% LXcorr(find(LXcorr>ma)) = ma;
% LXcorr(find(LXcorr<mi)) = mi;
% 
% 
% figure,
% imagesc(xdom,ydom,anatIm), colormap gray, axis image
% 
% %plotEXcirc(idExamp,xdom,ydom)
% %Need circles to be red, not black
% if ~isempty(idExamp)
%     for q = 1:length(idExamp)
%         hold on
%         plot(xdom(round(MK.CoM(idExamp(q),2))),ydom(round(MK.CoM(idExamp(q),1))),'or','MarkerSize',15)
%     end
% end
% 
% imagesc(xdom,ydom,LXcorr), colormap gray, axis image
% if ~isempty(idExamp)
%     for q = 1:length(idExamp)
%         hold on
%         plot(xdom(round(MK.CoM(idExamp(q),2))),ydom(round(MK.CoM(idExamp(q),1))),'or','MarkerSize',15)
%     end
% end
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