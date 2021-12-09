function X = plot2pColorTuning

%%
global cellS

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;

KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;


Ncell = length(cellS.mukern);

clrdim = 1; %color dimension of ori v color matrix
oridim = 3-clrdim;

Nc = size(cellS.mukern{1},1);
dc = 360/Nc;
cdom = 0:dc:(360-dc);

clear X
figure
for pid = 2:Ncell
    p = pid-1;
    
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    %imagesc(cellS.mukern{i})
    
    kern = cellS.mukern{pid};
    oritc = mean(kern,clrdim);
    oritc = oritc/sum(oritc);
    if clrdim == 1
        oritcmat = ones(size(kern,1),1)*oritc;
    else
        oritcmat = oritc*ones(1,size(kern,2));
    end
    colortc = mean(oritcmat.*kern,oridim);
    colortc = colortc(1:Nc/2) + colortc((Nc/2+1):end);
    colortc = [colortc(:)' colortc(:)'];
    polarplot([cdom cdom(1)]*pi/180,[colortc colortc(1)],'.-'),
    
    dvec(p) = sum(colortc.*exp(1i*2*cdom*pi/180));
    dir(p) = angle(dvec(p))*180/pi/2;
    if dir(p)<0
        dir(p) = dir(p)+180;
    end
    title(num2str(round(dir(p))))
    
    Smag(p) = abs(cos(dir(p)*pi/180+pi*.5));
    Mmag(p) = abs(cos(dir(p)*pi/180));
    pS(p) = Smag(p)/(Mmag(p)+Smag(p)); %percent S
    
    Colormag(p) = abs(cos(dir(p)*pi/180-pi*.75));
    Lummag(p) = abs(cos(dir(p)*pi/180-pi*.25));
    pS(p) = Smag(p)/(Mmag(p)+Smag(p)); %percent S
    pC(p) = Colormag(p)/(Lummag(p)+Colormag(p)); %percent color
    
    X.colortc(p,:) = colortc;
   
    
end
X.colordom = cdom;
X.pS = pS;
X.pC = pC;
X.dvec = dvec;
X.dir = dir;

%%
global maskS ACQinfo

[xmicperpix ymicperpix] = getImResolution;

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
celldom = celldom(2:end);

clear CoM
for p = 1:length(celldom)
    [idcelly idcellx] = find(masklabel == celldom(p));
    CoM(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
end


%% Plot %S map

%pS = dir(1:end)/180;


SMmap = jet;
SMmap(:,1) = .4;
SMmap = SMmap(1:40,:);
SMmap = interp1((1:40)',SMmap,linspace(1,40,64)');
SMmapdum = SMmap;
SMmap(:,2) = SMmapdum(:,3);
SMmap(:,3) = SMmapdum(:,2);

pS(find(pS>1 | pS<0)) = NaN;
pSidx = ceil(pS*64);
cdom = SMmap;


figure
for p = 1:length(pS)
    
    if ~isnan(pS(p))
        
        vc = cdom(pSidx(p),:);
        
        %ax2 = subplot(2,2,3);
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*xmicperpix,'.','Color',vc,'MarkerSize',30)
        hold on
        xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
    end
    
end

xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij
title('%S')
vR = [0 1];
colorbar('Ticks',[0 .5 1],'TickLabels',[0 .5 1])
colormap(cdom)

