function getOriSfdata(kernPop,kernSigPop,blank,blankSig)


global Analyzer maskS ACQinfo

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
[nID] = getNeuronMask;

matDim = size(kernPop);

Nsym = length(Analyzer.loops.conds{1}.symbol);

oridom = getdomain('ori');
sfdom = getdomain('sf_freq');

for i = 1:2
    
    switch matDim(i)
        case length(oridom)
            oriID = i;
        case length(sfdom)
            sfID = i;
    end        
end

%Get image dimensions
xVolt = ACQinfo.scanAmplitudeX/ACQinfo.zoomFactor;
yVolt = ACQinfo.scanAmplitudeY/ACQinfo.zoomFactor;
xmic = 94*xVolt-2;  %Kristina fit these lines
ymic = 135*yVolt+0.5;
xmicperpix = xmic/ACQinfo.pixelsPerLine;
ymicperpix = ymic/ACQinfo.linesPerFrame;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions (used to find maxima)

kernsmoother = getSmoother([.5 1 .5],[.5 1 .5],oridom,sfdom,oriID);

%Smoother before taking ori curve
orismoother = getSmoother([.1 1 .1],[.1 1 .1],oridom,sfdom,oriID);

%Smoother before taking sf curve
sfsmoother = getSmoother([.2 1 .2],[.2 1 .2],oridom,sfdom,oriID);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Ncell = matDim(end);

for p = 1:Ncell

    kdum = kernPop(:,:,p);
    kdumSig = kernSigPop(:,:,p);

    ma = max(kdum(:));
    [idmay idmax] = find(kdum == ma);

    %dprime = (ma-blank(i))/sqrt(kernSig(lcID,idma)^2/Nreps + blankSig(i)^2/Nblankreps);
    dprime = (ma-blank(p))/(kdumSig(idmay,idmax)+blankSig(p));

    if dprime > 3

        %get ori/sf curves
        kerndum = ifft2(fft2(kdum).*abs(fft2(kernsmoother)));
        [bestid1 bestid2] = find(kerndum == max(kerndum(:)));

        kernplotsf = ifft2(fft2(kdum).*abs(fft2(sfsmoother)));
        kernplotori = ifft2(fft2(kdum).*abs(fft2(orismoother)));
        if sfID == 1
            tcsf = squeeze(kernplotsf(:,bestid2));
            tcori = squeeze(kernplotori(bestid1,:));
        else
            tcsf = squeeze(kernplotsf(bestid1,:));
            tcori = squeeze(kernplotori(:,bestid2));
        end

        %process sf curves


        %     [dum sfpref(p)] = max(tcsf);
        %     sfpref(p) = sfdom(sfpref(p));
        %     sfmag(p) = (max(tcsf)-min(tcsf));

        [param ffit varacc ffitI domI pk BW] = DoGfit(tcsf',sfdom);

        if varacc > .8
            sfpref(p) = pk;
            sfmag(p) = (max(ffitI)-min(ffitI))/(max(ffitI)+min(ffitI));
        else
            sfpref(p) = NaN;
            sfmag(p) = 0;
        end

        %     figure(65)
        %     subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
        %     semilogx(sfdom,tcsf), hold on
        %     plot(domI,ffitI,'k')
        %     axis tight


        %process ori curves
        tcori = tcori(1:length(oridom)/2) + tcori(length(oridom)/2+1:end);
        tcori = phi(tcori-prctile(tcori,1));


        [OMag(p) OAng(p)] = orifind(tcori',oridom(1:end/2));


%         [param ffit varacc] = Gaussfit(oridom,tcori,1);
% 
%         param(1) = param(1)+oridom(1);

%         if varacc < .8
%             param = param*NaN;
%             OMag(p) = 0;
%         end
%         OAng(p) = param(1);


        %     figure(66)
        %     subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
        %     plot(oridom,tcori), hold on


        tcsfall(p,:) = tcsf;
        tcoriall(p,:) = tcori;
        
    else
        OAng(p) = NaN;
        OMag(p) = 0;
        
        sfpref(p) = NaN;
        sfmag(p) = 0;
    end
    
end


%%

%Now plot color image of tuning
mag = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
ang = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfmag2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfpref2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);


OMag = OMag-min(OMag);
OMag = OMag/max(OMag);
for p = 1:Ncell

    idcell = find(masklabel(:) == celldom(nID(p)));

    mag(idcell) = OMag(p);
    ang(idcell) = OAng(p);

    sfmag2(idcell) = sfmag(p);
    sfpref2(idcell) = sfpref(p);

end

%mag = log10(mag+.01);
mag = mag-min(mag(:));
mag = mag/max(mag(:));

figure,
imagesc(ang,'AlphaData',mag,[0 180]), colormap hsv, colorbar
axis image

figure
imagesc(log10(sfpref2),'Alphadata',sfmag2,log10([sfdom(1) sfdom(end)])), colorbar
for i = 1:length(sfdom)
    domcell{i} = round(sfdom(i)*100)/100;
end
axis image
sfvec = round(log10(sfdom)*100)/100;
sfvec(end) = floor(log10(sfdom(end))*100)/100;
sfvec(1) = ceil(log10(sfdom(1))*100)/100;
colorbar('YTick',sfvec,'YTickLabel',domcell)


%%
for p = 1:Ncell
    [idcelly idcellx] = find(masklabel == celldom(nID(p)));
    CoM(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
end

doriAll = []; dsfAll = []; doriEucAll = []; dsfEucAll = []; DistAll = [];
for i = 1:Ncell
    
    for j = i+1:Ncell
        
        dori = abs(oridiff(OAng(i)*pi/180,OAng(j)*pi/180)*180/pi);
        dsf = abs(log10(sfpref(i)/sfpref(j)));
        doriAll = [doriAll dori];
        dsfAll = [dsfAll dsf];
        
        v1 = tcoriall(i,:)/norm(tcoriall(i,:));  v2 = tcoriall(j,:)/norm(tcoriall(j,:));
        doriEuc = corrcoef(v1,v2); doriEuc = -doriEuc(1,2); doriEuc = (doriEuc+1)/2;
        v1 = tcsfall(i,:)/norm(tcsfall(i,:));  v2 = tcsfall(j,:)/norm(tcsfall(j,:));
        dsfEuc = corrcoef(v1,v2); dsfEuc = -dsfEuc(1,2); dsfEuc = (dsfEuc+1)/2;
        doriEucAll = [doriEucAll doriEuc];
        dsfEucAll = [dsfEucAll dsfEuc];
        
        dy = (CoM(i,1)-CoM(j,1))*ymicperpix;
        dx = (CoM(i,2)-CoM(j,2))*xmicperpix;
        Dist = sqrt(dy^2 + dx^2); %Dist between cells in microns
        DistAll = [DistAll Dist];
        
    end
end

Ddom = 0:50:250;
%Ddom = [0 50 90 187.5]
clear dori dsf doriEuc dsfEuc
for i = 1:length(Ddom)-1

    DLimitMin = Ddom(i);  %Limit pairs to be this distance (microns) apart
    DLimitMax = Ddom(i+1);  %Limit pairs to be this distance (microns) apart
    id = find(DistAll<DLimitMax & DistAll>DLimitMin & ~isnan(doriAll.*dsfAll) & ~isinf(doriAll.*dsfAll));

    dori{i} = doriAll(id);
    dsf{i} = dsfAll(id);

    doriEuc{i} = doriEucAll(id);
    dsfEuc{i} = dsfEucAll(id);

end
%%
plotD = 2;
[mat xdom ydom] = smoothscatter(dori{plotD},dsf{plotD},.8,.015);

Valstruct.dori = dori{plotD};
Valstruct.dsf = dsf{plotD};
Valstruct.Ddom = Ddom;

figure,
subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
xlabel('dori'), ylabel('dsf')
subplot(1,2,1),scatter(doriAll,dsfAll,'.k')
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
xlabel('dori'), ylabel('dsf')

%%
[mat xdom ydom] = smoothscatter(doriEuc{plotD},dsfEuc{plotD},.015,.015);

Valstruct.doriEuc = doriEuc{plotD};
Valstruct.dsfEuc = dsfEuc{plotD};

figure,
subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
xlabel('dori (Euc dist)'), ylabel('dsf (Euc dist)')
subplot(1,2,1),scatter(doriEucAll,dsfEucAll,'.k')
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
xlabel('dori (Euc dist)'), ylabel('dsf (Euc dist)')



function smoother = getSmoother(kori,ksf,oridom,sfdom,oriID)

if length(sfdom) == 1
    ksf = 1;
end
if length(oridom) == 1
    kori = 1;
end

ksf = [ksf zeros(1,length(sfdom)-length(ksf))];
kori = [kori zeros(1,length(oridom)-length(kori))];
smoother = kori'*ksf;
smoother = smoother/sum(smoother(:));
if oriID == 2
    smoother = smoother';
end

function [OMag OAng] = orifind(G,oridomain)

R = sum(G'.*exp(1i*oridomain*pi/90));
OAng = angle(R);                 %-pi to pi
OAng = OAng + pi*(1-sign(OAng+eps));  %0 to 2pi
OAng = OAng*90/pi;               %0 to 180
OMag = abs(R)/(sum(G));

function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;