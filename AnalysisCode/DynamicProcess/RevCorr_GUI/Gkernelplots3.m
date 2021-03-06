function Gkernelplots3

%Ian Nauhaus

global ACQinfo Analyzer cellS maskS G_RChandles Valstruct

%%%%

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
[nID] = getNeuronMask;

trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];

%%%%

%Get image dimensions
xVolt = ACQinfo.scanAmplitudeX/ACQinfo.zoomFactor;
yVolt = ACQinfo.scanAmplitudeY/ACQinfo.zoomFactor;
xmic = 94*xVolt-2;  %Kristina fit these lines
ymic = 135*yVolt+0.5;
xmicperpix = xmic/ACQinfo.pixelsPerLine;
ymicperpix = ymic/ACQinfo.linesPerFrame;

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string') ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with dtau spacing, and end at an estimate of kernDel(2)

delayWin = [100 500];  %assume the peak response is within this time window
[dum delayWinID(1)] = min(abs(delayWin(1)-taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-taudom));

logfileroot = get(G_RChandles.logfilePath,'string');
expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];

load([logfileroot Analyzer.M.anim '\' expt],'frate')

if ~exist('frate')
    frate = 60;
end

domains = getSeqInfo(trialdom);

%%%%%%%%%%%%%%%%%%%
oridom = domains{5}.oridom;
sfdom = domains{5}.sfdom;
phasedom = domains{5}.phasedom;
colordom = domains{5}.colordom;

for i = 1:length(colordom)
    for p = 1:length(cellS.muTime)
        kernC{i}{p} = squeeze(cellS.muTime{p}(:,:,:,i,:));
        kernSigC{i}{p} = squeeze(cellS.sigTime{p}(:,:,:,i,:));
    end    
end

Ncell = length(cellS.muTime);
NT = getnotrials;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions (used to find maxima)

kernsmooth = getSmoother([.5 1 .5],[.5 1 .5],[.2 1 .2],taudom,oridom,sfdom);

%Smoother before taking ori curve
orismoother = getSmoother([.1 1 .1],[1],[1],taudom,oridom,sfdom);

%Smoother before taking sf curve
sfsmoother = getSmoother([.5 1 .5],[1],[1],taudom,oridom,sfdom);

%Smoother before taking temporal curve
tsmoother = getSmoother([.1 1 .1],[1],[.3 1 .3],taudom,oridom,sfdom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorid = {'r','g','b'};
%plot sf curves
%%
figure
for c = 1:length(colordom)
    tcsfall{c} = zeros(Ncell,length(sfdom));
    
    for p = 1:Ncell
        subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
        kernplot = zeros(length(oridom),length(sfdom),length(phasedom),length(taudom));
        kernplot(:,:,:,:) = kernC{c}{p};
        kernplot = mean(kernplot,3); %average over phase
        kernplot = reshape(kernplot,length(oridom),length(sfdom),length(taudom)); %squeeze phase dimension
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        kernplot = ifftn(fftn(kernplot).*abs(fftn(sfsmoother)));        
        
        [ma idma] = max(kerndum,[],3); %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestsfid) + delayWinID(1) - 1;
        
        tcorisfraw = kernplot(:,:,tau);
        tcorisfraw = tcorisfraw-prctile(tcorisfraw(:),10);
        [u s v] = svd(tcorisfraw);                
        tcorisf = (v(:,1)*s(1,1)*u(:,1)')';  
        tcsf = (v(:,1)*u(1,1)*s(1,1))';
%         imagesc(tcorisf)
        
        %tcsf = squeeze(kernplot(bestoriid,:,tau));
        %tcsf = phi(tcsf-prctile(tcsf,10));
        tcsfall{c}(p,:) = tcsf;
        
        sfdomI = logspace(log10(sfdom(1)),log10(sfdom(end)),length(sfdom)*2);
        tcsfI = interp1(sfdom,tcsf,sfdomI,'spline');  %This is important to make the fits better
        
        %[dum sfpref{c}(p)] = max(tcsf);
        %sfpref{c}(p) = sfdom(sfpref{c}(p));
        sfmag{c}(p) = (max(tcsf)-min(tcsf));
        
        
        [param ffit varacc ffitI domI pk BW] = DoGfit(tcsfI,sfdomI);
        
        if varacc > .6
            sfpref{c}(p) = pk;        
            sfmag{c}(p) = (max(ffitI)-min(ffitI))/(max(ffitI)+min(ffitI));
        else
            sfpref{c}(p) = NaN;        
            sfmag{c}(p) = 0;
        end
        
%         tcsf = tcsf-min(tcsf);
%         tcsf = tcsf/sum(tcsf);
%         sfpref{c}(p) = sum(tcsf.*sfdom);
        
        semilogx(sfdom,tcsf,['.' colorid{c} '-']), hold on
        plot(domI,ffitI,'k')       
        axis tight
        
        %imagesc(tcorisf)
        
        
    end

end

%%
%plot ori curves
colorvec = zeros(Ncell,length(colordom));
figure
for c = 1:length(colordom)
    
    tcoriall{c} = zeros(Ncell,length(oridom));
    for p = 1:Ncell
        subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
        %mean over phase
        kernplot = zeros(length(oridom),length(sfdom),length(phasedom),length(taudom));
        kernplot(:,:,:,:) = kernC{c}{p};
        kernplot = mean(kernplot,3); %average over phase
        kernplot = reshape(kernplot,length(oridom),length(sfdom),length(taudom)); %squeeze phase dimension
        
        kernsigplot = zeros(length(oridom),length(sfdom),length(phasedom),length(taudom));
        kernsigplot(:,:,:,:) = kernSigC{c}{p};
        kernsigplot = mean(kernsigplot,3); %average over phase
        kernsigplot = reshape(kernsigplot,length(oridom),length(sfdom),length(taudom)); %squeeze phase dimension
        
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));  %smooth in all dimensions (used to find slice)
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        
        kernplot = ifftn(fftn(kernplot).*abs(fftn(orismoother))); %smooth over all dimensions but orientation
        kernsigplot = ifftn(fftn(kernsigplot).*abs(fftn(orismoother))); %smooth over all dimensions but orientation
        
        [ma idma] = max(kerndum,[],3);  %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestsfid) + delayWinID(1) - 1;
        
        tcori = squeeze(kernplot(:,bestsfid,tau));
        tcsigori = squeeze(kernsigplot(:,bestsfid,tau));
        
%         tcorisf = kernplot(:,:,tau);
%         tcorisf = tcorisf-prctile(tcorisf(:),10);
%         [u s v] = svd(tcorisf);
%         tcori = u(:,1)*s(1,1)*v(1,1);
        
%         tcsigorisf = kernsigplot(:,:,tau);
%         tcsigorisf = tcsigorisf-prctile(tcsigorisf(:),10);
%         [u s v] = svd(tcsigorisf);      
%         tcsigori = u(:,1)*s(1,1)*v(1,1);     
        
        %tcori = squeeze(mean(kernplot(:,1:2,tau),2));
        %tcsigori = squeeze(mean(kernsigplot(:,1:2,tau),2));
        
        orimuMat{c}(p,:) = tcori;
        orisigMat{c}(p,:) = tcsigori;        
        
        %tcori = phi(tcori-prctile(tcori,10));
        tcoriall{c}(p,:) = tcori;
        
        [OMag{c}(p) OAng{c}(p)] = orifind(tcori,oridom);
        
         [param ffit varacc] = Gaussfit(oridom,tcori',1);
         sigfit(p) = param(2);  %used later to compute dynamics
%  
%         param(1) = param(1)+oridom(1);
% %         
%         if varacc < .6
%             param = param*NaN;
%             OMag{c}(p) = 0;
%         end
%         OAng{c}(p) = param(1);
        
        
        plot(oridom,tcori,['.' colorid{c} '-']), hold on
 %       plot(oridom,ffit,'k')

        colorvec(p,c) = param(3);
    end
    
    if ~isempty(cellS.muTimeBlank{p})
        cellS.muBase(p) = mean(cellS.muTimeBlank{p}(tau:tau));
        cellS.sigBase(p) = mean(cellS.sigTimeBlank{p}(tau:tau));
    end
    
    OMag{c} = phi(OMag{c});
end

cellS.orimuMat = orimuMat{1};
cellS.orisigMat = orisigMat{1};


figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

    for c = 1:length(colordom)
        
        kernplot = zeros(length(oridom),length(sfdom),length(phasedom),length(taudom));
        kernplot(:,:,:,:) = kernC{c}{p};
        
        kerndum = kernplot;
        %kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));        
        tcourse = squeeze(mean(mean(mean(kerndum,1),2),3));
        [dum id] = max(tcourse);
        kerndum = kernplot(:,:,:,id);
        kerndum = reshape(kerndum,length(oridom),length(sfdom),length(phasedom)); %squeeze time dimension
        
        [ma idma] = max(kerndum,[],3);  %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        
        tcphase = squeeze(kerndum(bestoriid,bestsfid,:));
        
        plot(phasedom,tcphase,['.' colorid{c} '-']), hold on
        
    end

end

figure,semilogy(OAng{1},sfpref{1},'.'), xlabel('pref orientation'), ylabel('pref spatial frequency')
set(gca,'Ytick',round(sfdom*100)/100)
xlim([0 180])
ylim([0 sfdom(end)+.5])

%%

%     [dum ids] = sort(rand(1,length(sfpref{1})));
%     sfpref{1} = sfpref{1}(ids);

%Now plot color image of tuning
mag = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
ang = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfmag2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfpref2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);

for c = 1:length(colordom)
    OMag{c} = OMag{c}-min(OMag{c});
    OMag{c} = OMag{c}/max(OMag{c});
    for p = 1:Ncell

        idcell = find(masklabel(:) == celldom(nID(p)));

        mag(idcell) = OMag{c}(p);
        ang(idcell) = OAng{c}(p);
        
        if isnan(ang(idcell))
            ang(idcell) = 0;
        end

        sfmag2(idcell) = sfmag{c}(p);
        sfpref2(idcell) = sfpref{c}(p);

    end

    mag = log10(mag+.01);
    mag = mag-min(mag(:));
    mag = mag/max(mag(:));

    figure,
    imagesc(ang,'AlphaData',mag,[0 180]), colormap hsv, colorbar
    title(['color ' num2str(c)])
    axis image

    figure
    imagesc(1:length(sfpref2(:,1)),1:length(sfpref2(1,:)),log10(sfpref2),'Alphadata',sfmag2,log10([sfdom(1) sfdom(end)])), colorbar
    for i = 1:length(sfdom)
        domcell{i} = round(sfdom(i)*100)/100;
    end
    sfvec = round(log10(sfdom)*100)/100;
    sfvec(end) = floor(log10(sfdom(end))*100)/100;
    sfvec(1) = ceil(log10(sfdom(1))*100)/100;
    colorbar('YTick',sfvec,'YTickLabel',domcell)
    title(['color ' num2str(c)])
    axis image
    
end

%%
%Color analysis
if length(colordom) == 3
    colorSens = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
    for p = 1:Ncell
        idcell = find(masklabel(:) == celldom(nID(p)));

        for c = 1:length(colordom)

            kernplot = zeros(length(oridom),length(sfdom),length(phasedom),length(taudom));
            kernplot(:,:,:,:) = kernC{c}{p};
            kernplot = mean(kernplot,3); %average over phase
            kernplot = reshape(kernplot,length(oridom),length(sfdom),length(taudom)); %squeeze phase dimension
            kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
            kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window

            ma = var(kerndum(:)); %find maxima from smoothed version
            ma = (max(kerndum(:))-min(kerndum(:)))/(max(kerndum(:))+min(kerndum(:))); %find maxima from smoothed version

            colorvec(p,c) = ma;
        end
        colorvec(p,:) = phi(colorvec(p,:));
        colorvec(p,:) = colorvec(p,:)/max(colorvec(p,:));

        if length(colordom) == 3
            lumpref(p) = log((colorvec(p,2)/colorvec(p,3)));
            colorSens(idcell) = lumpref(p);  %(L-M)/(L+M)
        end
    end

    Colordir = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine,length(colordom));


    for c = 1:length(colordom)
        imdum = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
        for p = 1:Ncell
            idcell = find(masklabel(:) == celldom(nID(p)));
            imdum(idcell) = colorvec(p,c);

        end
        Colordir(:,:,c) = imdum;
    end

    figure,
    image(Colordir), axis image
    for i = 1:length(sfdom)
        domcell{i} = round(sfdom(i)*100)/100;
    end
    sfvec = round(log10(sfdom)*100)/100;
    sfvec(end) = floor(log10(sfdom(end))*100)/100;
    sfvec(1) = ceil(log10(sfdom(1))*100)/100;
    colorbar('YTick',sfvec,'YTickLabel',domcell)

    lo = prctile(lumpref,3); hi = prctile(lumpref,97);
    figure, imagesc(colorSens), axis image, colorbar
    title('color sensitivity')

    figure,scatter(sfpref{3},lumpref,'.k'), xlabel('Spat freq preference'),ylabel('Color Selectivity')


end
%%
for p = 1:Ncell
    [idcelly idcellx] = find(masklabel == celldom(nID(p)));
    CoM(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
end

doriAll = []; dsfAll = []; doriEucAll = []; dsfEucAll = []; DistAll = [];
for i = 1:Ncell
    
    for j = i+1:Ncell
        
        dori = abs(oridiff(OAng{1}(i)*pi/180,OAng{1}(j)*pi/180)*180/pi);
        dsf = abs(log10(sfpref{1}(i)/sfpref{1}(j)));
        doriAll = [doriAll dori];
        dsfAll = [dsfAll dsf];
        
        v1 = tcoriall{1}(i,:)/norm(tcoriall{1}(i,:));  v2 = tcoriall{1}(j,:)/norm(tcoriall{1}(j,:));
        doriEuc = corrcoef(v1,v2); doriEuc = -doriEuc(1,2); doriEuc = (doriEuc+1)/2;
        v1 = tcsfall{1}(i,:)/norm(tcsfall{1}(i,:));  v2 = tcsfall{1}(j,:)/norm(tcsfall{1}(j,:));
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
[mat xdom ydom] = smoothscatter(dori{1},dsf{1},.8,.015);

Valstruct.dori = dori;
Valstruct.dsf = dsf;
Valstruct.Ddom = Ddom;

figure,
subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
xlabel('dori'), ylabel('dsf')
subplot(1,2,1),scatter(doriAll,dsfAll,'.k')
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
xlabel('dori'), ylabel('dsf')

%%
[mat xdom ydom] = smoothscatter(doriEuc{1},dsfEuc{1},.015,.015);

Valstruct.doriEuc = doriEuc;
Valstruct.dsfEuc = dsfEuc;

figure,
subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
xlabel('dori (Euc dist)'), ylabel('dsf (Euc dist)')
subplot(1,2,1),scatter(doriEuc{1},dsfEuc{1},'.k')
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
xlabel('dori (Euc dist)'), ylabel('dsf (Euc dist)')

%%
ddom = oridom(2)-oridom(1);
%time courses
figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    for c = 1:length(colordom)

    kernplot = zeros(length(oridom),length(sfdom),length(phasedom),length(taudom));
    kernplot(:,:,:,:) = kernC{c}{p};
    
    kernplotSig = zeros(length(oridom),length(sfdom),length(phasedom),length(taudom));
    kernplotSig(:,:,:,:) = kernSigC{c}{p};
    
    kernplot = mean(kernplot,3); %average over phase    
    kernplot = reshape(kernplot,length(oridom),length(sfdom),length(taudom)); %squeeze phase dimension
    kernplotSig = mean(kernplotSig,3);
    kernplotSig = reshape(kernplotSig,length(oridom),length(sfdom),length(taudom)); %squeeze phase dimension
    
    kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
    kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
    kernplot = ifftn(fftn(kernplot).*abs(fftn(tsmoother)));
    kernplotSig = ifftn(fftn(kernplotSig).*abs(fftn(tsmoother)));
    
    [ma idma] = max(kerndum,[],3);
    [bestoriid bestsfid] = find(ma == max(ma(:)));
    [mi idmi] = min(kerndum,[],3);
    [worstoriid worstsfid] = find(mi == min(mi(:)));

    %x = 7;
    sigfit(p) = min([90 sigfit(p)]);
    Dflank = round(sigfit(p)/ddom);
    Dflank = max([Dflank 1]);
    
    flankoriid1 = bestoriid+Dflank;
    if flankoriid1 > length(oridom)
        flankoriid1 = flankoriid1 - length(oridom);
    end
    flankoriid2 = bestoriid-Dflank;
    if flankoriid2 < 1
        flankoriid2 = flankoriid2 + length(oridom);
    end

    tcourse_ma = squeeze(kernplot(bestoriid,bestsfid,:));
    tcourseSig_ma = squeeze(kernplotSig(bestoriid,bestsfid,:));
    tcourse_flank = (squeeze(kernplot(flankoriid1,bestsfid,:)) + squeeze(kernplot(flankoriid2,bestsfid,:)))/2;
    tcourseSig_flank = (squeeze(kernplotSig(flankoriid1,bestsfid,:)) + squeeze(kernplotSig(flankoriid2,bestsfid,:)))/2;
    tcourse_orth = squeeze(kernplot(worstoriid,worstsfid,:));
    tcourseSig_orth = squeeze(kernplotSig(worstoriid,worstsfid,:));
    tcourse_all = squeeze(mean(mean(kernplot(:,:,:),1),2));
    tcourseSig_all = squeeze(mean(mean(kernplotSig(:,:,:),1),2));
    
    

    tcoursemaMat{c}(p,:) = tcourse_ma;
    tcourseorthMat{c}(p,:) = tcourse_orth;
    
    tcourseflankMat{c}(p,:) = tcourse_flank;
    
    tcourseSigmaMat{c}(p,:) = tcourseSig_ma;
    tcourseSigorthMat{c}(p,:) = tcourseSig_orth;
    
    tcourseAllMat{c}(p,:) = tcourse_all;
    tcourseSigAllMat{c}(p,:) = tcourseSig_all;
    
%     tcourse_ma = mean(kernplot);
%     tcourse_orth = blankresp;
    
    plot(taudom,tcourse_ma,'.-'), hold on, plot(taudom,tcourse_orth,'.-r')
    xlim([taudom(1) taudom(end)])
    
    end

end

cellS.tcoursemaMat = tcoursemaMat{1};
cellS.tcourseflankMat = tcourseflankMat{1};
cellS.tcourseorthMat = tcourseorthMat{1};
cellS.tcourseSigmaMat = tcourseSigmaMat{1};
cellS.tcourseSigorthMat = tcourseSigorthMat{1};
cellS.tcourseAllMat = tcourseAllMat{1};
cellS.tcourseSigAllMat = tcourseSigAllMat{1};




function [OMag OAng] = orifind(G,oridomain)

R = sum(G'.*exp(1i*oridomain*pi/90));
OAng = angle(R);                 %-pi to pi
OAng = OAng + pi*(1-sign(OAng+eps));  %0 to 2pi
OAng = OAng*90/pi;               %0 to 180
OMag = abs(R)/(sum(G));


function smoother = getSmoother(ktau,kori,ksf,taudom,oridom,sfdom)

if length(sfdom) == 1
    ksf = 1;
end
if length(oridom) == 1
    kori = 1;
end

ktau = [ktau zeros(1,length(taudom)-length(ktau))];
ksf = [ksf zeros(1,length(sfdom)-length(ksf))];
kori = [kori zeros(1,length(oridom)-length(kori))];
kdum = kori'*ksf;
for i = 1:length(ktau)
    smoother(:,:,i) = kdum*ktau(i);
end
smoother = smoother/sum(smoother(:));


function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;
