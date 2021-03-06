function RCcolorplots

%Ian Nauhaus

global ACQinfo cellS MK DM kernC TC

delayWin = [100 500];  %assume the peak response is within this time window
[dum delayWinID(1)] = min(abs(delayWin(1)-DM.taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-DM.taudom));

kernsmooth = getSmoother(50,10,.1,DM.taudom,DM.oridom,DM.sfdom); %for establising time-to-peak

%Color analysis
if length(DM.colordom) == 3
    colorSens = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
    imMag = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
    for p = 1:MK.Ncell
        idcell = find(MK.masklabel(:) == MK.celldom(MK.nID(p)));

%         for c = 1:length(DM.colordom)
% 
%             kernplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
%             kernplot(:,:,:,:) = kernC{c}{p};
%             kernplot = mean(kernplot,3); %average over phase
%             kernplot = reshape(kernplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
%             kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
%             kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
% 
%             ma = var(kerndum(:)); %find maxima from smoothed version
%             %ma = (max(kerndum(:))-min(kerndum(:)))/(max(kerndum(:))+min(kerndum(:))); %find maxima from smoothed version
% 
%             colorvec(p,c) = ma;
%         end

        colorvec(p,:) = TC.tccolorall{1}(p,:);
        %colorvec(p,:) = phi(colorvec(p,:));
        colorvec(p,:) = colorvec(p,:)/max(colorvec(p,:));

        if length(DM.colordom) == 3
            %lumpref(p) = log((colorvec(p,2)/colorvec(p,3))); %log((L-M)/(L+M)) OR log(M/S)
            %lumpref(p) = log2((colorvec(p,1)/colorvec(p,2))); %log(S/(L-M)) OR log(L/M)
            lumpref(p) = (colorvec(p,1)-colorvec(p,2))/(colorvec(p,1)+colorvec(p,2));
            colorSens(idcell) = lumpref(p);
            if ~isnan(TC.phase{1}(p).*TC.phase{2}(p))
                imMag(idcell) = 1;
            end
        end
    end
    
    %Make plot for LM cone ratio or colorLum ratio
    lo = prctile(lumpref,3); hi = prctile(lumpref,97);
    figure, imagesc(colorSens,'AlphaData',imMag,[-1 1]), axis image, colorbar
    title('color sensitivity (L-M)')

    
    %Make plot for RGB
    Colordir = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine,length(DM.colordom));
    for c = 1:length(DM.colordom)
        imdum = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
        for p = 1:MK.Ncell
            idcell = find(MK.masklabel(:) == MK.celldom(MK.nID(p)));
            imdum(idcell) = colorvec(p,c);

        end
        Colordir(:,:,c) = imdum;
    end
    figure,
    image(Colordir), axis image, title('RGB')
    
    
    %Make plot for phase differential of cone inputs
    PhaseDiff = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
    magIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
    for p = 1:MK.Ncell
        idcell = find(MK.masklabel(:) == MK.celldom(MK.nID(p)));
        pd = abs(angle(exp(1i*pi/180*(TC.phase{1}(p)-TC.phase{2}(p))))*180/pi);
        if ~isnan(pd)
            PhaseDiff(idcell) = pd;
            magIm(idcell) = 1;
        else
            magIm(idcell) = 0;
        end
            
    end

    figure,
    imagesc(PhaseDiff,'AlphaData',magIm,[0 180]), axis image, title('LM phase difference')



end


function smoother = getSmoother(tausig,orisig,sfsig,taudom,oridom,sfdom)

if length(sfdom) == 1
    ksf = 1;
end
if length(oridom) == 1
    kori = 1;
end

ktau = exp(-(taudom-mean(taudom)).^2/(2*tausig^2)); %tausig is in ms

dom = log2(sfdom)-mean(log2(sfdom));
ksf = exp(-dom.^2/(2*sfsig^2)); %sfsig is in octaves

oridomdum = linspace(0,180,length(oridom)+1);
oridomdum = oridomdum(1:end-1);  %I use this one in case oridom wraps around
kori = exp(-(oridomdum-oridomdum(ceil(end/2))).^2/(2*orisig^2)); %orisig is in degrees

kdum = kori'*ksf;
smoother = zeros(length(kori),length(ksf),length(ktau));
for i = 1:length(ktau)
    smoother(:,:,i) = kdum*ktau(i);
end
smoother = smoother/sum(smoother(:));