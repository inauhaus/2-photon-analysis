function [pColorLowSFAll pColorHiSFAll] = organizeF1vectors(rst45_135,rst0_90)
%%


pColorLowSFAll = [];
pColorHiSFAll = [];
dvecLowAll = [];
dvecHiAll = [];
oriidAll = [];
projectionLowAll = [];
projectionHiAll = [];

hiSFID = 3;
loSFID = 1;

for i = 1:length(rst45_135)-1

    Ncell = length(rst45_135{i}.opref{1});
    
  
    %%%Get response at preferred SF of M+S
 
    [dum sfids] = max(rst45_135{i}.sftc{1}(:,:)'); %ID peak location of SF response
    clear RespHiSFColor RespHiSFLum
    for j = 1:Ncell
        sfpkid = sfids(j);
        RespHiSFColor(j) =  rst45_135{i}.sftc{2}(j,sfpkid);
        RespHiSFLum(j) =  rst45_135{i}.sftc{1}(j,sfpkid);
    end
    RespLowSFLum = rst45_135{i}.sftc{1}(:,1)';
    RespLowSFColor = rst45_135{i}.sftc{2}(:,1)';
    
 
    %%Get preferred Ori
    
%     ColorSel_lowSF = rst45_135{i}.sftc{2}(:,loSFID) ./ (rst45_135{i}.sftc{1}(:,loSFID) +  rst45_135{i}.sftc{2}(:,loSFID));
%     ColorSel_hiSF = rst45_135{i}.sftc{2}(:,hiSFID) ./ (rst45_135{i}.sftc{1}(:,hiSFID) + rst45_135{i}.sftc{2}(:,hiSFID));
    

    ColorSel_hiSF = RespHiSFColor./(RespHiSFColor+RespHiSFLum);
    ColorSel_lowSF = RespLowSFColor./(RespLowSFColor+RespLowSFLum);
    

    pColorLowSFAll = [pColorLowSFAll ColorSel_lowSF(:)'];
    pColorHiSFAll = [pColorHiSFAll ColorSel_hiSF(:)'];
    
   
    [dum oriids] = max(rst45_135{i}.oritc{1}(:,:)'); %ID peak location of Ori response
    oriidAll = [oriids oriidAll];
    
    oriids2 = oriids-4;
    oriids2(find(oriids2<1)) = oriids2(find(oriids2<1))+8;
    
    clear MF1 SF1 MF12 SF12 dvec projectionV
    for j = 1:Ncell
         
        MF1(j) = squeeze(rst0_90{i}.sforimatF1{1}(oriids(j),loSFID,j))  /  abs(squeeze(rst0_90{i}.sforimat{1}(oriids(j),loSFID,j))) ;
        SF1(j) = squeeze(rst0_90{i}.sforimatF1{2}(oriids(j),loSFID,j))  /  abs(squeeze(rst0_90{i}.sforimat{2}(oriids(j),loSFID,j))) ;
        
        MF12(j) = squeeze(rst0_90{i}.sforimatF1{1}(oriids2(j),loSFID,j))  /  abs(squeeze(rst0_90{i}.sforimat{1}(oriids2(j),loSFID,j)));
        SF12(j) = squeeze(rst0_90{i}.sforimatF1{2}(oriids2(j),loSFID,j))  /  abs(squeeze(rst0_90{i}.sforimat{2}(oriids2(j),loSFID,j)));
        
        MF1(j) = squeeze(rst0_90{i}.sforimatF1{1}(oriids(j),loSFID,j))  ;
        SF1(j) = squeeze(rst0_90{i}.sforimatF1{2}(oriids(j),loSFID,j)) ;
        
        MF12(j) = squeeze(rst0_90{i}.sforimatF1{1}(oriids2(j),loSFID,j)) ;
        SF12(j) = squeeze(rst0_90{i}.sforimatF1{2}(oriids2(j),loSFID,j)) ;
        
        dvec1 = (MF1(j)*conj(SF1(j)));
        dvec2 = (MF12(j)*conj(SF12(j)));        
        dvec(j) = dvec1+dvec2;
        
        projectionV1 = real(MF1(j))*real(SF1(j)) + imag(MF1(j))*imag(SF1(j));
        projectionV2 = real(MF12(j))*real(SF12(j)) + imag(MF12(j))*imag(SF12(j));        
        projectionV(j) = projectionV1 + projectionV2;
    
    end
    dvecLowAll = [dvecLowAll dvec(:).']; %Use nonconjugate transpose
    projectionLowAll = [projectionLowAll projectionV(:)'];
   


    clear MF1 SF1 MF12 SF12 dvec projectionV
    for j = 1:Ncell
            
        MF1(j) = squeeze(rst0_90{i}.sforimatF1{1}(oriids(j),hiSFID,j))  /  abs(squeeze(rst0_90{i}.sforimat{1}(oriids(j),hiSFID,j))) ;
        SF1(j) = squeeze(rst0_90{i}.sforimatF1{2}(oriids(j),hiSFID,j))  /  abs(squeeze(rst0_90{i}.sforimat{2}(oriids(j),hiSFID,j)));
        
        MF12(j) = squeeze(rst0_90{i}.sforimatF1{1}(oriids2(j),hiSFID,j))  /  abs(squeeze(rst0_90{i}.sforimat{1}(oriids2(j),hiSFID,j)));
        SF12(j) = squeeze(rst0_90{i}.sforimatF1{2}(oriids2(j),hiSFID,j))  /  abs(squeeze(rst0_90{i}.sforimat{2}(oriids2(j),hiSFID,j))) ;
        
        MF1(j) = squeeze(rst0_90{i}.sforimatF1{1}(oriids(j),hiSFID,j)) ;
        SF1(j) = squeeze(rst0_90{i}.sforimatF1{2}(oriids(j),hiSFID,j)) ;
        
        MF12(j) = squeeze(rst0_90{i}.sforimatF1{1}(oriids2(j),hiSFID,j)) ;
        SF12(j) = squeeze(rst0_90{i}.sforimatF1{2}(oriids2(j),hiSFID,j)) ;
        
        dvec1 = (MF1(j)*conj(SF1(j)));
        dvec2 = (MF12(j)*conj(SF12(j)));        
        dvec(j) = dvec1+dvec2;
        
        projectionV1 = real(MF1(j))*real(SF1(j)) + imag(MF1(j))*imag(SF1(j));
        projectionV2 = real(MF12(j))*real(SF12(j)) + imag(MF12(j))*imag(SF12(j)); 

        projectionV(j) = projectionV1 + projectionV2;

    
    end
    dvecHiAll = [dvecHiAll dvec(:).']; %Use nonconjugate transpose
    projectionHiAll = [projectionHiAll projectionV(:)'];
    
    
%     MF1 = squeeze(rst0_90{i}.sforimatF1{1}(oriid,loSFID,:));
%     SF1 = squeeze(rst0_90{i}.sforimatF1{2}(oriid,loSFID,:));
%     pdiff = angle(MF1.*conj(SF1))*180/pi;
%     phasediffLowAll = [phasediffLowAll pdiff(:).']; %Use nonconjugate transpose
%     
%     
%     %figure,hist(pdiff)
%     
%     MF1 = squeeze(rst0_90{i}.sforimatF1{1}(oriid,hiSFID,:));
%     SF1 = squeeze(rst0_90{i}.sforimatF1{2}(oriid,hiSFID,:));
%     pdiff = angle(MF1.*conj(SF1))*180/pi;
%     phasediffHiAll = [phasediffHiAll pdiff(:).']; %Use nonconjugate transpose
    
    
end

phasediffHiAll = abs(angle(dvecHiAll))*180/pi;
phasediffLowAll = abs(angle(dvecLowAll))*180/pi;

projectionHiAll = real(dvecHiAll);
projectionLowAll = real(dvecLowAll);

figure,
subplot(2,1,1)
hist(phasediffHiAll)
subplot(2,1,2)
hist(phasediffLowAll)

figure,
subplot(2,1,1)
domMax = .01
hdom = linspace(-domMax,domMax,50);
hist(projectionHiAll,hdom),title('hiSF')
xlim([-domMax domMax])
subplot(2,1,2)
hist(projectionLowAll,hdom),title('LowSF')
xlim([-domMax domMax])


figure
for i = 1:4
    
    id = find(oriidAll == i | oriidAll == i+4);
    
    subplot(4,2,(i-1)*2+1)
    hist(phasediffLowAll(id))
    
    subplot(4,2,(i)*2)
    hist(phasediffHiAll(id))
    
end

figure
for i = 1:4
    
    id = find(oriidAll == i | oriidAll == i+4);
    
    subplot(4,2,(i-1)*2+1)
    hist(projectionLowAll(id),hdom)
    xlim([-domMax domMax])
    
    subplot(4,2,(i)*2)
    hist(projectionHiAll(id),hdom)
    xlim([-domMax domMax])
    
end
    
    
