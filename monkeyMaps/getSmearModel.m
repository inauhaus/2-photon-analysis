function [CmplxModel smearsig] = getSmearModel(sfpref,Rsize,sfdom,sf2sigModel,size_BW)

if strcmp(size_BW,'size')
    %First get the smear parameter
    SizePred = getSizefromSFpref(sf2sigModel,sfpref);  %get size predictions of simple cells
    smearsig = sqrt(nanmean(Rsize.^2 - SizePred.^2)); %Get smear from actual size
    
    SizePred = getSizefromSFpref(sf2sigModel,sfdom); %get size predictions of simple cells from domain
    CmplxModel = sqrt(SizePred.^2 + smearsig^2); %Smear them.
    
elseif strcmp(size_BW,'BW')
    
    %Rsize input should be 1 sigma
     
    %First get the smear parameter
    SizePred = getSizefromSFpref(sf2sigModel,sfpref);  %get size predictions of simple cells
    SizePred = 1./(SizePred*2*pi); %1 sig bandwidth
   
    smearsig = sqrt(nanmedian(Rsize.^2 - SizePred.^2)); %Get smear from actual size
    
    SizePred = getSizefromSFpref(sf2sigModel,sfdom); %get size predictions of simple cells from domain
    SizePred = 1./(SizePred*2*pi); %1sig bandwidth
    CmplxModel = sqrt(SizePred.^2 + smearsig^2); %Smear them.
    
    
    
        %First get the smear parameter
%     SizePred = getSizefromSFpref(sf2sigModel,sfpref);  %get size predictions of simple cells
%     SizePred = 1./(SizePred*2*pi); %1 sig bandwidth
%     
%     SizePred = (SizePred./sfpref);
%     Rsize = (Rsize./sfpref);
%    
%     smearsig = sqrt(nanmean(Rsize.^2 - SizePred.^2)); %Get smear from actual size
%     
%     SizePred = getSizefromSFpref(sf2sigModel,sfdom); %get size predictions of simple cells from domain
%     SizePred = 1./(SizePred*2*pi); %1sig bandwidth
%     SizePred = (SizePred./sfdom);
%     
%     CmplxModel = sqrt(SizePred.^2 + smearsig^2).*sfdom; %Smear them.

    
    
end

     


