function [sfprefdom orisig sizesig BW] = getOriSigfromModel(sf2sigModel)


sfprefdom = logspace(log10(1/2),log10(8),30);
sizesig = getSizefromSFpref(sf2sigModel,sfprefdom);
BW = 2./(sizesig*2*pi); %2sig bandwidth



xdom = linspace(-3,3,200);
[xdom ydom] = meshgrid(xdom,xdom);

clear orisig
for i = 1:length(sfprefdom)
    
    xsig = sizesig(i)
    %xsig = 1/(4*sfprefdom(i));
    ysig = xsig*1.7;
    GxyWeight = exp(-xdom.^2/(2*xsig^2)) .* exp(-ydom.^2/(2*ysig^2));
    
    carr = cos(xdom*sfprefdom(i)*2*pi);
    carr90 = sin(xdom*sfprefdom(i)*2*pi);
    
    Gab = GxyWeight.*carr;
    Gab90 = GxyWeight.*carr90;
    
    %subplot(6,6,i), imagesc(Gab)
    
    rfF = abs(fft2(Gab));
    rfF = fftshift(fftshift(rfF,1),2);
    
    rfF90 = abs(fft2(Gab90));
    rfF90 = fftshift(fftshift(rfF90,1),2);

    sfmodel = 'LinGauss';
    
    %[param ffit varacc sigmaraw] = oriFitfrom2D(rfF);
    [sfpref sfsig sigmaraw varacc_sf varacc_ori] = sforiFitfrom2D_2(rfF+rfF90,1,0,sfmodel);
    
    orisig(i) = sigmaraw;
    
end

%%
figure,plot(log2(sfprefdom),(orisig))

ylim([0 45])
