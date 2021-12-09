function [percRod Isom LumDom pupDiDom] = getIsomerizationRate_LCD(pupilDi)

%This was puzzling me:
%Is the readout of the spectroradiometer 1) the amount of power in each
%wavelength "bin", or 2) a "sample point".  The former puts it in units of
%Watts/ster/m2.  The latter is Watts/ster/m2/nm, in which case I need to
%multiply by the bin width (nm) to get Watts in each bin. I figured this
%out two ways.  First, I compared to "integrated Watt" that shows up on the LCD
%screen.  I only got the correct answer when scaling by bin width.  Second,
%I computed luminance (see below) and compared to the machines readout. I
%only got the same answer if I scale by the bin width.  So, unit readout is
%Watts/ster/m2/nm.
% 
% [UVfile Greenfile] = getRawDataLocation(anim,unit)
% 
% % Load power measurements (%Units are Watts/steradian/m2 from spectroradiometer)
% %pupilDi = 2;
% 
% %root = '/Users/in2293/Documents/MATLAB/Stimulator_master/Calibration/UVprojector/06_07_18/';
% UVfile = [UVfile 'spectrum_UV128'];   
% load(UVfile,'I','dom')
% UV = I;
% %root = '/Users/in2293/Documents/MATLAB/Stimulator_master/Calibration/RGBprojector/06_07_18/';
% RGBfile = [Greenfile 'spectrum_Green128'];
% load(RGBfile,'I','dom')
% green = I;  %W/steradian/micron2/dlambda


load('C:\Stimulator_master\Calibration\Oldmeasurements\Salk monitors\CRT 6-9-10 PR701\spectrum_red','I','dom')
IR = I;
load('C:\Stimulator_master\Calibration\Oldmeasurements\Salk monitors\CRT 6-9-10 PR701\spectrum_green','I','dom')
IG = I;
load('C:\Stimulator_master\Calibration\Oldmeasurements\Salk monitors\CRT 6-9-10 PR701\spectrum_blue','I','dom')
IB = I;

Rad_128 = IR + IG + IB;
%%


%Compute luminance%%%%%%%%%%%%%%%%%
Vlam = 1.019*exp(-285.4*(dom(:)/1000-0.559).^2);
Km = 683;
specRad = Rad_128(:);
Luminance = Km*(Vlam'*specRad)*(dom(2)-dom(1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LumDom = 1:1:1000;

specRadMat = (specRad(:)/Luminance)*LumDom;

figure
plot(dom,specRadMat)
ylabel('W/m^2/steradian')


%% Make figure
domdum = 300:750;
[rod,Mcone,Scone] = photoreceptor_templates(500,508,360,domdum);
specdum = specRadMat(:,1);
specdum = specdum/max(specdum);
specdum   = interp1(dom',specdum,domdum'); 
specdum(find(isnan(specdum))) = 0;
figure,plot(domdum,specdum,'k')
hold on
plot(domdum, rod)
xlim([300 750])

%%



% Get absorption functions (Govardovskii)
lambda = 300:700;
id = find(lambda>=dom(1) & lambda<=dom(end));
lambda = lambda(id);

[rod,Mcone,Scone] = photoreceptor_templates(500,508,360,lambda);

% Scale absorption functions by end-on collection area at the peak (um^2) (see Naarendorp and Wang et al)
Sens = [Mcone'*1 Scone'*1 rod'*0.85]; 

% Interpolate power measurements.  UV has a sharp peak that we don't want to miss.
specRadMat_I   = interp1(dom',specRadMat,lambda','spline'); 

clear dom

% Account for lens attenuation%%%%%%%%%%%%%%%%%%%%
%Lei and Yao 2006
a = .01939; b = 2.426*10^-3; c = 2.7277*10^4; d = -2.270; 
lamo = 353.6; alpha = 60;
muo = .00821; 
muh = .00000;
D = 2.07; %Lens thickness
at = a./(1+b*(lambda-lamo).^2) + c*lambda.^d + muo + alpha*muh; %attenuation 
T = 10.^-(at*D); %Transmittance from Beers law
T = T*100;
%T = T/max(T);
%figure,plot(lambda,T), ylim([0 100]), ylabel('transmittance')

%Jacobs and Williams 2007
% T = [58.1 63.8 68 71.9 75.2 78.1 80.3 82.0 83.5 85.1 86.4 87.6 88.7 89.6 90.2 91.3 92.1 92.8 93.4 93.9 94.5 95 95.5 96.1 96.6 97.1 97.6 97.9 98.1 98.5 98.8 99.2 99.4 99.6 99.8 100];
% wl = 350:10:700;
% T  = interp1(wl,T,lambda,'linear');
% figure,plot(lambda,T), ylim([0 100]), %transmission

specRadMat_I = specRadMat_I.*(T'*ones(1,size(specRadMat_I,2)))/100;

% Convert to flux at the retina
%pupilDi = 2;
eyeDi = 3.4; % mm (Remtulla)

Apupil = pi*(pupilDi/2)^2;

converter = Apupil/(eyeDi)^2  %unitless

%Wperum2 = [green(:) UV(:)] * (pupilDi/eyeDi)^2 / (10^12);  %flux at the retina (W/um^2/dlambda)

Wperm2 = specRadMat_I * converter; %Convert display radiance to (irradiance) at the retina (W/m^2/dlambda)

Wperum2 = Wperm2 / (10^12);  %flux at the retina (W/um^2/dlambda)


%Convert to quanta
lambdaMat = lambda(:)*ones(1,size(specRadMat_I,2))*10^-9;  %units of meters
Photperum2perlam = (Wperum2.*lambdaMat)/(3*10^8)/(6.63*10^-34); %convert to quanta (photons/sec/um^2/dlambda)

%figure,plot(lambda,Photperum2perlam)
%ylabel('photons/sec/um2/dlambda')

dlam = lambda(2)-lambda(1);
IsomerizationRate_Mopsin = Sens(:,1)'*Photperum2perlam*dlam %R*/sec/receptor for M cones, S cones, and rods 
IsomerizationRate_Sopsin = Sens(:,2)'*Photperum2perlam*dlam %R*/sec/receptor for M cones, S cones, and rods 
IsomerizationRate_rods = Sens(:,3)'*Photperum2perlam*dlam %R*/sec/receptor for M cones, S cones, and rods 

figure,

subplot(3,1,2)
semilogy(LumDom,[IsomerizationRate_Mopsin' IsomerizationRate_Sopsin' IsomerizationRate_rods']), xlabel('luminance'), ylabel('isomerization rate')
legend('M-opsin','S-opsin','rods')

hold on
plot([LumDom(1) LumDom(end)],[10*10^3 10*10^3 ],'--k')
hold on
plot([LumDom(1) LumDom(end)],[400*10^3 400*10^3 ],'--k')

ylim([10^-3 10^7])
%set(gca,'XTick',round(elevmu))
set(gca,'YTick',[10^-3 10^0 10^3 10^6])


%% Now plug in %rod model

pupDiDom = [.5 1 2];


%xhat = [-.0067 -.3046];
param = [.5863 .0722];


subplot(3,1,3)

for i = 1:length(pupDiDom)
    
    converter = pupDiDom(i)^2/pupilDi^2;
    
    Isom(:,i) = IsomerizationRate_rods*converter;
    
    percRod(:,i) = exp(-param(2)*(Isom(:,i)/1000).^param(1));
    
    plot(LumDom,percRod(:,i),'MarkerSize',20)
    hold on
    
    legstr{i} = [num2str(pupDiDom(i)) ' mm'];

end

xlabel('Luminance (cd/m^2)')
ylabel('%rods')
legend(legstr)


%set(gca,'XTick',round(elevmu))
ylim([0 1])


% pRods = rodIsom2percentRods(IsomerizationRate_rods);
%
% figure,
% semilogy(elevmu,pRods,'o-');
% xlabel('solar elevation')
% 
% set(gca,'XTick',round(elevmu))


%%
% dom = lambda;
% 
% RGB = [green' UV'];
% %T = lms'*RGB*(dom(2)-dom(1));  %Total Watts/steradian/m2
% 
% 
% pupilDi = 2;
% eyeDi = 3.4;
% 
% Wperm2 = RGB * (pupilDi/eyeDi)^2;  %flux at the retina (W/m^2/dlambda)
% Photperm2perlam = Wperm2.*([dom' dom']*10^-9)/(3*10^8)/(6.63*10^-34); %flux at retina (photons/mm^2/dlambda)
% %Photperm2perlam = Wperm2.*500*(10^-9)/(3*10^8)/(6.63*10^-34);
% 
% Photperm2 = Photperm2perlam*(dom(2)-dom(1))/(10^12); %phot/um^2
% 
% IsomerizationRate = sum(Sens'*Photperm2,2) %R*/sec/receptor for M cones, S cones, and rods 

