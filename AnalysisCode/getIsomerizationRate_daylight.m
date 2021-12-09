function [percRod Isom elevmu pupDiDom] = getIsomerizationRate_daylight(pupilDi)

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


load('C:\2pScanboxAnalysis\mouse_color\solarIrradianceData')
%%

CIE_D65 = illuminant('d65',Rural.nm);  %daylight standard.  

clear PSDmu elevmu PSDmu_norm
elevationBinEdges = [-40:15:85];
%elevationBinEdges = [-30:15:90];
for i = 1:length(elevationBinEdges)-1
    
    id = find(Rural.SolarElev > elevationBinEdges(i) & Rural.SolarElev < elevationBinEdges(i+1));
    
    PSDmu(:,i) = mean(Rural.PSD(:,id),2);
    PSDmu_norm(:,i) = PSDmu(:,i)/max(PSDmu(:,i));
    
    elevmu(i) = mean(Rural.SolarElev(id));
    
end
dom = Rural.nm;

figure,plot(dom,PSDmu)

id455 = find(dom == 455);
for i = 1:size(PSDmu,2)  %loop each elevation
    if isnan(PSDmu(1,i))  %check to see if they made the measurement at low wavelengths
        PSDmu(:,i) = PSDmu(id455,i)*CIE_D65/CIE_D65(id455);
    end
    PSDmu_norm(:,i) = PSDmu(:,i)/max(PSDmu(:,i));
end
hold on
plot(dom,PSDmu)



figure,
subplot(2,1,1), 
semilogy(dom,PSDmu), ylabel('Watts/m^2'), xlabel('nm')
xlim([340 800])
subplot(2,1,2),
plot(dom,PSDmu_norm), ylabel('linear unitless'), xlabel('nm')
xlim([340 800])

%% For paper figure; plot the night and day spectra
% id = find(dom>300 & dom<750);
% figure, 
% plot(dom(id),PSDmu(id,2)/max(PSDmu(id,2)))
% hold on
% plot(dom(id),PSDmu(id,end)/max(PSDmu(id,end)))

%% Compute luminance at each elevation, from irradiance values

clear Lum
specRad = PSDmu/pi;  %convert to radiance W/m2/steradian, from irradiance. 
%I don't know for sure if this conversion is correct.  Its based on
%equating my two equations for retinal irradiance: 
%Irr_outside*Apupil/Aretina = Rad_outside*Apupil/(2*radius_retina)^2, gives
%Rad_outside = Irr_outside/pi
%Literature says daylight is 10^5, so these values seem low.  But these are
%"average luminance" across direction.  
for i = 1:size(PSDmu,2)
      
    Vlam = 1.019*exp(-285.4*(dom(:)/1000-0.559).^2);
    Km = 683;
    
    
    specRaddum = specRad(:,i);
    
    Lum(i) = Km*(Vlam'*specRaddum)*(dom(2)-dom(1))
    
end
figure, semilogy(elevmu,Lum)
xlabel('solar elevation')
ylabel('luminance')
set(gca,'YTick',[10^-3 10^0 10^3 10^6])

%%

% Get absorption functions (Govardovskii)
lambda = 300:750;
id = find(lambda>=dom(1) & lambda<=dom(end));
lambda = lambda(id);

[rod,Mcone,Scone] = photoreceptor_templates(500,508,360,lambda);


% Scale absorption functions by end-on collection area at the peak (um^2) (see Naarendorp and Wang et al)
Sens = [Mcone'*1 Scone'*1 rod'*0.85]; 

% Interpolate power measurements.  UV has a sharp peak that we don't want to miss.
PSDmu_I   = interp1(dom',PSDmu,lambda','spline'); 

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

PSDmu_I = PSDmu_I.*(T'*ones(1,size(PSDmu,2)))/100;

% Convert to flux at the retina
%pupilDi = 2;
eyeDi = 3.4; % mm (Remtulla)

Aretina = 4*pi*(eyeDi/2)^2;
Apupil = pi*(pupilDi/2)^2;


converter = Apupil/Aretina;  %unitless

%Wperum2 = [green(:) UV(:)] * (pupilDi/eyeDi)^2 / (10^12);  %flux at the retina (W/um^2/dlambda)

Wperm2 = PSDmu_I * converter; %flux (irradiance) at the retina (W/m^2/dlambda)

Wperum2 = Wperm2 / (10^12);  %flux at the retina (W/um^2/dlambda)


%Convert to quanta
lambdaMat = lambda(:)*ones(1,size(PSDmu,2))*10^-9;  %units of meters
Photperum2perlam = (Wperum2.*lambdaMat)/(3*10^8)/(6.63*10^-34); %convert to quanta (photons/sec/um^2/dlambda)

%figure,plot(lambda,Photperum2perlam)
%ylabel('photons/sec/um2/dlambda')

dlam = lambda(2)-lambda(1);
IsomerizationRate_Mopsin = Sens(:,1)'*Photperum2perlam*dlam; %R*/sec/receptor for M cones, S cones, and rods 
IsomerizationRate_Sopsin = Sens(:,2)'*Photperum2perlam*dlam; %R*/sec/receptor for M cones, S cones, and rods 
IsomerizationRate_rods = Sens(:,3)'*Photperum2perlam*dlam; %R*/sec/receptor for M cones, S cones, and rods 

figure,

subplot(3,1,1),
semilogy(elevmu,max(PSDmu),'o-'), xlabel('solar elevation'), ylabel('peak irradiance (W/m^2)')

subplot(3,1,2)
semilogy(elevmu,[IsomerizationRate_Mopsin' IsomerizationRate_Sopsin' IsomerizationRate_rods'] ,'.-');
xlabel('solar elevation')

hold on
plot([elevmu(1) elevmu(end)],[10*10^3 10*10^3 ],'--k')
hold on
plot([elevmu(1) elevmu(end)],[400*10^3 400*10^3 ],'--k')

ylim([10^-3 10^7])
%set(gca,'XTick',round(elevmu))
set(gca,'YTick',[10^-3 10^0 10^3 10^6])

ylabel('R*/receptor/sec')
legend('Mopsin','Sopsin','rods')

% %% Now plug in %rod model
% 
% subplot(3,1,3)
% 
% %xhat = [-.0067 -.3046];
% param = [.5863 .0722]
% percRod = exp(-param(2)*(IsomerizationRate_rods/1000).^param(1));
% 
% plot(elevmu,percRod,'.-k','MarkerSize',20)
% xlabel('solar elevation (deg)')
% ylabel('%rods')
% 
% %set(gca,'XTick',round(elevmu))
% ylim([0 1])


%% Now plug in %rod model

pupDiDom = [.5 1 2];


%xhat = [-.0067 -.3046];
param = [.5863 .0722];


subplot(3,1,3)

for i = 1:length(pupDiDom)
    
    converter = pupDiDom(i)^2/pupilDi^2;
    
    Isom(:,i) = IsomerizationRate_rods*converter;
    
    percRod(:,i) = exp(-param(2)*(Isom(:,i)/1000).^param(1));
    
    plot(elevmu,percRod(:,i),'.-')
    hold on
    
    legstr{i} = [num2str(pupDiDom(i)) ' mm'];

end
    
xlabel('solar elevation (deg)')
ylabel('%rods')
legend(legstr)


%set(gca,'XTick',round(elevmu))
ylim([0 1])


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

