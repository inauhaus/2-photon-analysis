function IsomerizationRate = getIsomerizationRate(anim,unit,pupilDi)

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

[UVfile Greenfile] = getRawDataLocation(anim,unit)

% Load power measurements (%Units are Watts/steradian/m2 from spectroradiometer)
%pupilDi = 2;

%root = '/Users/in2293/Documents/MATLAB/Stimulator_master/Calibration/UVprojector/06_07_18/';
UVfile = [UVfile 'spectrum_UV128'];   
load(UVfile,'I','dom')
UV = I;
%root = '/Users/in2293/Documents/MATLAB/Stimulator_master/Calibration/RGBprojector/06_07_18/';
RGBfile = [Greenfile 'spectrum_Green128'];
load(RGBfile,'I','dom')
green = I;  %W/steradian/micron2/dlambda

%Compute luminance%%%%%%%%%%%%%%%%%
Vlam = 1.019*exp(-285.4*(dom(:)/1000-0.559).^2);
Km = 683;
specRad = green(:);
Luminance = Km*(Vlam'*specRad)*(dom(2)-dom(1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get absorption functions (Govardovskii)
lambda = 380:700;
[rod,Mcone,Scone] = photoreceptor_templates(500,508,360,lambda);
id = find(lambda>=dom(1) & lambda<=dom(end));
lambda = lambda(id);
Scone = Scone(id);
Mcone = Mcone(id);
rod = rod(id);




% Scale absorption functions by end-on collection area at the peak (um^2) (see Naarendorp and Wang et al)
Sens = [Mcone'*1 Scone'*1 rod'*0.85]; 

% Interpolate power measurements.  UV has a sharp peak that we don't want to miss.
UV   = interp1(dom,UV,lambda,'spline'); 
green = interp1(dom,green,lambda,'spline');

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

UV = UV.*T/100;
green = green.*T/100;

% Convert to flux at the retina
%pupilDi = 2;
eyeDi = 3.4; % mm (Stone & Pinto 93; Lyubursky 2004)

Apupil = pi*(pupilDi/2)^2;


%Convert spectroradiometer measure (W/steradian/m2/dlambda)
%Steradians = Apupil/monitorDistance^2.  Multiplying by this gives the flux through the pupil per um2 off the display
%Areal image magnification is eyeDi^2/monitorDistance^2
%To get flux/um/dlambda, where um is on the retina, we multiply by
%steradians/magnification, or Apupil/eyeDi^2.

converter = Apupil/(eyeDi)^2  %unitless

%Wperum2 = [green(:) UV(:)] * (pupilDi/eyeDi)^2 / (10^12);  %flux at the retina (W/um^2/dlambda)

Wperm2 = [green(:) UV(:)] * converter; %flux at the retina (W/m^2/dlambda)

Wperum2 = Wperm2 / (10^12);  %flux at the retina (W/um^2/dlambda)


%Convert to quanta
Photperum2perlam = Wperum2.*([lambda' lambda']*10^-9)/(3*10^8)/(6.63*10^-34); %convert to quanta (photons/sec/um^2/dlambda)

%figure,plot(lambda,Photperum2perlam)
%ylabel('photons/sec/um2/dlambda')

dlam = lambda(2)-lambda(1);
IsomerizationRate = sum(Sens'*Photperum2perlam*dlam,2) %R*/sec/receptor for M cones, S cones, and rods 

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


function [UVfile Greenfile] = getRawDataLocation(anim,unit)

%Calibrated 4/17/17
%Keynote UV and 100% offset RGB.
%UV current at 250. Green current at 160.
root.a.UV = 'C:\Stimulator_master\Calibration\UVprojector\4_17_17\';
%root.a.UV = [root.a.UV 'spectrum_UV128'];
root.a.green = 'C:\Stimulator_master\Calibration\RGBprojector\4_17_17\';
%root.a.green = [root.a.green 'spectrum_G128'];
ID.nz4_U2 = 'a';

%Calibrated 6/23/17
%UV projector with peak at 385. This is the UV projector from Keynote (Dallas).
%Moved projectors about 8" closer to make it brigher. Screen size is now 25x15.5.
%Teflon is .01" thick.-UV current set to 100. Green current set to 64
%In 1.164 (ISI) room. Using the 100% offset RGB for the green.
root.b.UV = 'C:\Stimulator_master\Calibration\UVprojector\6_23_17\';
%root.b.UV = [root.a.UV 'spectrum_UV250'];
root.b.green = 'C:\Stimulator_master\Calibration\RGBprojector\6_23_17\';
%root.b.green = [root.a.green 'spectrum_G160'];
ID.ra5_U2 = 'b';
ID.ra6_U1 = 'b';
ID.ra5_U3 = 'b';
ID.ra6_U2 = 'b';

%Calibrated 7/14/17
%Old EKB Projector, 405 LED. Fly eye burnt on the Keynote Projector
%UV current 128, Green Current 64. 
%0.01" Teflon. Green 100% offset
root.c.UV = 'C:\Stimulator_master\Calibration\UVprojector\7_14_17\';
root.c.green = 'C:\Stimulator_master\Calibration\RGBprojector\7_14_17\';
ID.rc2_U2 = 'c';
ID.rf7_U0 = 'c';
ID.rf5_U0 = 'c'; 
ID.rf6_U0 = 'c';
ID.rf4_U0 = 'c';
ID.rf3_U0 = 'c';
ID.rf9_U0 = 'c';
ID.rb8_U2 = 'c';
ID.rc4_U5 = 'c';

%Calibrated 10/16/17. 
%Old EKB Projector, 405 LED. 
%UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. 
%Projector was pushed up by 21.8cm new screen size is 46.7 by 29.8 cm
root.d.UV = 'C:\Stimulator_master\Calibration\UVprojector\10_16_17\';
root.d.green = 'C:\Stimulator_master\Calibration\RGBprojector\10_16_17\';
ID.rc3_U6 = 'd'; 
ID.rk4_U2 = 'd'; 
ID.rk5_U3 = 'd'; 
ID.rl0_U3 = 'd'; 
ID.rl1_U3 = 'd'; 
ID.rk6_U1 = 'd'; 
ID.rk7_U2 = 'd'; 
ID.rl5_U4 = 'd'; 
ID.rl7_U3 = 'd'; 
ID.rm0_U2 = 'd'; 
ID.rm8_U2 = 'd'; 
ID.rm8_U3 = 'd'; 
ID.rn1_U3 = 'd'; 
ID.rm2_U3 =  'd';
ID.rm2_U4 =  'd';


%Calibrated 3/01/18.
%Old EKB Projector, 405 LED.
%UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset.
%Sanity check recalibration.
root.e.UV = 'C:\Stimulator_master\Calibration\UVprojector\03_01_18\';
root.e.green = 'C:\Stimulator_master\Calibration\RGBprojector\03_01_18\';
ID.rm1_U1 = 'e';

%Calibrated 3/22/18
%Old EKB Projector, 405 LED.
%UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. 
%Recalibrated because new mobile screen was attached to the projector. Screen is now 1 inch farther than it used to be.
%Recalibrated because green projector was pushed up closer to the UV projector on the bretboard
root.f.UV = 'C:\Stimulator_master\Calibration\UVprojector\03_22_18\';
root.f.green = 'C:\Stimulator_master\Calibration\RGBprojector\03_22_18\';
ID.rm1_U3 = 'f';

%Calibrated 3/26/18. 
%Old EKB Projector, 405 LED.
%UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. 
%Recalibrated because projectors were moved by allison, unclear what the movement was...
root.g.UV = 'C:\Stimulator_master\Calibration\UVprojector\03_26_18\';
root.g.green = 'C:\Stimulator_master\Calibration\RGBprojector\03_26_18\';
ID.rm5_U3 = 'g'; 
ID.rm0_U3 = 'g'; 
ID.rl5_U7 = 'g'; 
ID.rl9_U5 = 'g'; 
ID.rm4_U4 = 'g'; 
ID.rm5_U4 = 'g'; 
ID.rr0_U2 = 'g'; 
ID.rr1_U2 = 'g'; 
ID.rr5_U1 = 'g'; 
ID.rm4_U5 = 'g'; 
ID.rr1_U3 = 'g'; 
ID.rr8_U2 = 'g'; 
ID.rs9_U2 = 'g'; 
ID.rw3_U1 = 'g'; 
ID.rw3_U2 = 'g'; 
ID.rt0_U2 = 'g'; 
ID.rt2_U2 = 'g'; 
ID.rs6_U2 = 'g'; 
ID.rs9_U3 = 'g'; 

% Calibrated 6/07/18. 
% Old EKB Projector, 405 LED. 
% UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. Buffer values are both 128. 
% Recalibrated as a sanity check.
root.h.UV = 'C:\Stimulator_master\Calibration\UVprojector\06_07_18\';
root.h.green = 'C:\Stimulator_master\Calibration\RGBprojector\06_07_18\';
ID.rw5_U3 = 'h';
ID.rw6_U3 = 'h';


%Calibrated 7/23/18. 
%Old EKB Projector, NEW!!! 405 LED. New EKB fly eye, and condensor from Dallas projector. 
%Replaced burnt fly eye from 7/06/18. 0.01" Teflon. Green 100% offset. 
%Buffer values both at 128';
%Tmatrix values [0.1823 0.0898; 0.0013 0.0587; 0.1600 0.0985];  %outer product of cones x gun spectra in dichromat: [Mcone Scone rod]'*[Green UV];
root.i.UV = 'C:\Stimulator_master\Calibration\UVprojector\7_23_18\';
root.i.green = 'C:\Stimulator_master\Calibration\RGBprojector\07_23_18\';

ID.sb6_U2 = 'i';


% Calibrated 9/20/17
%2p room
% .01" Teflon
% Proj 43cm from screen (different from previous; moved closer towards screen)
% RGB at 128 current
% 0% offset Green/RGB projector UV at 128 current
% Matrix values [0.2013 0.1233; 0.0002 0.0517; 0.1763 0.1375] %buffer values g128, uv128
root.j.UV = 'C:\Stimulator_master\Calibration\UVprojector\9_20_17\';
root.j.green = 'C:\Stimulator_master\Calibration\RGBprojector\9_20_17\';

ID.rk3_U1 =  'j';
ID.rk4_U1 =  'j';
ID.rk5_U1 =  'j';
ID.rk7_U1 =  'j';
ID.rk3_U2 =  'j';
ID.rl0_U1 =  'j';
ID.rl2_U0 =  'j';
ID.rl1_U2 =  'j';
ID.rk9_U1 =  'j';
ID.rl0_U2 =  'j';
ID.rk5_U2 =  'j';
ID.rm0_U0 =  'j';
ID.rl9_U1 =  'j';
ID.rl8_U1 =  'j';
ID.rl9_U2 =  'j';
ID.rl7_U1 =  'j';

% Calibrated 1/5/2018
%2p
% .01" Teflon
% Proj 43cm from screen (different from previous; moved closer towards screen)
% RGB at 244 current UV at 250 current
% 0% offset Green/RGB projector
%Baselines 6,5,4,3,2,1,0.5 were all run using the values below, B .25 was
%at buffer 128 for green and UV.
%Tmatrix Values [0.5680 0.3576; 0.0002 0.1393; 0.4961 0.4003] %baseline .25 at buffer 128 for green and UV 
root.k.UV = 'C:\Stimulator_master\Calibration\UVprojector\1_5_18\'; % For UV250
root.k.green = 'C:\Stimulator_master\Calibration\RGBprojector\1_5_18\'; % For RGB244
%
root.k2.UV = 'C:\Stimulator_master\Calibration\UVprojector\1_2_18\'; % For UV128
root.k2.green = 'C:\Stimulator_master\Calibration\UVprojector\1_2_18\'; % For RGB128

%.25 baseline comparison: UV128/UV250 = 0.5115 RGB128/RGB244 = 0.5083

ID.rl5_U1 =  'k'; % This doesn’t have B .25
ID.rl7_U2 =  'k'; % This doesn’t have B .25
ID.rm4_U1 =  'k';
ID.rm3_U0 =  'k';
ID.rn2_U0 =  'k';
ID.rm2_U0 =  'k'; % Start of B .25 experiments
ID.rl7_U2 =  'k';
ID.rl5_U1 =  'k';



% Calibrated 1/15/18. 
%2p
%.01" Teflon, 43cm from screen, same UV and RGB projectors from 7/13/17 calibration. 
%Green at 244, UV at 250 current
%Baselines 6,5,4,3,2,1,0.5 were all run using the values below, B .25 was
%at buffer 128 for green and UV.
%Tmatrix values [0.5660 0.3453; 0.0002 0.1312; 0.4938 0.3870] %baseline .25 at buffer 128 for green and UV

% root.l.UV = 'C:\Stimulator_master\Calibration\UVprojector\1_15_18\'; % For UV250
% root.l.green = 'C:\Stimulator_master\Calibration\RGBprojector\1_15_18\'; % For RGB244
%
root.l.UV = 'C:\Stimulator_master\Calibration\UVprojector\1_15_18\'; % For UV128
root.l.green = 'C:\Stimulator_master\Calibration\RGBprojector\1_15_18\'; % For RGB128 

%Baseline .25 Comparison of UV128/UV250 = 0.5149 RGB128/RGB244 = 0.5154

ID.rm5_U1 =   'l';
ID.rn1_U1 =  'l';
ID.rm5_U2 =  'l';
ID.rm4_U2 =  'l';
ID.rm3_U1 =  'l';
ID.rm3_U2 =  'l';
ID.rm6_U1 =  'l';
ID.rm2_U1 =  'l';
ID.rm2_U2 =  'l';
ID.rn2_U1 =  'l';
ID.rn5_U1 =  'l';
ID.rn1_U2 =  'l';
ID.rn4_U2 = 'l';
ID.rr1_U1 = 'l';
ID.rr0_U1 = 'l';
ID.rr2_U1 = 'l';
ID.rr9_U1 = 'l';
ID.rr8_U1 = 'l';

% Calibrated 2/27/18
%2p
%.01" Teflon, 43cm from screen, same UV and RGB projectors from 7/13/17 calibration but UV LED swapped for 'UV405-keynote' from the Keynote (Dallas?) projector.
%Green buffer at 220, UV buffer at 225 current
%Baselines 6,5,4,3,2,1,0.5 were all run using the values below, B .25 was
%at buffer 128 for green and UV.
%Tmatrix values  [0.3430 0.1859; 0.0002 0.0800; 0.3003 0.2070]; % baseline 0.25 current values are 128 for green and UV
root.m.UV = 'C:\Stimulator_master\Calibration\UVprojector\2_27_18\';
root.m.green = 'C:\Stimulator_master\Calibration\RGBprojector\2_27_18\';

%Baseline .25 comparison: UV128/UV225 = 0.5656  RGB128/RGB220 = 0.5686

ID.rr3_U1 = 'm';
ID.rs7_U1 = 'm';
ID.rs6_U1 = 'm';
ID.rt2_U1 = 'm';


%%



unit = num2str(unit);
flds = fields(ID);
for i = 1:length(flds)
   
    animdum = flds{i}(1:3);
    unitdum = flds{i}(6);
    if strcmp(anim,animdum) & strcmp(unit,unitdum)
        letterID = eval(['ID.' anim '_U' unit]);
        break
    end
    
end

UVfile = eval(['root.' letterID '.UV;']);
Greenfile = eval(['root.' letterID '.green;']);

