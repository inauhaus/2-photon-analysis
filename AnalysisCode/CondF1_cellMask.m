function F1 = CondF1_cellMask

global ACQinfo cellS Analyzer

Flim = getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning
Flim(2) = Flim(1) + 1000*getparam('stim_time');

tdom = (0:size(cellS.cellMat{1},2)-1)*ACQinfo.msPerLine*ACQinfo.linesPerFrame;
[dum Flim(1)] = min(abs(tdom-Flim(1)));
[dum Flim(2)] = min(abs(tdom-Flim(2)));


try
    STIMfrate = Analyzer.framerate; %ms/period
catch
    STIMfrate = 60;
    'Frame rate does not exist in Analyzer. Setting to 60Hz'
end

T = getParamVal('t_period')/STIMfrate*1000; %ms/cycle

tdom = tdom - tdom(Flim(1));  %Need phase to start at zero


maxCycle = max(unique(floor(tdom(1:Flim(2))/T)));
maxID = find(tdom/T>=maxCycle);
maxID = maxID(1);
Flim(2) = maxID;



%%%%
%Flim(1) = 1;
%Flim(2) = length(tdom);
%%%%



tdom = tdom(Flim(1):Flim(2));

phi = 2*pi*tdom/T;

ce = exp(1i*phi);

fs = 1000/(ACQinfo.msPerLine*ACQinfo.linesPerFrame); %samples/sec
fdom = linspace(0,fs,length(tdom)+1); 
fdom = fdom(1:end-1);

[dum F1id] = min(abs(fdom-1000/T));

%Now correct for phase shift of the bar at t = 0 from perp bisector:
dxperp_cm = getparam('dx_perpbis');
dyperp_cm = getparam('dy_perpbis');
dxperp_deg = atan(dxperp_cm/Analyzer.M.screenDist)*180/pi %Shift in deg of visual field
dyperp_deg = atan(dyperp_cm/Analyzer.M.screenDist)*180/pi; %Shift in deg of visual field
sPerx = getparam('x_size');
sPery = getparam('y_size');
dphaseX = dxperp_deg/sPerx*360; 
dphaseY = dyperp_deg/sPery*360; 

%%

hannWin = hann(length(ce));
%ce = ce.*hannWin'; %Smooth in temporal frequency a bit
for c = 1:length(cellS.cellMat)
    
    ori = Analyzer.loops.conds{c}.val{1};
    
   if ori == 0
       ce_shift = ce*exp(1i*dphaseX*pi/180);
   elseif ori == 180
       ce_shift = ce*exp(-1i*dphaseX*pi/180);
   elseif ori == 90 
       ce_shift = ce*exp(1i*dphaseY*pi/180);
   elseif ori == 270
       ce_shift = ce*exp(-1i*dphaseY*pi/180);
   end
       
    
    dumMat = cellS.cellMat{c}(:,Flim(1):Flim(2),:);  %Grab chunk in stimulus window
    
    R{c} = mean(dumMat,3);
    
    muR{c} = mean(R{c},2);
    
    R{c} = R{c}-muR{c}*ones(1,length(tdom));
    
    

    
%     for p = 1:size(R{c},1)
%         
%         dum = squeeze(dumMat(p,:,:));
%         phimat = repmat(phi(:),[1 5]);
%         
%         [mu sigma ffit(p,:) varacc{c}(p) param{c}(p,:)] = CircGaussFit2(R{c}(p,:),rem(phi,2*pi));
%         %[mu sigma ffit(p,:) varacc{c}(p) param{c}(p,:)] = CircGaussFit2(dum(:),rem(phimat(:),2*pi));
%     end
% %     

    
    %F1{c} = R{c}*exp(1i*phi(:));
   
%     F1{c} = ffit*ce_shift(:);

%     Rw = fft(R{c}');
%     Rw(1:2,:) = 0;
%     Rw((F1id-1)*8:end,:) = 0;
%     R{c} = real(ifft(Rw)');

    F1{c} = R{c}*ce(:);  %Don't use transpose b/c it takes complex conjugate
    
    F1{c} = 2*F1{c}/sum(hannWin);
    
    Rw = fft(R{c}');
    %     Rw(1:F1id/2,:) = 0;
    %     Rw((F1id-1)*2:end,:) = 0;
    %     harmvec = F1id:(F1id-1):2*F1id+1;
    
    harmE = (abs(Rw(F1id,:)));
    W = floor(F1id/2);
    localID = [F1id-W:F1id-1 F1id+1:F1id+W];
    LocalE = abs(Rw(localID,:)); %Get local energy around the harmonic
    LocalEmu = mean(LocalE);
    LocalEsig = std(LocalE);
    SNR = (harmE - LocalEmu)./LocalEsig;
    
    id = find(SNR<5)
    F1{c}(id) = NaN;
    
 
    %F1F0 = abs(F1{c})./sum(R{c},[],2);
    
    
%     figure,plot(angle(F1{c}),angle(F1x{c}),'.')
%     hold on, plot([-pi pi],[-pi pi])
%     
%     F1{c} = exp(1i*param{c}(:,1));
     
end

