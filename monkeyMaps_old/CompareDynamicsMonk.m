function [tt100 tt50 tt50back] = CompareDynamicsMonk(tcourseMat)

global G_RChandles ACQinfo DM

matDim = size(tcourseMat');

tauN = str2num(get(G_RChandles.kernelLength,'string'));

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);

taudom = DM.taudom;  %it will start at exactly kernDel(1) with acqPeriod spacing, and end at an estimate of kernDel(2)
taudomI = DM.taudom;

for i = 1:matDim(1)   %loop through each neuron
    
        tc = squeeze(tcourseMat(:,i));        
        
        %time to peak
        tcI = interp1(taudom,tc,taudomI,'spline');
        [ma idpk] = max(tcI);
        tt100(i) = taudomI(idpk); 
        
        
        %FRONT END
        tcdum = tcI(1:idpk);
        
        %time to 50%        
        thresh = (ma+tcI(1))*.5;
        [dum id] = min(abs(tcdum-thresh));
        tt50(i) = taudomI(id);
        
        %time to 10%
        thresh = (ma+tcI(1))*.1;
        [dum id] = min(abs(tcdum-thresh));
        tt10(i) = taudomI(id);       
        
        
        %BACK END
        tcdum = tcI(idpk:end);
        
        %time to 50%        
        thresh = (ma+tcI(1))*.5;
        [dum id] = min(abs(tcdum-thresh));
        id = id + idpk - 1;
        tt50back(i) = taudomI(id);
        
        %time to 10%
        thresh = (ma+tcI(1))*.1;
        [dum id] = min(abs(tcdum-thresh));
        id = id + idpk - 1;
        tt10back(i) = taudomI(id);    
     
        
end

id = find(tt10>500 | tt10<0);
tt10(id) = NaN;

id = find(tt50>400 | tt50<50);
tt50(id) = NaN;

id = find(tt50back>1800 | tt50back<200);
tt50back(id) = NaN;

id = find(tt100>600 | tt100<0);
tt100(id) = NaN;

FWHM = tt50back - tt50;
figure, hist(FWHM,[0:150:1400])
xlabel('FWHM (ms) ')
title(['mean = ' num2str(nanmean(FWHM))])
nanstd(FWHM)
nanmedian(FWHM)
%%
figure


subplot(3,1,1)
hist(tt100,[-100:100:600])
xlabel('change in delay @ peak (ms) ')
title(['mean = ' num2str(nanmean(tt100))])
subplot(3,1,2)
hist(tt50,[-100:100:600])
xlabel('change in delay @ 50% rise (ms) ')
title(['mean = ' num2str(nanmean(tt50))])
subplot(3,1,3)
hist(tt50back,[-100:100:1800])
xlabel('change in delay @ 50% fall (ms) ')
title(['mean = ' num2str(nanmean(tt50back))])

