function experimentMotionAnalyzer

global maskS cellS

temp = maskS.anatomyTemplate;

mbest = 0;
nbest = 0;

cellS.motionInfo = struct;

for t = 1:getnotrials
    
    %Apply slow movement correction
    
    endT = getparam('stim_time')+getparam('predelay')+getparam('postdelay'); %Grab more than trial length
    CH_raw = GetTrialData([1 0],t,endT*1000);    
    CH = CH_raw{1};
    CH = smoothTensor(CH,1);

    for n = 1:size(CH_raw{1},3)
        imdum = CH(:,:,n);
        [mbest nbest] = getShiftVals(imdum.^2,temp.^2,[0 0]);  %squaring seems to really help sometimes
        CH(:,:,n) = circshift(CH(:,:,n),[round(-mbest) round(-nbest) 0]);
        
        cellS.motionInfo.mbest{t}(n) = mbest;
        cellS.motionInfo.nbest{t}(n) = nbest;
    end
    
    cellS.motionInfo.rTemp{t} = findBadframes(CH);
    
    t
end


function To = smoothTensor(To,tsig)

dim = size(To);
To = reshape(To,[dim(1)*dim(2) dim(3)])';

if tsig>0
    h = fspecial('gaussian',[dim(3) 1],tsig);
    h = h*ones(1,size(To,2));
    To = ifft(abs(fft(h)).*fft(To));
end

To = reshape(To',[dim(1) dim(2) size(To,1)]);
    