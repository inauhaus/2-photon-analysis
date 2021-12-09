function CHs = Load2phImage(frame,trial)

%trial and frame number are 1 to N
%if you want to use the preset condition/repeat to define the trial, leave out varargin.


global twophDATADIR ACQinfo

%global tf filepath

if ACQinfo.SBInfo.scanbox_version==3
    trialstartID = find(ACQinfo.SBInfo.event_id == 1); %Find events corresponding to rising edge from Stimulus%rising edge ID
    %RTO 8/17/21- changed from 3(rise) 2 (fall) for sbx 3.0
elseif  ACQinfo.SBInfo.scanbox_version==2
    trialstartID = find(ACQinfo.SBInfo.event_id == 3);
else
    disp('ERROR: Unknown Scanbox version')
    
end
firstFrames = ACQinfo.SBInfo.frame(trialstartID) + 1; %frame indices of trial beginning.

fno = firstFrames(trial)+frame-1; 

imgs = sbxread([twophDATADIR],fno,1);
dim = size(imgs);

CHs{1} = squeeze(imgs(1,:,ACQinfo.unblanked));

if dim(1) == 2
    CHs{2} = squeeze(imgs(2,:,ACQinfo.unblanked));
else
    CHs{2} = 0;
end

