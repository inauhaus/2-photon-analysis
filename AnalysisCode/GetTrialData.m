function [CHs f1 f2] = GetTrialData(chvec,tcr,varargin)

%example: GetTrialData([1 1 0 0],[cond rep]); or GetTrialData([1 1 0 0],trial)

%if tcr has one element, then it is the trial no.  if it has two
%elements then it is the cond/repeat

global twophDATADIR ACQinfo


if length(tcr) == 1
    trial = tcr;
elseif length(tcr) == 2
    cond = tcr(1);
    rep = tcr(2);
    trial = gettrial(cond,rep);
end

%TTL0 is configured to generate an event on the rising edge 
%while TTL1 is configured to generate an event on either edge.
%event_id = 3 means start of trial (rising edge to TTL0 and TTL1); 
%event_id = 2 means end of trial (falling edge to TTL0 and TTL1);

if ACQinfo.SBInfo.scanbox_version==3
trialstartID = find(ACQinfo.SBInfo.event_id == 1); %Find events corresponding to rising edge from Stimulus
%RTO 8/17/21- changed from 3(rise) 2 (fall) for sbx 3.0
firstFrames = ACQinfo.SBInfo.frame(trialstartID);
trialendID = find(ACQinfo.SBInfo.event_id == 2); %Find events corresponding to falling edge from Stimulus
lastFrames = ACQinfo.SBInfo.frame(trialendID);
elseif ACQinfo.SBInfo.scanbox_version==2
    trialstartID = find(ACQinfo.SBInfo.event_id == 3); %Find events corresponding to rising edge from Stimulus
    %RTO 8/17/21- changed from 3(rise) 2 (fall) for sbx 3.0
    firstFrames = ACQinfo.SBInfo.frame(trialstartID);
    trialendID = find(ACQinfo.SBInfo.event_id == 2); %Find events corresponding to falling edge from Stimulus
    lastFrames = ACQinfo.SBInfo.frame(trialendID);
else 
    disp('ERROR: Unknown Scanbox version')
end

%This conditional indicates that a falling edge started the experiment for
%some reason... Maybe last experiment was aborted?
% if ACQinfo.SBInfo.event_id(1) == 2 & length(trialendID)>length(trialstartID)
%     trialendID = trialendID(2:end);
%     lastFrames = lastFrames(2:end);
%     'Extra sync.  Trying to correct this'
% end

f1 = firstFrames(trial);
f2 = lastFrames(trial);

if isempty(varargin)
    Nf = f2-f1+1;
else %Sometimes, we want to specify a different end time than the falling edge
    T = varargin{1};
    Nf = round(T/(ACQinfo.msPerLine*ACQinfo.linesPerFrame));
    f2 = f1+Nf-1; %redefine for the output.
end

imgs = single(sbxread(twophDATADIR,f1,Nf));
dim = size(imgs);

if chvec(1)
    CHs{1} = squeeze(imgs(1,:,:,:));    
    CHs{1} = CHs{1}(:,ACQinfo.unblanked,:);
end

if chvec(2)
    if dim(1) == 2
        CHs{2} = squeeze(imgs(2,:,:,:));
        CHs{2} = CHs{2}(:,ACQinfo.unblanked,:);
    else
        CHs{2} = 0;
    end
end


