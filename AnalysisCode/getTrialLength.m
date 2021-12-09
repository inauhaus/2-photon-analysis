function Nf = getTrialLength(trial)

global ACQinfo

trialstartID = find(ACQinfo.SBInfo.event_id == 3); %Find events corresponding to rising edge from Stimulus
firstFrames = ACQinfo.SBInfo.frame(trialstartID);
trialendID = find(ACQinfo.SBInfo.event_id == 2); %Find events corresponding to rising edge from Stimulus
lastFrames = ACQinfo.SBInfo.frame(trialendID);

%This conditional indicates that a falling edge started the experiment for
%some reason... Maybe last experiment was aborted?
if ACQinfo.SBInfo.event_id(1) == 2 & length(trialendID)>length(trialstartID)
    trialendID = trialendID(2:end);
    lastFrames = lastFrames(2:end);
    'Extra sync.  Trying to correct this'
end

f1 = firstFrames(trial);
f2 = lastFrames(trial);
Nf = f2-f1+1;