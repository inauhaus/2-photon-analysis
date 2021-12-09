function parseToCellS(cellTraces)

global cellS ACQinfo


trialstartID = find(ACQinfo.SBInfo.event_id == 3); %Find events corresponding to rising edge from Stimulus
firstFrames = ACQinfo.SBInfo.frame(trialstartID);
trialendID = find(ACQinfo.SBInfo.event_id == 2); %Find events corresponding to falling edge from Stimulus
lastFrames = ACQinfo.SBInfo.frame(trialendID);
Nmin = min(lastFrames-firstFrames);


nc = getnoconditions;
nr = getnorepeats(1);

cellS.cellMat = cell(1,nc);

Nshort = 10^10;
%for c = 1:nc
for c = 1:15
    for r = 1:nr
        trialBlock = GetTrialData_fromCellTraces(cellTraces,[c r]);
        trialBlock = trialBlock(1:Nmin,:);
        cellS.cellMat{c}(:,:,r) = trialBlock';
    end
end
