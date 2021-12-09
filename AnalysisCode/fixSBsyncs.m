function fixSBsyncs

global ACQinfo

%% If there are more than 2^16 frames, scanbox cycles back to zero. This
%unwraps it...
if ~isempty(find(diff(ACQinfo.SBInfo.frame)<0))
    id = find(diff(ACQinfo.SBInfo.frame)<0); %Indicates last frame before returning to zero    
    ACQinfo.SBInfo.frame(id+1:end) = ACQinfo.SBInfo.frame(id+1:end)+(2^16-1);
end

%% Sometimes the first sync indicates a falling edge, and its also an extra one:

trialstartID = find(ACQinfo.SBInfo.event_id == 3); %rising edge ID
trialendID = find(ACQinfo.SBInfo.event_id == 2); %falling edgeID
if ACQinfo.SBInfo.event_id(1) == 2 & length(trialendID)>length(trialstartID)
    
    ACQinfo.SBInfo.event_id(1) = [];
    ACQinfo.SBInfo.line(1) = [];
    ACQinfo.SBInfo.frame(1) = [];
    
    %lastFrames = lastFrames(2:end);
    'Falling edge sync at beginning.  Trying to correct this'
end

% 
%% Sometimes it does not read a rising edge as a synchronous event on
% TTL0 and TTL1. The stimulus TTL (which is high during the
% stimulus) is connected to both TTL0 and TTL1.  imask = 3.  I use event_id
% of 3 and 2 to get the beginning and end of the trial, respectively.
% Occassionally (~1/1000 trials), the eventid = 3 gets replaced with two
% eventids of 1 and 2 that are very close in time... ~ 1 single scan line
% apart... i.e. they id the same frame but not the same line.
% 
% In other words,TTL0 and TTL1 are receiving the same input, but its not
% registering the rising edge as occuring simultaneously... on rare
% occasion.

Nev1 = length(find(ACQinfo.SBInfo.event_id == 1)); 
Nev2 = length(find(ACQinfo.SBInfo.event_id == 2));
Nev3 = length(find(ACQinfo.SBInfo.event_id == 3));

if (Nev2-Nev1 == Nev3+Nev1) && (Nev3+Nev1 == getnotrials) %My logic behind what is happening only makes sense if this is true
    
    'WTF, messed up syncs. There are event_id(s) == 1'
    
    idR = find(diff(ACQinfo.SBInfo.frame)==0);
    id1 = find(ACQinfo.SBInfo.event_id == 1);
    
    if ~isempty(id1) & ~isempty(idR)
        for i = length(id1):-1:1  %idR and id1 should be the same length. Work backwards so that indices don't change with the update in each iteration
            RedundFr = [idR(i) idR(i)+1]; %event_id's that have same frame number
            idx = find(id1(i) == RedundFr); %One of the two Redundant Frames elements should be coincident with eventid = 1
            removeID = RedundFr(idx); %event id that should be removed
            changeID = RedundFr(3-idx); %event id that should be changed from 2 to 3
            if ACQinfo.SBInfo.event_id(changeID) == 2 %My logic behind what is happening only makes sense if this is true
                ACQinfo.SBInfo.event_id(changeID) = 3;
                
                ACQinfo.SBInfo.event_id(removeID) = [];
                ACQinfo.SBInfo.line(removeID) = [];
                ACQinfo.SBInfo.frame(removeID) = [];
                
            end
            
        end
    end
    
end


