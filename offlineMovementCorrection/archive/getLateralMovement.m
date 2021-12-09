function getLateralMovement

global ACQinfo 

nc = getnoconditions;

NT = getnotrials;

W = 4;
for t = 1:nt
    localT = (t-W/2):(t+W/2);
    idbad = find(localT > nT | localT < 1);
    localT(idbad) = [];
    
    for i = 1:length(localT)
        CHs = GetTrialData([1 1 0 1],localT(i)); 
        getTTLmovie
    end
    
end

