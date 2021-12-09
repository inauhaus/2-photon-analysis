function [c r] = getcondrep(ttag)

%ttag is 0 to N-1

global looperInfo

nc = length(looperInfo.conds);
nr = length(looperInfo.conds{1}.repeats);

for c = 1:nc
    for r = 1:nr
        ttagMat(c,r) = looperInfo.conds{c}.repeats{r}.timetag;
    end
end

[c r] = find(ttagMat == ttag);  %get cond and rep for this timetag