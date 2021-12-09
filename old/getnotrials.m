function nt = getnotrials

global looperInfo

nc = getnoconds;

nt = 0;
for c = 1:nc
    nr = length(looperInfo.conds{1}.repeats);
    nt = nt+nr;
end