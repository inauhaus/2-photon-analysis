function ttag = gettimetag(c,r)

global looperInfo

ttag = looperInfo.conds{c}.repeats{r}.timetag;
