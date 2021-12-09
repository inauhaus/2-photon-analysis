function [lco hco] = gethcolco(dom,y,cothresh)

%get highpass cutoff
[dum idmax] = max(y);
dumtc = y(idmax:end); dumdom = dom(idmax:end);
dumtc = dumtc/dumtc(1);
if dumtc(end)>cothresh*1.05
    hco = NaN;
else
    [dum idco] = min(abs(dumtc-cothresh));
    hco = dumdom(idco); %high pass cutoff
end

%get lowpass cutoff
[dum idmax] = max(y);
dumtc = y(1:idmax); dumdom = dom(1:idmax);
dumtc = dumtc/dumtc(end);
if dumtc(1)>cothresh*1.05
    lco = NaN;
else
    [dum idco] = min(abs(dumtc-cothresh));
    lco = dumdom(idco); %low pass cutoff
end
            