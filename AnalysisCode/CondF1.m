function F1image = CondF1(Flim)

global ACQinfo Tens

Flim = getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning
Flim(2) = Flim(1) + 1000*getparam('stim_time');

tdom = (0:size(Tens{1},3)-1)*ACQinfo.msPerLine*ACQinfo.linesPerFrame;
[dum Flim(1)] = min(abs(tdom-Flim(1)));
[dum Flim(2)] = min(abs(tdom-Flim(2)));

tdom = tdom-tdom(Flim(1));

try
    STIMfrate = Analyzer.framerate; %ms/period
catch
    STIMfrate = 60;
    'Frame rate does not exist in Analyzer. Setting to 60Hz'
end
T = getParamVal('t_period')/STIMfrate*1000; %ms/cycle
ce = exp(1i*2*pi*tdom/T);

maxCycle = max(unique(floor(tdom/T)));
maxID = find(tdom/T>=maxCycle);
maxID = maxID(1);
Flim(2) = maxID;

for c = 1:length(Tens)
    c
    muImage = mean(Tens{c},3);
    
    F1image{c} = 0;
    for i = Flim(1):Flim(2)
        F1image{c} = F1image{c} + (Tens{c}(:,:,i)-muImage)*ce(i);
    end
    F1image{c} = (F1image{c})/length(Flim(1):Flim(2)); %peak amplitude

end
