function OnePeriod = getOnePeriod

global ACQinfo cellS Analyzer

Flim = getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning
Flim(2) = Flim(1) + 1000*getparam('stim_time');

tdom = (0:size(cellS.cellMat{1},2)-1)*ACQinfo.msPerLine*ACQinfo.linesPerFrame;
[dum Flim(1)] = min(abs(tdom-Flim(1)));
[dum Flim(2)] = min(abs(tdom-Flim(2)));


try
    STIMfrate = Analyzer.framerate; %ms/period
catch
    STIMfrate = 60;
    'Frame rate does not exist in Analyzer. Setting to 60Hz'
end

T = getParamVal('t_period')/STIMfrate*1000; %ms/cycle

tdom = tdom - tdom(Flim(1));  %Need phase to start at zero


maxCycle = max(unique(floor(tdom(1:Flim(2))/T)));
maxID = find(tdom/T>=maxCycle);
maxID = maxID(1);
Flim(2) = maxID;


tdom = tdom(Flim(1):Flim(2));

phi = 2*pi*tdom/T;

phiWrap = angle(exp(1i*phi))*180/pi;
id = find(phiWrap<0);
phiWrap(id) = phiWrap(id)+360;
phiWrap = round(phiWrap/5)*5;
phidom = unique(phiWrap);


clear OnePeriod
for c = 1:length(cellS.cellMat)
    Rmat = cellS.cellMat{c}(:,Flim(1):Flim(2),:);
    Rmat = mean(Rmat,3); %Mean over repeats
    mu = mean(Rmat,2);
    Rmat = Rmat-mu*ones(1,length(tdom));
    Rmat = Rmat./(mu*ones(1,length(tdom)));

    for i = 1:length(phidom)
        
        id = find(phiWrap == phidom(i));
        OnePeriod{c}(i,:) = mean(Rmat(:,id),2);
        
    end
end

