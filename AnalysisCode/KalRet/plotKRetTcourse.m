function plotKRetTcourse

global ACQinfo cellS Analyzer idExamp

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


fs = 1000/(ACQinfo.msPerLine*ACQinfo.linesPerFrame); %samples/sec
fdom = linspace(0,fs,length(tdom)+1); 
fdom = fdom(1:end-1);

[dum F1id] = min(abs(fdom-1000/T));

%%

for c = 1:length(cellS.cellMat)
    
    dumMat = cellS.cellMat{c}(:,Flim(1):Flim(2),:);  %Grab chunk in stimulus window
    
    R{c} = mean(dumMat,3);
    
end

%%
fr = ACQinfo.msPerLine*ACQinfo.linesPerFrame;
tdom =(0:(size(R{1},2)-1))*fr/1000; 
figure
for i = 1:length(idExamp)
     
    Right = squeeze(R{1}(idExamp(i),:));
    Up = squeeze(R{2}(idExamp(i),:));
    Left = squeeze(R{3}(idExamp(i),:));
    Down = squeeze(R{4}(idExamp(i),:));
    
    mu = median(Right);
    Right = (Right-mu)/mu;
    
    
    phidom = linspace(0,360,length(Right));
    
    subplot(length(idExamp),1,i)
%     plot(tdom,zscore(Right),'k'), hold on
%     plot(tdom,zscore(Up)+4,'k'), hold on
%     plot(tdom,zscore(Left)+8,'k'), hold on
%     plot(tdom,zscore(Down)+12,'k')
    
    plot(abs(fft(zscore(Right))),'k'), hold on
    plot(abs(fft(zscore(Up)+4)),'k'), hold on
    plot(abs(fft(zscore(Left)+8)),'k'), hold on
    plot(abs(fft(zscore(Down)+12)),'k')
    
    axis tight
    
end
%%


