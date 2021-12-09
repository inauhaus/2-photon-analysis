function [sfpref_raw ocdommap_raw dprime2 magodgrad dphaseROI dphaseROISEL] = CombineOcdomSFexpt(anim,dpthresh)

%This has to be used instead of SfOcdom.m if the SF and OD experiment were
%not part of the same looper.  i.e. it loads two experiments

global f0m f0m_var funcmap symbolInfo Analyzer G_handles flipeyebit

%% Load ORI OD expt

switch  anim
    case 'ab8'
        expt = 'u001_017';   
    case 'ab9'
        expt = 'u000_102';
    case 'ac0'        
        expt = 'u000_056';      
end

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels
load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

[xmicperpix ymicperpix] = getImResolution(1);
set(G_handles.Lwidth,'string',num2str(1/xmicperpix*3.823));
hh = makeMapFilter;
plotMapExamples_oriocdom(hh) %shows ori and od curves

%%


set(G_handles.HPflag,'Value',0);
set(G_handles.LPflag,'Value',1);

set(G_handles.Lwidth,'string',num2str(9/xmicperpix*3.823))
hh = makeMapFilter;
set(G_handles.Lwidth,'string',num2str(.8/xmicperpix*3.823));
hh2 = makeMapFilter;
anatomyflag = 0;
bwCellPlot = ones(size(funcmap));



%% get ocdom map

for i = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp(Analyzer.loops.conds{1}.symbol{i},'Leye_bit');
        idsym = i;
        break
    end
end

symbolInfo.ID(1) = idsym;
set(G_handles.primSymbol,'value',idsym); 
if idsym == 1
    set(G_handles.secSymbol,'value',2);
else
    set(G_handles.secSymbol,'value',1);
end
    
setsymbolstruct

ocdommap = GprocessBinary2(f0m,hh,hh);
[ocdommap_raw odmask] = GprocessBinary2(f0m,hh2,hh);

%mid = (prctile(ocdommap(find(dprimemask)),90) + prctile(ocdommap(find(dprimemask)),10))/2;
mid = nanmean(ocdommap(:));
ocdommap = ocdommap-mid;
ocdommap_raw = ocdommap_raw-mid;

%This is to make it so that positive is for the contralateral eye
if flipeyebit
    ocdommap = -ocdommap;
    ocdommap_raw = -ocdommap_raw;
end

%% Load ORI SF expt

switch  anim
    case 'ab8'
        expt = 'u001_007';        
    case 'ab9'
        expt = 'u000_084';        
    case 'ac0'        
        expt = 'u000_047';  %ori/sf expt
end

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])


%% get sfmap

sfdom = getdomain('s_freq');
for i = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp(Analyzer.loops.conds{1}.symbol{i},'s_freq');
        idsym = i;
        break
    end
end

symbolInfo.ID(1) = idsym;
set(G_handles.primSymbol,'value',idsym); 
if idsym == 1
    set(G_handles.secSymbol,'value',2);
else
    set(G_handles.secSymbol,'value',1);
end
    
setsymbolstruct

funcmap = GprocessLog(f0m,bwCellPlot,hh);   %output is complex
sfmag = real(funcmap);
sfpref = imag(funcmap);

funcmap = GprocessLog(f0m,bwCellPlot,hh2);   %output is complex
sfmag_raw = real(funcmap);
sfpref_raw = imag(funcmap);



edgeTrunc = 8;
dim = size(ocdommap);
trY = (edgeTrunc+1):(dim(1)-edgeTrunc);
trX = (edgeTrunc+1):(dim(2)-edgeTrunc);

%% Get mask for data selection

[dprime2 dprimemask] = getMapMask(dpthresh,hh2,hh,edgeTrunc);
dprimemask = dprimemask.*round(odmask);
dprime2 = dprime2.*odmask;

%% Plot stuff
[magodgrad dphaseROI dphaseROISEL] = plotSFODstuff(sfpref,sfpref_raw,ocdommap,ocdommap_raw,dprime2,dprimemask,trX,trY);