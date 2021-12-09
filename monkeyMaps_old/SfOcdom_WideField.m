function [sfpref_raw ocdommap_raw dprimemask magodgrad dphaseROI dphaseROISEL] = SfOcdom_WideField(dpthresh)

global f0m f0m_var funcmap symbolInfo Analyzer G_handles flipeyebit

[xmicperpix ymicperpix] = getImResolution(1);

set(G_handles.HPflag,'Value',0);
set(G_handles.LPflag,'Value',1);

set(G_handles.Lwidth,'string',num2str(9/xmicperpix*3.823));
hh = makeMapFilter;
set(G_handles.Lwidth,'string',num2str(.8/xmicperpix*3.823));
hh2 = makeMapFilter;

bwCellPlot = ones(size(funcmap));


%% get sfmap

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

%Compare sfmap for each eye
% funcmap = GprocessLog(f0m,bwCellPlot,hh2,0);   %output is complex
% sfmagA = real(funcmap);
% sfprefA = imag(funcmap);
% funcmap = GprocessLog(f0m,bwCellPlot,hh2,1);   %output is complex
% sfmagB = real(funcmap);
% sfprefB = imag(funcmap);
% figure,scatter(sfprefB(10:10:end),sfprefA(10:10:end),'.'), hold on, plot([0 6],[0 6],'r')
% odmap = sfmagB-sfmagA;



% sfpref_raw = medfilt2(sfpref_raw,[3 3]);
% sfpref = medfilt2(sfpref,[3 3]);
% 
% id = find(sfpref_raw<0.5);
% sfpref_raw(id) = 0.5;
% id = find(sfpref_raw>8);
% sfpref_raw(id) = 8;
% 
% id = find(sfpref<0.5);
% sfpref(id) = 0.5;
% id = find(sfpref>8);
% sfpref(id) = 8;



%% get ocdom map (set symbol order to eye/sfreq/misc)


for i = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp(Analyzer.loops.conds{1}.symbol{i},'Leye_bit');
        idsym = i;
        break
    end
end

for i = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp(Analyzer.loops.conds{1}.symbol{i},'s_freq');
        idsym2 = i;
        break
    end
end

symbolInfo.ID(1) = idsym;
set(G_handles.primSymbol,'value',idsym); 

symbolInfo.ID(2) = idsym2;
set(G_handles.secSymbol,'value',idsym2); 

if idsym + idsym2 == 3 
    set(G_handles.tertSymbol,'value',3);
elseif idsym + idsym2 == 4 
    set(G_handles.tertSymbol,'value',2);
elseif idsym + idsym2 == 5
    set(G_handles.tertSymbol,'value',1);
end

setsymbolstruct

[ocdommap] = GprocessBinary2(f0m,hh,hh);
[ocdommap_raw odmask] = GprocessBinary2(f0m,hh2,hh);

%mid = (prctile(ocdommap(find(dprimemask)),90) + prctile(ocdommap(find(dprimemask)),10))/2;
mid = nanmean(ocdommap(:));
ocdommap_raw = ocdommap_raw - mid;
ocdommap = ocdommap - mid;

%This is to make it so that positive is for the contralateral eye
if flipeyebit
    ocdommap = -ocdommap;
    ocdommap_raw = -ocdommap_raw;
end

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


