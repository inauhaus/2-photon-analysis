pRev

global G_RChandles G_handles

set(G_RChandles.kernelLength,'string','[-500 2000]');
set(G_RChandles.LPflag,'value',0);
set(G_RChandles.HPflag,'value',1);
set(G_RChandles.LPWind,'value',1);
set(G_RChandles.HPWind,'value',1);
set(G_RChandles.Lwidth,'string',50);
set(G_RChandles.Hwidth,'string',5000);
set(G_RChandles.blankNorm,'value',0);

%set(G_handles.datadir,'string','C:\2p_data\')
%set(G_handles.analyzedir,'string','C:\2p_data\AnalyzerFiles\')

set(G_handles.datadir,'string','e:\2p_data\')
set(G_handles.analyzedir,'string','e:\2p_data\AnalyzerFiles\')

%%  
getAllorisfinfo_color

getAllorisfinfo_DKLcolor

%%  Determine preferred axes of functional maps

sfOriaxes2(dori,dsf,ax,Dist,animID_d)

%%
pairWisePlots(abs(dori),abs(dsf),(doriEuc*2-1),(dsfEuc*2-1),abs(Dist))

%%
%CompareMapVals
orisfcolorHists
colorFormScatter
%PopSummary

%%
c = 4;
%BSclustering4(opref{c},sfpref{c},dphaseApair{c},tcoriall_fit{c},tcsfall_fit{c},tccolorall{c},abs(dori{c}),abs(dsf{c}),abs(dphaseA{c}),doriEuc{c},dsfEuc{c},dcolorEuc{c},Dist{c},animID,animID_d);  %Boot-strap clustering

BSclustering_color(lumpref{1},LMphaseDiff{1},abs(dlumpref{1}),abs(dLMphaseDiff{1}),Dist{1},animID,animID_d)
