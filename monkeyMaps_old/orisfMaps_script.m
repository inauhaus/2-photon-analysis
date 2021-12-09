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

set(G_handles.datadir,'string','f:\2p_data\')
set(G_handles.analyzedir,'string','f:\2p_data\AnalyzerFiles\')

%%  
getAllorisfinfo2

%%  Determine preferred axes of functional maps

sfOriaxes2(dori,dsf,ax,Dist,animID_d)

%%
pairWisePlots(abs(dori),abs(dsf),-(doriEuc*2-1),-(dsfEuc*2-1),abs(Dist))

pairWisePlots_sfF1F0(dsf,dF1F0,abs(Dist))

%%
%CompareMapVals
%orisfHists
PopSummary

%%
id = find(abs(dori)>20 | abs(dsf)>1);
dphaseA(id) = NaN;
dphaseApair(id,:) = NaN;

%BSclustering2(opref,sfpref,dphaseApair,tcoriall_fit,tcsfall_fit,abs(dori),abs(dsf),abs(dphaseA),doriEuc,dsfEuc,Dist,animID,animID_d);  %Boot-strap clustering
BSclustering5(opref,sfpref,F1F0,dphaseApair,tcoriall_fit,tcsfall_fit,abs(dori),abs(dsf),abs(dF1F0),abs(dphaseA),doriEuc,dsfEuc,Dist,animID,animID_d);  %Boot-strap clustering
