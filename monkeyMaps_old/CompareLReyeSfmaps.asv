global bw f0m f0m_var funcmap ACQinfo symbolInfo Analyzer G_handles idExamp kernelsIm 

set(G_handles.datadir,'string','e:\2p_data\')
set(G_handles.analyzedir,'string','e:\2p_data\AnalyzerFiles\')

%set(G_handles.datadir,'string','C:\2p_data\')
%set(G_handles.analyzedir,'string','C:\2p_data\AnalyzerFiles\')

%%%%%%%First, left hemisphere%%%%%%%%%%%

idExamp = [];  %left hemi

%make sure to hit "Set Directory" in GUI

%% ac1; u000_115

anim = 'ac1';
expt = 'u000_115';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 3;
Sf_eachEye_WideField(dthresh)

%% ac1; u001_011 

anim = 'ac1';
expt = 'u001_011';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 2;
Sf_eachEye_WideField(dthresh)


%% ac1; u002_019 

anim = 'ac1';
expt = 'u002_019';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 3;
Sf_eachEye_WideField(dthresh)