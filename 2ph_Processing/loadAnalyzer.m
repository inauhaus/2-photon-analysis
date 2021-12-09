function loadAnalyzer(ue)

global anadir Analyzer G_handles
anadir
Anim = get(G_handles.loadana,'string');
fname = [Anim '_' ue];
path = [anadir fname '.analyzer']

load(path,'Analyzer','-mat')