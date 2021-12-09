function loadAnalyzer


global G_handles Analyzer

%Get analyzer file
anapath = get(G_handles.analyzedir,'String'); %partial path for analyzer file
anapath = [anapath '.analyzer'];
load(anapath,'Analyzer','-mat')
