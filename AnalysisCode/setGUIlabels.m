function setGUIlabels

global Analyzer G_handles

Gsetdirectories

conds = getnoconditions;
reps = getnorepeats(1);
set(G_handles.nocond,'string',num2str(conds))
set(G_handles.norep,'string',num2str(reps))
set(G_handles.dirstatus,'string','Loaded')

set(G_handles.predelay,'string',['predelay=' num2str(getParamVal('predelay'))])
set(G_handles.postdelay,'string',['postdelay=' num2str(getParamVal('postdelay'))])
set(G_handles.trialtime,'string',['trialtime=' num2str(getParamVal('stim_time'))])
    
Nsym = length(Analyzer.loops.conds{1}.symbol);

set(G_handles.primSymbol,'string',Analyzer.loops.conds{1}.symbol);  %Set the strings in drop down menu
set(G_handles.primSymbol,'value',1)

if Nsym > 1
    set(G_handles.secSymbol,'string',Analyzer.loops.conds{1}.symbol);  %Set the strings in drop down menu
    set(G_handles.secSymbol,'value',2)
    set(G_handles.secSymbol,'enable','on')
    set(G_handles.secCollapse,'enable','on')
    dom = getdomain(Analyzer.loops.conds{1}.symbol{2});
    domstr{1} = 'mean';
    for i = 2:length(dom)+1
        domstr{i} = dom(i-1);
    end
    set(G_handles.secCollapse,'string',domstr)
else
    set(G_handles.secSymbol,'enable','off')
    set(G_handles.tertSymbol,'enable','off')
    set(G_handles.secCollapse,'enable','off')
    set(G_handles.tertCollapse,'enable','off')
end

if Nsym > 2
    set(G_handles.tertSymbol,'string',Analyzer.loops.conds{1}.symbol);  %Set the strings in drop down menu
    set(G_handles.tertSymbol,'value',3)
    set(G_handles.tertSymbol,'enable','on')
    set(G_handles.tertCollapse,'enable','on')
    dom = getdomain(Analyzer.loops.conds{1}.symbol{3});
    for i = 1:length(dom)
        domstr{i} = dom(i);
    end
    set(G_handles.tertCollapse,'string',domstr)
else
    set(G_handles.tertSymbol,'enable','off')
    set(G_handles.tertCollapse,'enable','off')
end

