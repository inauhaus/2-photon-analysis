function setsymbolstruct

%Put all the symbol information into global structure

global symbolInfo Analyzer G_handles

Nsym = length(Analyzer.loops.conds{1}.symbol);
symbolInfo = struct;
Fsymbol = get(G_handles.primSymbol,'string'); %primary parameter symbol in looper to analyze
symbolInfo.ID(1) = get(G_handles.primSymbol,'value');  %The index with respect to the looper
symbolInfo.str{1} = Fsymbol{symbolInfo.ID(1)};  %Selected string
symbolInfo.domType = get(G_handles.domType,'value');  %Type of domain for primary symbol... .e.g. circular 'Axis'

if Nsym > 1
    Fsymbol = get(G_handles.secSymbol,'string'); %secondary symbol
    symbolInfo.ID(2) = get(G_handles.secSymbol,'value');
    symbolInfo.str{2} = Fsymbol{symbolInfo.ID(2)};
    symbolInfo.Collapse(1) = get(G_handles.secCollapse,'value');  %Describes how to collapse across secondary loop domains
    
    if Nsym > 2
        
        Fsymbol = get(G_handles.tertSymbol,'string'); %tertiary symbol
        symbolInfo.ID(3) = get(G_handles.tertSymbol,'value');
        symbolInfo.str{3} = Fsymbol{symbolInfo.ID(3)};
        symbolInfo.Collapse(2) = get(G_handles.tertCollapse,'value');  %Describes how to collapse across tertiary loop domains
    end
    
end