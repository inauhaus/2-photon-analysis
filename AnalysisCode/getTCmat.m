function [kernPop kernSigPop blankMu blankSig] = getTCmat

%Just get the f0 kernel for each neuron in the mask

global cellS maskS Analyzer

getCellStats  %reset the mean based on time window

masklabel1 = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel1);
nID = celldom(2:end);


Nsym = length(Analyzer.loops.conds{1}.symbol);  %number of looping params
for idsym = 1:Nsym
    if strcmp('ori',Analyzer.loops.conds{1}.symbol{idsym})
        oriid = idsym;
    elseif strcmp('greengain',Analyzer.loops.conds{1}.symbol{idsym})
        greenid = idsym;
    elseif strcmp('bluegain',Analyzer.loops.conds{1}.symbol{idsym})
        blueid = idsym;
    else
        noc = idsym;  %to identify a third condition, such as sfreq
    end
end

Nsym = length(Analyzer.loops.conds{1}.symbol);  %number of looping parameters

for i = 1:Nsym
    allDom{i} = getdomain(Analyzer.loops.conds{1}.symbol{i});
end
if Nsym == 1
    allDom{2} = NaN;
end

nc = getnoconditions;
bflag = stimblank(getnoconditions); %if a blank exists in this experiment
if bflag
    nc = nc-1;
end

for i = 1:length(nID)   %loop through each neuron
    
    kern = zeros(length(allDom{1}),length(allDom{2}),length(allDom{3}));
    kernSig = zeros(length(allDom{1}),length(allDom{2}),length(allDom{3}));

    for c = 1:nc
        for s = 1:Nsym       
            val = Analyzer.loops.conds{c}.val{s};
            idsym(s) = find(val == allDom{s});            
        end
        if Nsym == 1
            idsym(2) = 1;
        end
        
        if Nsym == 2
            kern(idsym(1),idsym(2)) = cellS.mu{c}(nID(i),:);
            kernSig(idsym(1),idsym(2)) = cellS.sig{c}(nID(i),:);
        elseif Nsym == 3
            kern(idsym(1),idsym(2),idsym(3)) = cellS.mu{c}(nID(i),:);
            kernSig(idsym(1),idsym(2),idsym(3)) = cellS.sig{c}(nID(i),:);
        end
    end
    
    if Nsym == 4
        kern = squeeze(mean(kern,noc));
        kernSig = squeeze(mean(kernSig,noc))/sqrt(size(kern,noc));
    end

%     if oriid == 1
%         kern = kern';
%         kernSig = kernSig';
%         dim = size(kern);
%     end
    
    if bflag
        blankMu(i) = cellS.mu{end}(nID(i));
        blankSig(i) = cellS.sig{end}(nID(i));
    else
        blankMu(i) = NaN;
        blankSig(i) = NaN;
    end

    for c = 1:dim(1)  %each contrast

        kernPop(c,:,i) = kern(c,:);
        kernSigPop(c,:,i) = kernSig(c,:);

    end
            
end
