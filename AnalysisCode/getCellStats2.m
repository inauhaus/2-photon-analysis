function getCellStats2

global ACQinfo cellS bsflag G_handles Analyzer

bsflag = get(G_handles.basesub,'value');  %Frame start in ms (to average)

Flim = str2double(get(G_handles.epistart,'String'));  %Frame start in ms (to average)
Flim(2) = str2double(get(G_handles.epistop,'String')); %Frame stop in ms (to average)

framePer = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %frame period in ms
Flim = Flim+getParamVal('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning

Frame1 = floor(Flim(1)/framePer) + 1;
Frame2 = ceil(Flim(2)/framePer) + 1;

Nt = length(cellS.cellMat{1}(1,:,1));%N time points
Ncell = length(cellS.cellMat{1}(:,1,1)); 
dom = linspace(-framePer*Nt/2,framePer*Nt/2,Nt);
sig = 300;  %ms
smoother = exp(-dom.^2/(2*sig^2));
smoother = smoother/sum(smoother);
smoother = ones(Ncell,1)*smoother;
smoother = abs(fft(smoother,[],2));

nf = length(Frame1:Frame2);

%Baseline normalization
if bsflag
    blim = str2double(get(G_handles.bstart,'String')); %in msec as well
    blim(2) = str2double(get(G_handles.bstop,'String'));
    blim = blim+getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning
    bframe1 = floor(blim(1)/framePer) + 1;
    bframe2 = ceil(blim(2)/framePer) + 1;
    for i = 1:length(cellS.cellMat)  %loop through each condition
        cellS.cellMat_norm{i} = zeros(size(cellS.cellMat{i})); %Preallocate
        for j = 1:length(cellS.cellMat{i}(:,1,1)) %loop through each cell
            condcellMat = squeeze(cellS.cellMat{i}(j,:,:)); %[time repeats]
            dim = size(condcellMat);
%             if dim(1)<dim(2) %to account for conditions with 1 repeat
%                 condcellMat = condcellMat';
%             end            
            blank = mean(condcellMat(bframe1:bframe2,:));
            blank = ones(length(condcellMat(:,1)),1)*blank;
            condcellMat = (condcellMat-blank)./blank;
            cellS.cellMat_norm{i}(j,:,:) = condcellMat;
        end
    end
    
end

cdom = [];
for i = 1:length(cellS.cellMat) 
    if ~isempty(cellS.cellMat{i})
        cdom = [cdom i];
    end
end


%%


%Stats
for q = 1:length(cdom)  %loop through each condition
    i = cdom(q);

    nr = getnorepeats(i);
    
    dim = size(cellS.cellMat{i});
    
    if bsflag
        cellS.muTime{i} = squeeze(mean(cellS.cellMat_norm{i},3)); %mean across repeats
        cellS.muRep{i} = squeeze(mean(cellS.cellMat_norm{i}(:,Frame1:Frame2,:),2)); %mean across time
    else
        cellS.muTime{i} = squeeze(mean(cellS.cellMat{i},3)); %mean across repeats
        cellS.muRep{i} = squeeze(mean(cellS.cellMat{i}(:,Frame1:Frame2,:),2)); %mean across time
    end
    
    %remove artifact
    
    %cellS.muTime{i} = cellS.muTime{i} - DCArtifact;

    
    %%%%%%%%%%%%%%%
    
    
    cellS.mu{i} = mean(cellS.muTime{i}(:,Frame1:Frame2),2);  %mean across repeats and time window
    
    if strcmp(Analyzer.P.type,'PG') %requires t_period to compute F1
        for p = 1:length(cellS.muTime{i}(:,1))
            cellS.F1{i}(p) = getF1(cellS.muTime{i}(p,:));
        end
    end
    
    
    %%%%%%%%%%%
    Nt = length(cellS.cellMat{i}(1,:,1));
    Ncell = length(cellS.cellMat{1}(:,1,1));
    dom = linspace(-framePer*Nt/2,framePer*Nt/2,Nt);
    sig = 300;  %ms
    smoother = exp(-dom.^2/(2*sig^2));
    smoother = smoother/sum(smoother);
    smoother = ones(Ncell,1)*smoother;
    smoother = abs(fft(smoother,[],2));
    %%%%%%%%%%%%%%%%
    
    
    muTimesmooth = ifft(fft(cellS.muTime{i},[],2).*smoother,[],2);  %smooth before taking max
    cellS.maxi{i} = max(muTimesmooth(:,Frame1:Frame2),[],2);  %mean across repeats; max over time window

%     dum = reshape(cellS.cellMat{i}(:,Frame1:Frame2,:),dim(1),dim(3)*length(Frame1:Frame2));
%     cellS.sig{i} = std(dum,[],2)/sqrt(nr*nf); %std error across time and repeats
    
    %I think computing the mean across time first (as done below) is a
    %better way to compute the standard error for each condition.
    %Otherwise, the standard error becomes really low due to the number of
    %time samples (as above). It also doesn't really make sense to compute
    %a standard deviation across two dimensions (i.e. repeats and time
    %points)
    if bsflag
        cellS.sigTime{i} = std(cellS.cellMat_norm{i},[],3)/sqrt(nr); %std error across repeats
        dum = squeeze(mean(cellS.cellMat_norm{i}(:,Frame1:Frame2,:),2));
    else
        cellS.sigTime{i} = std(cellS.cellMat{i},[],3)/sqrt(nr); %std error across repeats
        dum = squeeze(mean(cellS.cellMat{i}(:,Frame1:Frame2,:),2));
    end
    cellS.sig{i} = std(dum,[],2)/sqrt(nr); %std error across repeats, for each condition
    
end


%% Make an N-dim kernel for each cell.  Each dimension is a looping variable

cellS.F1kern = cell(1,size(cellS.mu{1},1));  %each element is a neuron.
cellS.mukern = cell(1,size(cellS.mu{1},1));  %each element is a neuron.
cellS.sigkern = cell(1,size(cellS.mu{1},1));  %each element is a neuron.
cellS.mukernTime = cell(1,size(cellS.mu{1},1));
cellS.sigkernTime = cell(1,size(cellS.mu{1},1));
cellS.Repkern = cell(1,size(cellS.mu{1},1));

nr = getnorepeats(1);
blank = strcmp('blank',Analyzer.loops.conds{end}.symbol{1});
for c = 1:(getnoconditions-blank)  %loop each condition
    for p = 1:length(cellS.mu{c}) %loop each cell
        if strcmp(Analyzer.P.type,'PG') %requires t_period to compute F1
            cellS.F1kern{p}(c) = cellS.F1{c}(p);
        end
        cellS.mukern{p}(c) = cellS.mu{c}(p);
        cellS.mukernTime{p}(c,:) = cellS.muTime{c}(p,:);
        cellS.sigkern{p}(c) = cellS.sig{c}(p);
        cellS.sigkernTime{p}(c,:) = cellS.sigTime{c}(p,:);
        for r = 1:nr
            cellS.Repkern{p}(c,r) = cellS.muRep{c}(p,r); %This will be the response to each condition and repeat, averaged over the time window.
        end
    end
end

for i = 1:length(Analyzer.L.param)  %loop each looping parameter
    kdim(i) = length(eval(Analyzer.L.param{i}{2}));
end

if size(kdim) == 1
    kdim = [kdim 1];
end

for p = 1:length(cellS.mukern)
    if strcmp(Analyzer.P.type,'PG') %requires t_period to compute F1
        cellS.F1kern{p} = reshape(cellS.F1kern{p}(:),kdim);
    end
    cellS.mukern{p} = reshape(cellS.mukern{p}(:),kdim);
    cellS.mukernTime{p} = reshape(cellS.mukernTime{p}(:),[kdim size(cellS.mukernTime{p},2)]) ;
    
    cellS.sigkern{p} = reshape(cellS.sigkern{p}(:),kdim);
    cellS.sigkernTime{p} = reshape(cellS.sigkernTime{p}(:),[kdim size(cellS.sigkernTime{p},2)]) ;
    
    %This will be the response to each condition and repeat, averaged over the time window. 
    %repeats will be last dimension.  Averaging over repeats will give the
    %same thing as 'mukern'
    cellS.Repkern{p} = reshape(cellS.Repkern{p}(:),[kdim nr]); 
    
    
end
    
if blank
    for p = 1:length(cellS.mukern)
        cellS.mublank{p} = cellS.mukern{p}(end);
    end
end

