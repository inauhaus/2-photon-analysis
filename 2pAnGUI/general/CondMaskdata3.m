function CondMaskdata3

%Ian Nauhaus

%3 allows for downsampling

%This is an adaptation of CondTensor.  It computes the cellS structure, but
%not the pF0 stuff... 'Tensor', 'f0m'

global ACQinfo repDom G_handles maskS cellS


%%%Get the flags for image movement correction%%%%

slowMoFlag = get(G_handles.slowMotionFlag,'Value');  %lateral movement correction
fastMoFlag = get(G_handles.fastMotionFlag,'Value');  %lateral movement correction

frameRemovalFlag = get(G_handles.frameRemoval,'Value');  %lateral movement correction
deletionThresh_Z = str2num(get(G_handles.deletionThreshold,'string'));

SVDremovalFlag = get(G_handles.PCremoval,'Value');  %lateral movement correction
N_PC = str2num(get(G_handles.nPC,'string'));

opticFlowXsig = 2; %smooth in X before computing optic flow
rthresh_floor = .3; %absolute floor; correlation coefficient with the template.
%deletionThresh_Z = 1.5; %Threshold below the median in SD
%N_PC = 1; %number of principal components to remove

paramVec = [opticFlowXsig rthresh_floor deletionThresh_Z N_PC];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DSflag = 1;

if isfield(maskS,'maskTens')
    masklabel = zeros(size(maskS.maskTens(:,:,1)));
else
    masklabel = bwlabel(maskS.bwCell{1},4);
    
    if DSflag
       
        masklabel = Tensor_2x2bin(masklabel);  %2x2 binning
        masklabel = sign(masklabel);
        
    end
    
end

celldom = unique(masklabel);
Ncell = length(celldom);

nt = getnotrials;
nc = getnoconditions;

if isfield(cellS,'motionInfo')
    motionInfo = cellS.motionInfo;  
else
    motionInfo = [];
end

cellS = struct; %clear it....
cellS.motionInfo = motionInfo; %...but keep the motion info if that was loaded.
cellS.cellMat = cell(1,getnoconditions);


%masklabel = circshift(masklabel,[-50 50]);  %For the control

%Preallocate cell response matrix

cellS.cellMat = cell(1,getnoconditions);

%It saves time to create the following once and provide as an input to
%getcelltimecourse

for p = 1:Ncell
    
    [idcelly idcellx] = find(masklabel == celldom(p));
    
    CoM(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
    cellWidth(p,:) = [std(idcelly) std(idcellx)];  %std... cell width
    idcell{p} = find(masklabel(:) == celldom(p));
end


trialdom = 1:getnotrials;


k = 0;

for tdum = 1:length(trialdom)

    t = trialdom(tdum);

    percentdone = round(t/nt*100);
    set(G_handles.status,'string',[num2str(percentdone) '%']), drawnow
    
    [c r] = getcondrep(t);
    
    if stimblank(c);      %Identify blank condition since it will have a unique number of repeats
        repeatDomain = 1:getnorepeats(nc);
        %The following is to take only the blank repeats that are contained
        %within the selected repeats
        if repDom(end) < getnorepeats(1)
            bperr = getnorepeats(nc)/(getnorepeats(1));    %blanks per repeat
            repeatDomain = 1:floor(repDom(end)*bperr);
        end
    else
        repeatDomain = repDom;
    end
    
    if ~isempty(find(r == repeatDomain))
        
        k = k+1;
        
        if tdum < (length(trialdom)-2)
            endT = getparam('stim_time')+getparam('predelay')+getparam('postdelay')+2; %Grab more than trial length
            %endT = endT*2;
            [CH_raw f1 fN] = GetTrialData([1 0],t,endT*1000);
            %CH_raw = GetTrialData([1 0],t);
        else  %Handle the last trial
            [CH_raw f1 fN] = GetTrialData([1 0],t);
            Napp = round(1000*endT/(ACQinfo.msPerLine*ACQinfo.linesPerFrame)) - size(CH_raw{1},3);
            app = zeros([size(CH_raw{1}(:,:,1)) Napp]);
            dum = median(CH_raw{1},3);
            for x = 1:Napp
                app(:,:,x) = dum;
            end
            CH_raw{1}(:,:,end+1:end+Napp) = app;

        end
        
        t
        
        CH = CH_raw{1};
    
        fN = f1+size(CH,3)-1;
        
        CH = TensorCleaner(CH,slowMoFlag,fastMoFlag,frameRemovalFlag,SVDremovalFlag,DSflag,[f1 fN],paramVec);
        
        %Get waveforms from cell mask
        
        
        dum = getCelltimecourse(CH,idcell,CoM,cellWidth);
        
        %         figure,
        %         subplot(1,2,1), plot(cellS.motionInfo.rTemp{t})
        %         subplot(1,2,2), plot(dum')
        % drawnow
        
        %Preallocate
        if isempty(cellS.cellMat{c})
            cellS.cellMat{c} = zeros(length(dum(:,1)),length(dum(1,:)),getnorepeats(c));
        end        
        
        %Some trials are shorter longer.  This truncates to the shortest trial
        try
            cellS.cellMat{c}(:,:,r) = dum;
        catch
            newL = min([size(dum,2) size(cellS.cellMat{c}(:,:,r),2)]);
            dum = dum(:,1:newL);
            cellS.cellMat{c} = cellS.cellMat{c}(:,1:newL,:);
            
            cellS.cellMat{c}(:,:,r) = dum;
        end
        

    end
    
end


