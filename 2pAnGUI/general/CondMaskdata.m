function CondMaskdata(slowMo,fastMo)

%Ian Nauhaus

%This is an adaptation of CondTensor.  It computes the cellS structure, but
%not the pF0 stuff... 'Tensor', 'f0m'

global ACQinfo repDom G_handles maskS cellS mbestall nbestall maskS

GcampFlag = 0;

alignCh_fast = 1;
alignCh_slow = 1;
if slowMo 
    
    temp = maskS.anatomyTemplate; %This was created during Set directories. It can be recreate via button in GUI
    
else
    temp = [];
end

if isfield(maskS,'maskTens')
    masklabel = zeros(size(maskS.maskTens(:,:,1)));
else
    masklabel = bwlabel(maskS.bwCell{1},4);
end

celldom = unique(masklabel);
Ncell = length(celldom);

nt = getnotrials;
nc = getnoconditions;

cellS = struct; %clear it
cellS.cellMat = cell(1,getnoconditions);

mbestall = [];
nbestall = [];
mbest = 0;
nbest = 0;

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
            CH_raw = GetTrialData([1 0],t,endT*1000);
        else  %Handle the last trial
            CH_raw = GetTrialData([1 0],t);        
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
               
                
        cleanflag = 0;
        tic
        if cleanflag
            CHdum = cleanTensor(CH,0); %This will downsample it.
            dim = size(CHdum);
            dimI = [ACQinfo.SBInfo.sz(1) length(ACQinfo.unblanked)];
            for i = 1:size(CH,3)
                imdum = CHdum(:,:,i);
                imdum = interp1(1:dim(1),imdum,linspace(1,dim(1),dimI(1)));
                imdum = interp1(1:dim(2),imdum',linspace(1,dim(2),dimI(2)))';
                CH(:,:,i) = imdum;
            end
        end
        toc
        
        %Apply slow movement correction
        if slowMo
            
            imdum = median(CH_raw{alignCh_slow}(:,:,2:end-2),3);
            
            
            [mbest nbest] = getShiftVals(imdum.^2,temp.^2,[mbest nbest])  %squaring seems to really help sometimes
            
            mbestall = [mbestall mbest];
            nbestall = [nbestall nbest];
            
            CH = circshift(CH,[round(-mbest) round(-nbest) 0]);
        end
        
        %Get waveforms from cell mask
        
        if fastMo
            
            %Get fast shift
            [Px_fast Py_fast] = getTrialMotion_tracker(CH_raw{alignCh_fast},0);
            
            %Get slow shift
            [Px_slow Py_slow] = getTrialMotion_tracker(CH_raw{alignCh_fast},1);
            
            %%%%Look for bad alignment%%%%%
            n1 = round(.05*length(Px_slow)); n2 = round(.95*length(Px_slow));
            if std(Px_slow(n1:n2))>2 | std(Py_slow(n1:n2))>2 | std(Px_fast(n1:n2))>2 | std(Py_fast(n1:n2))>2
                ['Warning:  Bad alignment on trial ' num2str(trialdom(tdum))]
            end
            %%%%%%%%%%%%%%%%%%%%%%%
            
            dum = getGeoTrxTimeCourse2(CH_raw{alignCh_fast},Px_slow+Px_fast+nbest,Py_slow+Py_fast+mbest);
            
        else
            dum = getCelltimecourse(CH,idcell,CoM,cellWidth);
        end
        
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
        
        %Get sync times from channel 3
%         if get(G_handles.loadSyncs,'value')
%             
%             Fs = 1000*ACQinfo.pixelsPerLine/ACQinfo.msPerLine;
%             
%             CH_raw{1} = []; %Helps keep it from running out of memory;
%             CH_raw{2} = [];
%             
%             [frameid lineid synctimes] = getsyncFrameTime(CH_raw{3});
%             cellS.synctimes{c,r} = synctimes;
%             cellS.frameid{c,r} = frameid;
%             
%             %Get Synctimes
%             dim = size(CH_raw{3});
%             for i = 1:length(CH_raw{3}(1,1,:))
%                 dum = CH_raw{3}(:,:,i)';
%                 CH_raw{3}(:,:,i) = reshape(dum(:),dim(1),dim(2));  %Don't want to create new variable, because of memory issues
%             end
%             
%             %             cellS.synctimes{c,r} = getFlashGraterSynctimes_CRT(CH_raw{3}(:),Fs);
%             %
%             %             id = find(cellS.synctimes{c,r} < getparam('predelay')/2);
%             %             cellS.synctimes{c,r}(id) = [];
%             %             id = find(cellS.synctimes{c} > getparam('predelay')+getparam('stim_time')+getparam('postdelay')/2);
%             %             cellS.synctimes{c,r}(id) = [];
%             
%             clear CH
%             
%         end
        
    end
    
end


