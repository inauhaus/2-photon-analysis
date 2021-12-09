function [y y_var] = CondTensor5(b)

%Compute the tensor for each condition

%'4' checks for trials that don't match the template, and throws them out.

%The new version builds the entire Tensor, as opposed to truncating within a time interval
%
%
%b is a 2D vector corresponding the the beginning and end of
%the baseline subtraction images, in milliseconds. e.g. varargin = {[0 500]} sums
%the images from 0 to .5 seconds for each repetition and then subtracts it
%from the mean response in the repeat.
%
%slowMo performs movement correction
%Rflag fits a line to the red/green scatter plot and subtracts this trend
%from the data in the green channel

global ACQinfo bsflag repDom G_handles mbestall nbestall twophDATADIR maskS cellS



%%%Get the flags for image movement correction%%%%

slowXcorrFlag = get(G_handles.slowXcorr,'Value');  %lateral movement correction
fastXcorrFlag = get(G_handles.fastXcorr,'Value');  %lateral movement correction
OpticFlowFlag = get(G_handles.OpticFlowCorrection,'Value');  %lateral movement correction

frameRemovalFlag = get(G_handles.frameRemoval,'Value');  %lateral movement correction
deletionThresh_Z = str2num(get(G_handles.deletionThreshold,'string'));

SVDremovalFlag = get(G_handles.PCremoval,'Value');  %lateral movement correction
N_PC = str2num(get(G_handles.nPC,'string'));

opticFlowXsig = 2; %smooth in X before computing optic flow
rthresh_floor = .4; %absolute floor; correlation coefficient with the template.
%deletionThresh_Z = 1.5; %Threshold below the median in SD
%N_PC = 1; %number of principal components to remove

DSflag = get(G_handles.DSflag,'Value');

paramVec = [opticFlowXsig rthresh_floor deletionThresh_Z N_PC];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for c = 1:getnoconditions
    nReps(c) = getnorepeats(c);
end

%%%Make the mask for taking the baseline
b = b+getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning

if ACQinfo.SBInfo.scanbox_version==3
    trialstartID = find(ACQinfo.SBInfo.event_id == 1); %rising edge ID %RTO 8/17/21- changed from 3(rise) 2 (fall) for sbx 3.0
    trialendID = find(ACQinfo.SBInfo.event_id == 2); %falling edgeID
elseif ACQinfo.SBInfo.scanbox_version==2
    trialstartID = find(ACQinfo.SBInfo.event_id == 3); %rising edge ID %RTO 8/17/21- changed from 3(rise) 2 (fall) for sbx 3.0
    trialendID = find(ACQinfo.SBInfo.event_id == 2); %falling edgeID
else 
    disp('ERROR: Unknown Scanbox version')
end
FramesPerTrial = min(ACQinfo.SBInfo.frame(trialendID) - ACQinfo.SBInfo.frame(trialstartID));


dim = [ACQinfo.linesPerFrame ACQinfo.pixelsPerLine FramesPerTrial];
tdomdum = (0:prod(dim)-1) * (ACQinfo.msPerLine/ACQinfo.pixelsPerLine);
tdomdum = reshape(tdomdum,dim(2),dim(1),dim(3));
tdom = zeros(size(tdomdum,2),size(tdomdum,1),size(tdomdum,3));
for i = 1:size(tdom,3)
    tdom(:,:,i) = tdomdum(:,:,i)';
end
id = find(tdom > b(1) & tdom < b(2));
blim_mask = zeros(size(tdom));
blim_mask(id) = 1;

tdom1D = squeeze(tdom(1,1,:));
baseIDX = find(tdom1D > b(1) & tdom1D < b(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bimg_mu = zeros(dim(1),dim(2));
baseCounter = 0;

nt = getnotrials;
nc = getnoconditions;
y = cell(1,nc);
y2 = cell(1,nc);
ybase = cell(1,nc);

Px = 0;
Py = 0;

ExptMu = 0;

h = fspecial('gaussian',[length(ACQinfo.ydom) length(ACQinfo.xdom)],2);
h = abs(fft2(h));

mbestall = [];
nbestall = [];

varflag = 0;

k = 0;
for t = 1:nt
    
    percentdone = round(t/nt*100);
    set(G_handles.status,'string',[num2str(percentdone) '%']), drawnow
    
    [c r] = getcondrep(t);
    
    if stimblank(c);      %Identify blank condition since it will have a unique number of repeats
        repeatDomain = 1:getnorepeats(nc);
        %The following is to take only the blank repeats that are contained
        %within the selected repeats
        if repDom(end) < getnorepeats(c)
            bperr = getnorepeats(nc)/(getnorepeats(1));    %blanks per repeat
            repeatDomain = 1:floor(repDom(end)*bperr);
        end
    else
        repeatDomain = repDom;
    end
    
    if ~isempty(find(r == repeatDomain))
        
        k = k+1;
        
        if repeatDomain(1) == r
            y{c} = 0;
            y2{c} = 0;
        end
        
        T = str2num(G_handles.epistop.String);
        %CH_raw = GetTrialData([1 0],t,T);
        [CH_raw f1 fN] = GetTrialData([1 0],t);
        
        CH = CH_raw{1}(:,:,1:FramesPerTrial); %Truncate to the shortest trial
        fN = f1+FramesPerTrial-1;
        
        if r == 1
            TensCounter{c} = zeros(1,FramesPerTrial);
        end
        CHdim = size(CH);
        
        %% This takes a while        

        CH = TensorCleaner(CH,slowXcorrFlag,fastXcorrFlag,OpticFlowFlag,frameRemovalFlag,SVDremovalFlag,DSflag,[f1 fN],paramVec);
        idgoodFrames = find(~isnan(CH(1,1,:)));
        
        %%
%         figure(20)
%         for i = 1:size(CH,3)
%             cla
%             imagesc(CH(:,:,i))
%             pause(.4)
%             drawnow
%             
%         end
        %%
        %Baseline normalization within trial
        if bsflag == 1
            
            % bimg1 = nansum(CH.*blim_mask,3)./nansum(blim_mask,3); %mean
            bimg1 = nanmean(CH(:,:,baseIDX),3);
            %bimg1 = ifft2(fft2(bimg1).*h); %smooth it for stability
            
            if ~isnan(bimg1(1,1))
                
                
                for z = 1:length(CH(1,1,:))
                    CH(:,:,z) = CH(:,:,z) - bimg1;   %% baseline subtraction
                end
                
                bimg_mu = bimg_mu+bimg1;
                baseCounter = baseCounter+1;
                
                %             id = find(bimg1(:) == 0);
                %             bimg1(id) = NaN;
                %             for z = 1:length(CH(1,1,:))
                %                 CH(:,:,z) = CH(:,:,z)./bimg1;   %% baseline division
                %             end
                
            else
                
                CH = CH*NaN; %kill this trial if there is not a usable blank.
                idgoodFrames = [];
                'hi'
            end
            
        end
        
        TensCounter{c}(idgoodFrames) = TensCounter{c}(idgoodFrames) + 1;
        
        
        if r == 1
            y{c} = zeros(size(CH));
        end
        
        y{c}(:,:,idgoodFrames) = y{c}(:,:,idgoodFrames) + CH(:,:,idgoodFrames);
        
        
        if varflag
            y2{c}(:,:,idgoodFrames) = y2{c}(:,:,idgoodFrames) + (CH(:,:,idgoodFrames).^2);  %to compute variance later
        end
        
        t
        
    end
    
    clear CHs
    
    
end


%Baseline division using an average across all repeats and conditions
if bsflag
    bimg_mu = bimg_mu/baseCounter;
    
    for c = 1:length(y)
        
        for n = 1:size(y{c},3)
            y{c}(:,:,n) = y{c}(:,:,n)./bimg_mu;
        end
        
    end
    
end


for c = 1:length(y)

    for n = 1:size(y{c},3)

        y{c}(:,:,n) = y{c}(:,:,n)/TensCounter{c}(n);

        if varflag
            y2{c}(:,:,n) = y2{c}(:,:,n)/TensCounter{c}(n);
        end

    end

end


y_var = cell(1,length(y));

if varflag
    for i = 1:length(y)
       y_var{i} = y2{i} - y{i}.^2;
    end
end


% bflag = stimblank(getnoconditions); %if a blank exists in this experiment
% if ~bsflag && bflag  %if baseline subtraction not checked and blanks are in this experiment    
%     dum = mean(y{end},3);
%     f0blank = zeros(size(y{1}));
%     for z = 1:length(y{1}(1,1,:))
%         f0blank(:,:,z) = dum;
%     end
%     
%     for c = 1:length(y)
%         y{c} = (y{c} - f0blank)./f0blank;
%     end
%     
% end

