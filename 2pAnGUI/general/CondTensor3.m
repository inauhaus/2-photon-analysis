function [y y_var] = CondTensor3(b,slowMo,fastMo)

%Compute the tensor for each condition

%The new version builds the entire Tensor, as opposed to truncating within a time interval 
%
%b is a 2D vector corresponding the the beginning and end of
%the baseline subtraction images, in milliseconds. e.g. varargin = {[0 500]} sums
%the images from 0 to .5 seconds for each repetition and then subtracts it
%from the mean response in the repeat.
%
%slowMo performs movement correction 
%Rflag fits a line to the red/green scatter plot and subtracts this trend
%from the data in the green channel

global ACQinfo bsflag repDom G_handles mbestall nbestall twophDATADIR

GcampFlag = 1;
tempTrial = 1;

alignCh = 1;
if slowMo || ZmovementDeletion
    chvec = [0 0 0 0];
    chvec(alignCh) = 1;
    
    frame1 = ACQinfo.SBInfo.frame((tempTrial*2)-1); %first frame of trial is every other sync
    frameN =  length(frame1:ACQinfo.SBInfo.frame((tempTrial*2)));
    dum = sbxread(twophDATADIR,frame1,frameN);
    
    dum = double(squeeze(dum(1,:,:,:)));
    
    dum = dum(:,ACQinfo.unblanked,:);
    
    temp = median(dum,3);
    
else
    temp = [];
end

%%%Make the mask for taking the baseline 
b = b+getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning

T = str2num(G_handles.epistop.String);
%CHdum = GetTrialData([1 0],1,T); 
CHdum = GetTrialData([1 0],1); 
FramesPerTrial = size(CHdum{1},3);

trialstartID = find(ACQinfo.SBInfo.event_id == 3); %rising edge ID
trialendID = find(ACQinfo.SBInfo.event_id == 2); %falling edgeID
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
mbest = 0;
nbest = 0;

varflag = 0;
mbest = 0;
nbest = 0;
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
        CH_raw = GetTrialData([1 0],t); 
        t
        if fastMo
            [Px_fast Py_fast] = getTrialMotion3(CH_raw{alignCh});
            
            CH = makeGeoTrx(CH_raw{1},Px_fast,Py_fast);
            
%             figure(4)
%             for j = 1:3
%                 for i = 1:length(CH_raw{1}(1,1,:))
%                     subplot(1,2,1)
%                     imagesc(CH(:,:,i)), colormap gray
%                     subplot(1,2,2)
%                     imagesc(CH_raw{1}(:,:,i)), colormap gray
%                     pause(.2)
%                 end
%             end
%                figure,imagesc(mean(CH(:,:,2:end-1),3)), colormap gray
%                figure,imagesc(mean(CHfilt(:,:,2:end-1),3)), colormap gray



%             if slowMo && k ~= 1               
%                 [Vx(k-1) Vy(k-1)] = getTrialShift(mean(CH(:,:,2:end-1),3),mean(ImLast(:,:,2:end-1),3));
%                 Px(k) = Px(k-1) + Vx(k-1);
%                 Py(k) = Py(k-1) + Vy(k-1); %update the position
%                 
%                 
%                 figure(23),
%                 plot([Py(:) Px(:)])
%                 
%                 CH = makeGeoTrx(CH,Px(k)*ones(1,length(CH(:))),Py(k)*ones(1,length(CH(:))));
%                 
% 
%                 %figure(22),plot([Vx(:) Vy(:)])
%             end 
        else 
            CH = CH_raw{1};
        end
        
        CH = CH(:,:,1:FramesPerTrial); %Truncate to the shortest trial
        
        %ImLast = CH;
        
        ExptMu = ExptMu+CH;
        
        %Apply movement correction
        
        if slowMo
    
            if GcampFlag == 1
                imdum = median(CH_raw{alignCh}(:,:,2:end-2),3);     
            else
                imdum = mean(CH_raw{alignCh}(:,:,2:end-2),3); 
            end
            
            [mbest nbest] = getShiftVals(imdum.^2,temp.^2,[mbest nbest])  %squaring seems to really help sometimes
            
            mbestall = [mbestall mbest];
            nbestall = [nbestall nbest];

            CH = circshift(CH,[-round(mbest) -round(nbest) 0]);

        end
        
        %Baseline normalization
        if bsflag == 1
                 
            bimg1 = sum(CH.*blim_mask,3)./sum(blim_mask,3);
            bimg1 = ifft2(fft2(bimg1).*h);
            
            for z = 1:length(CH(1,1,:))
                CH(:,:,z) = CH(:,:,z) - bimg1;   %% baseline subtraction
            end
            
            id = find(bimg1(:) == 0);
            bimg1(id) = NaN;
            for z = 1:length(CH(1,1,:))
                CH(:,:,z) = CH(:,:,z)./bimg1;   %% baseline division
            end
            
        end
        
        
        %Average repeats:  /nr is important for when blanks have different
        %number of reps
        
%         for k = 1:length(CH(1,1,:))
%             hh = makeMapFilter;
%             CH(:,:,k) = ifft2(fft2(CH(:,:,k)).*abs(fft2(hh)));
%         end
        
        y{c} = y{c} + CH/length(repeatDomain);
        
        if varflag
            y2{c} = y2{c} + (CH.^2)/length(repeatDomain);  %to compute variance later
        end
               
  
        clear CHs
        
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

