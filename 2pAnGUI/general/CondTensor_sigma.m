function y = CondTensor_sigma(b,slowMo,fastMo,F0_1)

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

global ACQinfo bsflag repDom G_handles

alignCh = 1;
if slowMo
    chvec = [0 0 0 0];
    chvec(alignCh) = 1;
    CHtemp = GetTrialData(chvec,1); 
    if fastMo
        [Px_fast Py_fast] = getTrialMotion3(CHtemp{1});        
        CHtemp{1} = makeGeoTrx(CHtemp{1},Px_fast,Py_fast);
    end
    temp = mean(CHtemp{1}(:,:,2:end-1),3);
else
    temp = [];
end


framePer = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %frame period in ms
b = b+getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning

bframe1 = floor(b(1)/framePer) + 1;
bframe2 = ceil(b(2)/framePer) + 1;

nt = getnotrials;
nc = getnoconditions;
y = cell(1,nc);

HH = fspecial('gaussian', [ACQinfo.linesPerFrame ACQinfo.pixelsPerLine], .2);
HH = HH/sum(HH(:));
HH = abs(fft2(HH));

Px = 0;
Py = 0;

ExptMu = 0;

k = 0;
for t = 1:nt
    
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
        
        if repeatDomain(1) == r
            y{c} = 0;     
        end        
        
        CH_raw = GetTrialData([1 0 0 0],t);  
        t
        if fastMo
            [Px_fast Py_fast] = getTrialMotion3(CH_raw{1});
            
            CH = makeGeoTrx(CH_raw{1},Px_fast,Py_fast);
            
        else 
            CH = CH_raw{1};
        end
        
        %ImLast = CH;
        
        ExptMu = ExptMu+CH;
        
        %Apply movement correction
        if slowMo
            imdum = mean(CH(:,:,2:end-2),3);
            [mbest nbest] = getShiftVals(imdum,temp);  %get the transformation for this trial
            for z = 1:length(CH(1,1,:))
                CH(:,:,z) = circshift(CH(:,:,z),[-mbest -nbest]); %transform
            end
        end
        
        
        if bsflag == 1
            
            bimg1 = mean(CH(:,:,bframe1:bframe2),3);
            

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
        
        y{c} = y{c} + ((CH-F0_1{c+1}).^2)/length(repeatDomain);
        y{c} = y{c} + CH/length(repeatDomain);
        
        clear CHs
        
    end
    
end

bflag = stimblank(getnoconditions); %if a blank exists in this experiment
if ~bsflag && bflag  %if baseline subtraction not checked and blanks are in this experiment    
    dum = mean(y{end},3);
    f0blank = zeros(size(y{1}));
    for z = 1:length(y{1}(1,1,:))
        f0blank(:,:,z) = dum;
    end
    
    for c = 1:length(y)
        y{c} = (y{c} - f0blank)./f0blank;
    end
end



