function experimentMotionAnalyzer2

global maskS cellS ACQinfo twophDATADIR

tsig = 2;  %temporal smoothing sigma, in "frames"

temp = maskS.anatomyTemplate;  %This could be from a different experiment if desired when setting directory

cellS.motionInfo = struct;

Nframes = ACQinfo.SBInfo.frame(end);

Win = 30; %sec  %load blocks in time that are this long

acqPeriod = ACQinfo.msPerLine*ACQinfo.linesPerFrame; %ms
framesPerWin = round(Win*1000/acqPeriod);
fstarts = 1:framesPerWin:Nframes;
fstarts = fstarts(1:end-1);  %Last window will be longer.
fends = [fstarts(2:end)-1 Nframes];
Nfs = fends-fstarts+1; %number of frames in each window

mbest = 0;
nbest = 0;

cellS.motionInfo.rTempFrames = [];

mbestin = 0;
nbestin = 0;

lpsig = .5;
hpsig = inf;

trunc = 0;
thresh = .2; %This should be really low. Its only for severe fuck ups

temp = temp(trunc+1:end-trunc,trunc+1:end-trunc);

k = 1;
for Wid = 1:length(fstarts) %loop each window in the experiment
    Wid
    imgs_uncorrected = single(sbxread(twophDATADIR,fstarts(Wid),Nfs(Wid)));

    imgs_uncorrected = squeeze(imgs_uncorrected(1,:,:,:));
    imgs_uncorrected = imgs_uncorrected(:,ACQinfo.unblanked,:);
    
    %imgs_uncorrected = smoothTensorTime(imgs_uncorrected,tsig); %smooth in time;
    imgs_uncorrected = TensorFilter_space(imgs_uncorrected,lpsig,hpsig); %filter in space
    
    imgs_uncorrected = imgs_uncorrected(trunc+1:end-trunc,trunc+1:end-trunc,:);
    
    imgs = imgs_uncorrected; %preallocate
    
    clear rdum
     for n = 1:size(imgs,3)
        imdum = imgs_uncorrected(:,:,n);
        [mbest nbest] = getShiftVals(imdum.^2,temp.^2,[mbestin nbestin]);  %squaring seems to really help sometimes
        %imgshift = circshift(imgs(:,:,n),[round(-mbest) round(-nbest) 0]);  %need to shift to find the "bad frames"
        imgshift = circshift_continous2(imgs(:,:,n),-nbest,-mbest);
        
        [r p] = corrcoef(temp(:),imdum(:));
        rdum(n) = r(1,2);
        pdum(n) = p(1,2);
        if rdum(n)<thresh  %if its a terrible match after shift, don't shift
            imgs(:,:,n) = imgs_uncorrected(:,:,n);
            mbest = mbestin;
            nbest = nbestin;         
        else
            imgs(:,:,n) = imgshift;
            mbestin = mbest;
            nbestin = nbest;
        end
        
        cellS.motionInfo.mbestFrame(k) = mbest;
        cellS.motionInfo.nbestFrame(k) = nbest;
        
        imdum = imgs(:,:,n);
        
        mbestdum(n) = mbest;
        nbestdum(n) = nbest;
        


            
        cellS.motionInfo.rTempFrames(k) = rdum(n);
        
        k = k+1;
            
     end

    %figure,plot(rdum)
    
    %%
%     mbestdum = medfilt1((mbestdum),5);
%     nbestdum = medfilt1((nbestdum),5);
%     mbestdum = round(mbestdum);
%     nbestdum = round(nbestdum);
         [xgrid ygrid] = meshgrid(20:50:600,20:50:600);
    %%
%      mi = prctile(imgs(:),1);
%     ma = prctile(imgs(:),99.5);
    
    %%
%     figure(95)
%     for n = 1:size(imgs,3)
%         clf
%         subplot(1,2,1)
%         imagesc(imgs(:,:,n),[mi ma]), 
%         %imagesc(ToI(:,:,n))
%         hold on, plot(xgrid,ygrid,'.r')
%         title(['corrected ' num2str(mbestdum(n))])
%         
%         %         if rdum(n)<thresh
%         %             plot(xgrid,ygrid,'xr')
%         %             drawnow
%         %         end
%         subplot(1,2,2)
%         imagesc(imgs_uncorrected(:,:,n),[mi ma]), colormap gray
%         hold on, plot(xgrid+nbestdum(n),ygrid+mbestdum(n),'.r')
%         drawnow
%         
%         
%         pause(.001)
%         
%         n
%     end
    
    
    %%
    
    
    
    
     
end


%cellS.motionInfo.mbestFrame = medfilt1(cellS.motionInfo.mbestFrame,5);
%cellS.motionInfo.nbestFrame = medfilt1(cellS.motionInfo.nbestFrame,5);

%%



%% Organize by trials as well


trialstartID = find(ACQinfo.SBInfo.event_id == 3); %Find events corresponding to rising edge from Stimulus
firstFrames = ACQinfo.SBInfo.frame(trialstartID);
trialendID = find(ACQinfo.SBInfo.event_id == 2); %Find events corresponding to falling edge from Stimulus
lastFrames = ACQinfo.SBInfo.frame(trialendID);

for t = 1:getnotrials
    
    f1 = firstFrames(t);
    f2 = lastFrames(t);
    
    cellS.motionInfo.mbest{t} = cellS.motionInfo.mbestFrame(f1:f2);
    cellS.motionInfo.nbest{t} = cellS.motionInfo.nbestFrame(f1:f2);
    cellS.motionInfo.rTemp{t} = cellS.motionInfo.rTempFrames(f1:f2);
    
end

    



