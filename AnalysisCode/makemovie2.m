function makemovie2(trial,sig,twin_ms) 

%Directories and mask should be set from pF0


CHs = GetTrialData([1 0],trial);
 
  %%

global ACQinfo maskS cellS

fp = ACQinfo.linesPerFrame*ACQinfo.msPerLine;

xwin = 10:660; ywin = 10:500; 
%xwin = 1:size(CH,2); ywin = 1:size(CH,1); 
twin = round(twin_ms(1)/fp):round(twin_ms(2)/fp); %convert to frame window
twin_ms = twin*fp;
CH2 = CHs{1}(ywin,xwin,twin);

nID = getNeuronMask;  %get the index values for the neurons
masklabel = bwlabel(maskS.neuronmask(ywin,xwin),4);
celldom = unique(masklabel);

[xmicperpix ymicperpix] = getImResolution

%% Spatial filtering
if sig>0
    [x y] = meshgrid(1:length(xwin),1:length(ywin));
    x = x-mean(x(:)); y = y-mean(y(:)); r = sqrt(x.^2 + y.^2);
    hh = exp((-r.^2)/(2*sig^2));
    hh = hh/sum(hh(:));
    hh = abs(fft2(hh));
    
    for i = 1:length(twin)
        CH2(:,:,i) = ifft2(fft2(CH2(:,:,i)).*hh);
    end
end

mi = prctile(CH2(:),2);
ma = prctile(CH2(:),99.9);

%% get time courses from mask
dim = size(CH2(:,:,1));
Nt = length(CH2(1,1,:));
cellMat = zeros(length(nID),Nt);
for p = 1:length(nID)
   
    idlocs = find(masklabel == p);
    idlocs = idlocs(:)*ones(1,Nt) + ones(length(idlocs),1)*((0:Nt-1)*dim(1)*dim(2));
    
    cellMat(p,:) = mean(CH2(idlocs))';
    
end

clear spiketrain cellMat_norm
for i = 1:length(cellMat(:,1))
   
    trace = cellMat(i,:);
    id = find(trace>prctile(trace,0) & trace<prctile(trace,50));
    trace = (trace-median(trace(id)))/median(trace(id));   
    trace = phi(trace)
    %trace = trace/max(trace);
    
    %Deconvolve
    H = [trace(1:end-1)'];
    y = trace(2:end)';
    x = inv(H'*H)*H'*y;
    spiketrain(i,:) = y-H*x;
    spiketrain = phi(spiketrain);
    
    traceSeparation = 3;
    cellMat_norm(i,:) = trace+(i-1)*traceSeparation;
        
end


%%

global Analyzer
anim = Analyzer.M.anim;
expt = Analyzer.M.expt;
unit = Analyzer.M.unit;
fname = [anim '_' unit '_' expt '_trial' num2str(trial)];
mov = VideoWriter(['C:\Movies\'  fname '.avi','Grayscale AVI']);

set(mov,'FrameRate',15)


%%

 tdom = twin*fp/1000 - getparam('predelay');

%pID = nID([1:50]); %ab2 000_014

pID = nID

%pID = nID([35 54]); %ab2 000_014

twin_ms = twin_ms-twin_ms(1);

scaleBar = 200; %microns

fps = 1000/fp;
clear F
open(mov)
figure(20)
for f = 1:length(twin)
    
    
    tensdum = CH2(:,:,f);  
    
    
    figure(21)
    subplot(1,2,1)
    cla
    
    imagesc(tensdum,[mi ma]), colormap gray
    xlim manual, ylim manual
    hold on
    contour(maskS.neuronmask(ywin,xwin),[.5 .5],'r','LineWidth',2)
    hold on
    plot([0 scaleBar/xmicperpix],[size(tensdum,1)+10 size(tensdum,1)+10],'k','LineWidth',2) %plot scale bar
    text(0,size(tensdum,1)+30,[num2str(scaleBar) ' microns'])
    %hold on
    
%     
    axis image
%     xlim([xwin(1) xwin(end)])
%     ylim([ywin(1) ywin(end)])    
    axis off
    %drawnow
    
    
    subplot(1,2,2)
    xlim manual
    ylim manual
    if f == 1
        cla
        
        plot(twin_ms(1:end)/1000, cellMat_norm(:,1:length(twin))','k')
        
        hold on
        plot([twin_ms(1) twin_ms(1)]/1000 - 1 ,[0 1],'k','LineWidth',2)  %plot y-axis scale bar: 1 dF/F 
        text(-1.3,.5,'1.0 dF/F','HorizontalAlignment','right')
        
        hold on
        plot([0 5] ,[-2 -2],'k','LineWidth',2)  %plot x-axis scale bar: 5 sec
        text(0,-3,'5 sec','HorizontalAlignment','left')
        axis off
        
        
        xlabel('seconds')
        set(gca,'YTickLabels',[])
        
        
        xlim([twin_ms(1)-1500 twin_ms(end)]/1000)
        
        set(gca,'Color',[.94 .94 .94],'box','off')
    end
    
    hold on
    plot([twin_ms(1) twin_ms(f)]/1000,[-1 -1],'r','LineWidth',10) %plot the progress bar    
    ylim([-3 max(cellMat_norm(:))])
    
    %%%%%If I've loaded the image of the traces    
%     subplot(3,1,3)
%     hold on, plot([tdom(f) tdom(f)],[-4 -3],'-b','linewidth',2)
%     subplot(3,1,2)
%     hold on, plot([tdom(f) tdom(f)],[-4 -3],'-r','linewidth',2)
%     drawnow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    F(f) = getframe(gcf);
    %%mov = addframe(mov,F);
    
    %pause(1/fps)
end
writeVideo(mov,F)
close(mov)
%movie2avi(F,'RawTrial7.avi','fps',round(fps*2),'compression','Indeo3')

