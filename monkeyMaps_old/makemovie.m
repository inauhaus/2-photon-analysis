%%%%%%%%%
%To get this running first get the CoMx2 and CoMy2 from
%getGeoTrxTimeCourse2 by globalizing them within the file.  Then generate
%the figure of traces from revCorrPrediction.  Then execute these lines...
%%%%%%%%

CHs = GetTrialData([1 0],1);

%%

global ACQinfo maskS

fp = ACQinfo.linesPerFrame*ACQinfo.msPerLine;

xwin = 100:400; ywin = 100:400; 
%xwin = 1:size(CH,2); ywin = 1:size(CH,1); 
twin = round(10000/fp):round(51000/fp);
CH2 = CHs{1}(ywin,xwin,twin);

nID = getNeuronMask;  %get the index values for the neurons
masklabel = bwlabel(maskS.neuronmask,4);
celldom = unique(masklabel);


%%
[x y] = meshgrid(1:length(xwin),1:length(ywin));
x = x-mean(x(:)); y = y-mean(y(:)); r = sqrt(x.^2 + y.^2);
sig = .07;
hh = exp((-r.^2)/(2*sig^2));
hh = hh/sum(hh(:));
hh = abs(fft2(hh));

for i = 1:length(twin)
    CH2(:,:,i) = ifft2(fft2(CH2(:,:,i)).*hh);
end

mi = prctile(CH2(:),2);
ma = prctile(CH2(:),99.8)

%%

tdom = twin*fp/1000 - getparam('predelay');

%pID = nID([21 35 54]); %ab2 000_014

%pID = nID([35 54]); %ab2 000_014

fps = 1000/fp;

clear F
%mov = avifile('RawTrial.avi')
figure(20)
for f = 1:length(twin)/10
    tensdum = CH2(:,:,f);  
    
    
    figure(20)
    subplot(3,2,1)
    cla
    
    imagesc(tensdum,[mi ma]), colormap gray
    

    %hold on
    
%     id = find(CoMx2{twin(f)} >= xwin(1) & CoMy2{twin(f)} >= ywin(1) & CoMx2{twin(f)} <= xwin(end) & CoMy2{twin(f)} <= ywin(end));
%     plot(CoMx2{twin(f)}(id)-xwin(1)+1,CoMy2{twin(f)}(id)-ywin(1)+1,'.r','markersize',5)
    
%     hold on
%     plot(CoMx2{twin(f)}(pID(1))-xwin(1)+1,CoMy2{twin(f)}(pID(1))-ywin(1)+1,'ob','markersize',20)
%     hold on
%     plot(CoMx2{twin(f)}(pID(2))-xwin(1)+1,CoMy2{twin(f)}(pID(2))-ywin(1)+1,'or','markersize',20)
%     
    
    axis image
%     xlim([xwin(1) xwin(end)])
%     ylim([ywin(1) ywin(end)])    
    axis off
    %drawnow
    
    %%%%%If I've loaded the image of the traces    
    subplot(3,1,3)
    hold on, plot([tdom(f) tdom(f)],[-4 -3],'-b','linewidth',2)
    subplot(3,1,2)
    hold on, plot([tdom(f) tdom(f)],[-4 -3],'-r','linewidth',2)
    drawnow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %F(f) = getframe(gcf);
    %%mov = addframe(mov,F);
    
    pause(1/fps/8)
end
%movie2avi(F,'RawTrial7.avi','fps',round(fps*2),'compression','Indeo3')

