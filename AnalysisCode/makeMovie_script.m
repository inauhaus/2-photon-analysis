pF0

%initialize the Gui

global G_handles Analyzer cellS maskS idExamp ACQinfo


anim = 'fry';
expt = 'u005_004';
dataRoot = 'f:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end) ''])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

global twophDATADIR ACQinfo

%%%%%%%%%
%To get this running first get the CoMx2 and CoMy2 from
%getGeoTrxTimeCourse2 by globalizing them within the file.  Then generate
%the figure of traces from revCorrPrediction.  Then execute these lines...

%Also, for the movie in the NN paper, I was using ab2_000_014, trial 13,
%cells 35 and 54.
%%%%%%%%

%CHs = GetTrialData([1 0 0 0],3);

%% Get the raw data

%imgs = single(sbxread(twophDATADIR,1,2000));
fps = 1000/(ACQinfo.msPerLine*ACQinfo.linesPerFrame);
Nsec  = 30; 
Nframes = round(30*fps);
imgs = single(sbxread([dataRoot anim '\' anim '_' expt(2:end) ''],1,Nframes));
CHs{1} = squeeze(imgs);
CHs{1} = CHs{1}(:,ACQinfo.unblanked,:);



%% Redefine

dims = size(CHs{1});
xwin = 1:dims(2); ywin = 1:dims(1); twin = 1:dims(3);
CH2 = CHs{1}(ywin,xwin,twin);

%% Smooth each frame

[x y] = meshgrid(1:length(xwin),1:length(ywin));
x = x-mean(x(:)); y = y-mean(y(:)); r = sqrt(x.^2 + y.^2);
sig = .3;
hh = exp((-r.^2)/(2*sig^2));
hh = hh/sum(hh(:));
hh = abs(fft2(hh));

for i = 1:length(twin)
    CH2(:,:,i) = ifft2(fft2(CH2(:,:,i)).*hh);
end
%% Template image
tempIm = median(CH2,3);

%%
mi = prctile(tempIm(:),1);
ma = prctile(tempIm(:),99.999);

%mi = prctile(CH2(:),5);
%ma = prctile(CH2(:),99.9);


tdom = twin/fps - getparam('predelay');

%pID = nID([21 35 54]); %ab2 000_014

%pID = nID([35 54]); %ab2 000_014


clear F

figure(20)
%%
clear F

for f = 1:length(twin)
    
    tensdum = CH2(:,:,f);  
    
    tensdum = (tensdum-mi)/(ma-mi);
    tensdum(find(tensdum>1)) = 1;
    tensdum(find(tensdum<0)) = 0;
    %im3D(:,:,2) = tensdum;    
    
    %subplot(2,2,[1 3])
    cla    

%     image(tensdum*64),    
%     colormap([zeros(64,1) (0:63)'/63 zeros(64,1)])
    
    imagesc(tensdum), colormap gray
    
    title(num2str(f))

    %imagesc(tensdum,[mi ma]), colormap gray
    

%     hold on
%     
%     id = find(CoMx2{twin(f)} >= xwin(1) & CoMy2{twin(f)} >= ywin(1) & CoMx2{twin(f)} <= xwin(end) & CoMy2{twin(f)} <= ywin(end));
%     plot(CoMx2{twin(f)}(id)-xwin(1)+1,CoMy2{twin(f)}(id)-ywin(1)+1,'.r','markersize',5)
    
%     hold on
%     plot(CoMx2{twin(f)}(pID(1))-xwin(1)+1,CoMy2{twin(f)}(pID(1))-ywin(1)+1,'oy','markersize',20)
%     hold on
%     plot(CoMx2{twin(f)}(pID(2))-xwin(1)+1,CoMy2{twin(f)}(pID(2))-ywin(1)+1,'or','markersize',20)
%     
    
    axis image
%     xlim([xwin(1) xwin(end)])
%     ylim([ywin(1) ywin(end)])    
    axis off
    %drawnow
    
    %%%%%If I've loaded the image of the traces    
%     subplot(2,2,2)
%     hold on, plot([tdom(f) tdom(f)],[-3 -2],'-r','linewidth',2)
%     subplot(2,2,4)
%     hold on, plot([tdom(f) tdom(f)],[-3 -2],'-y','linewidth',2)
%     drawnow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Take screen shot 
    F(f) = getframe(gcf); %Take screen shot   
    
   
    pause(1/fps/8)
end

%% Convert to "indexed" from "truecolor" for compression

v = VideoWriter('cartoonMovie.avi')
open(v)

[fum clrmap] = rgb2ind(rand(1024,1024,3),2^8);
clear Fidx
for f = 1:length(F)
    [Fidx(f).cdata] = rgb2ind(F(f).cdata,clrmap);
    Fidx(f).colormap = [];
    
    %[Fidx(f).cdata] = rgb2ind(F(f).cdata,2^8);  %This doesn't work when I compress
    
    writeVideo(v,F(f))
    
end
close(v)

%%

v = VideoWriter('revcoMovie.avi')
open(v)
writeVideo(v,CH2)

%%
movie2avi(Fidx,'NHPMovie16.avi','fps',round(fps*2),'colormap',clrmap,'compression','MSVC')


movie2avi(Fidx,'NHPMovie.avi','fps',round(fps*2),'colormap',clrmap,'compression','MSVC')


%%
global botpred botactual toppred topactual xdom
%You can be more independent of revcorrprediction.m if you extract the time
%courses here and then run the rest of this code.

botpred = get(gco,'Ydata');
botactual = get(gco,'Ydata');

toppred = get(gco,'Ydata');
topactual = get(gco,'Ydata');

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;
xdom = (0:length(botpred)-1)*acqPeriod/1000;

%%

figure(20)

subplot(2,2,2)
plot(xdom,toppred,'k')
hold on
plot(xdom,topactual,'r')
axis off
xlabel('seconds'), ylabel('standard deviation')
text(30,4,'prediction','Color',[0 0 0],'FontSize',12)
text(30,3,'actual','Color',[1 0 0],'FontSize',12)
plot([0 5],[2 2],'k')
text(0,2.5,'5 sec','Color',[0 0 0],'FontSize',10)
ylim([-3 4])

subplot(2,2,4)
plot(xdom,botpred,'k')
hold on
plot(xdom,botactual,'y')
axis off
set(gca,'XTick',[0 60])
xlabel('seconds'), ylabel('standard deviation')
text(30,4,'prediction','Color',[0 0 0],'FontSize',12)
text(30,3,'actual','Color',[1 1 0],'FontSize',12)
ylim([-3 4])
%%
% tdom = twin*fp/1000 - getparam('predelay');
% 
% %pID = nID([21 35 54]); %ab2 000_014
% 
% pID = nID([35 54]); %ab2 000_014
% 
% fps = 1000/fp;
% 
% clear F
% %mov = avifile('RawTrial.avi')
% figure(20)
% for f = 1:length(twin)
%     tensdum = CH2(:,:,f);  
%     
%     
%     figure(20)
%     subplot(2,2,[1 3])
%     cla
%     
%     imagesc(tensdum,[mi ma]), colormap gray
%     
% 
%     %hold on
%     
% %     id = find(CoMx2{twin(f)} >= xwin(1) & CoMy2{twin(f)} >= ywin(1) & CoMx2{twin(f)} <= xwin(end) & CoMy2{twin(f)} <= ywin(end));
% %     plot(CoMx2{twin(f)}(id)-xwin(1)+1,CoMy2{twin(f)}(id)-ywin(1)+1,'.r','markersize',5)
%     
%     hold on
%     plot(CoMx2{twin(f)}(pID(1))-xwin(1)+1,CoMy2{twin(f)}(pID(1))-ywin(1)+1,'ob','markersize',20)
%     hold on
%     plot(CoMx2{twin(f)}(pID(2))-xwin(1)+1,CoMy2{twin(f)}(pID(2))-ywin(1)+1,'or','markersize',20)
%     
%     axis image
% %     xlim([xwin(1) xwin(end)])
% %     ylim([ywin(1) ywin(end)])    
%     axis off
%     %drawnow
%     
%     %%%%%If I've loaded the image of the traces    
%     subplot(2,2,2)
%     hold on, plot([tdom(f) tdom(f)],[-4 -3],'-b','linewidth',2)
%     subplot(2,2,4)
%     hold on, plot([tdom(f) tdom(f)],[-4 -3],'-r','linewidth',2)
%     drawnow
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     
%     F(f) = getframe(gcf);
%     %%mov = addframe(mov,F);
%     
%     pause(1/fps/8)
% end
% movie2avi(F,'RawTrial7.avi','fps',round(fps*2)')
