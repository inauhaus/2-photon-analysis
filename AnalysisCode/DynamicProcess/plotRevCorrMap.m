function orimap = plotRevCorrMap

%Ian Nauhaus

global kernelsIm G_handles maskS G_RChandles ACQinfo

[xmicperpix ymicperpix] = getImResolution;
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  

%Compute the orimap

dim = size(kernelsIm);
dim2 = size(kernelsIm{1,1,1}); 

eval(['t1_t2 =' get(G_RChandles.STMovieTimeWindow,'string') ';']);
t1 = t1_t2(1); t2 = t1_t2(2);
[dum id1] = min(abs(taudom-t1));
[dum id2] = min(abs(taudom-t2));

%Get image from time zero
k = 0;
imT0 = 0;
for ori = 1:dim(1)    
    for sf = 1:dim(2)
        for phase = 1:dim(3)
            if ~isempty(kernelsIm{ori,sf,phase})
                imT0 = imT0 + kernelsIm{ori,sf,phase}(:,:,1);
                k = k+1;
            end
        end
    end
end
imT0 = squeeze(imT0/k);

Tens = zeros(dim2(1),dim2(2),dim(1));
oricounter = zeros(1,dim(1));
for ori = 1:dim(1)    
    for sf = 2:dim(2)-1
        for phase = 1:dim(3)
            if ~isempty(kernelsIm{ori,sf,phase})
                Tens(:,:,ori) = Tens(:,:,ori) + mean(kernelsIm{ori,sf,phase}(:,:,id1:id2),3);
                %Tens(:,:,ori) = Tens(:,:,ori)- mean(kernelsIm{ori,sf,phase}(:,:,1),3);
                oricounter(ori) = oricounter(ori)+1;
               
            end
        end
    end    
end
for ori = 1:dim(1)   
    Tens(:,:,ori) = Tens(:,:,ori)/oricounter(ori);
end

% for ori = 1:dim(1)   
%     Tens(:,:,ori) = (Tens(:,:,ori)-imT0)./imT0;
% end

% mi = min(Tens,[],3);
% for ori = 1:dim(1)   
%     Tens(:,:,ori) = Tens(:,:,ori)-mi;
% end
% su = sum(Tens,3);
% for ori = 1:dim(1)   
%     Tens(:,:,ori) = Tens(:,:,ori)./su;
% end

dori = 180/dim(1);
oridom = 0:dori:180-dori;
orimap = zeros(dim2(1),dim2(2));

for i = 1:dim(1)    
    orimap = orimap + Tens(:,:,i)*exp(1i*2*oridom(i)*pi/180);    
end

dum = sum(Tens,3);
%orimap = orimap./dum;

id = find(isnan(orimap(:)));
orimap(id) = 0;

h = fspecial('gaussian', size(orimap), 1);
h = abs(fft2(h));
orimap = ifft2(fft2(orimap).*h);
dum = ifft2(fft2(dum).*h);


mag = abs(orimap);
ang = angle(orimap)*180/pi;
ang = (ang + (1-sign(ang))*180)/2;

% [ma id] = max(Tens,[],3);
% [mi dum] = min(Tens,[],3);
% mag = (ma-mi);

%This is because the edges are wierd...
%mag = mag(3:end-2,3:end-2); ang = ang(3:end-2,3:end-2);

mag = mag.^2;

ma = prctile(mag(:),99.5);
mag(find(mag>ma)) = ma;
mi = prctile(mag(:),.5);
mag(find(mag<mi)) = mi;

mag = mag-min(mag(:));
mag = mag/max(mag(:));

%mag = ones(size(mag));

%Plot the map with anatomy
figure
dim = size(ang);
set(gcf,'Color',[1 1 1]);

anatflag = 0;

if anatflag
    
    imanat = maskS.im{1};
    imanat = imanat(3:end-2,3:end-2);
   
    mi = prctile(imanat(:),0);    
    imanat = phi(imanat-mi);
    ma = prctile(imanat(:),100);
    id = find(imanat>ma);
    imanat(id) = ma;
    imanat = imanat/ma;
    
    mag = sqrt(imanat.*mag);

    imfunc = ang;
    imfunc = imfunc/180;
    imfunc = round(imfunc*63+1);
    %imanat = round(imanat*63+1);

    hsvid = hsv;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)            
            imout(i,j,:) = mag(i,j)*hsvid(imfunc(i,j),:);
        end
    end
    
    imanat(:,:,2) = imanat;
    imanat(:,:,3) = imanat(:,:,1);
    
    imout = imout+(imanat).^.3;

    imout = imout/max(imout(:));
    
    %imout = imout(5:end,1:end-7,:);
    
    x = image(xdom,ydom,imout,'CDataMapping','direct','AlphaDataMapping','none');

else
    imout = ang;
    imout = imout/180;
    imout = round(imout*63+1);
    x = image(xdom,ydom,imout,'CDataMapping','direct','AlphaData',mag,'AlphaDataMapping','none');
    
end

axis image;

fh = gcf;

colormap hsv
%colorbar('YTick',[1 16:16:64],'YTickLabel',{'0','45','90','135','180'})

makeOriLegend(xdom,ydom)


function makeOriLegend(xdom,ydom)
%%
hold on

%Create the orientation legend%%%%%%%%%%%%%%%%%%
legdom = 0:30:180;
hsvdom = hsv;
id = round(linspace(1,64,length(legdom)));
hsvdom = hsvdom(id,:);
R = 14;
rid = linspace(ydom(1),ydom(end),length(legdom));
cid = xdom(end)+15;
xpts_o = [0 0];
ypts_o = [1-R 1+R];

for i = 1:length(legdom)
   
    xpts = xpts_o*cos(legdom(i)*pi/180) + ypts_o*sin(legdom(i)*pi/180);
    ypts = xpts_o*sin(legdom(i)*pi/180) - ypts_o*cos(legdom(i)*pi/180);
    ypts = ypts + rid(i);
    xpts = xpts + cid;
    hold on
    line(xpts,ypts,'Color',hsvdom(i,:),'Clipping','off','LineWidth',3);
    
end

xlim([xdom(1) xdom(end)])
ylim([ydom(1) ydom(end)])
hold off

%%%%%%%%%%%%%%%