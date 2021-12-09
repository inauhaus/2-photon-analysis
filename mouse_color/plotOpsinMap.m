function plotOpsinMap(percS,rfit)

%Ian Nauhaus

global TC DM MK ACQinfo maskS idExamp cellS

MK = struct; %Reset these

%Make/store another structure related to the maskS, 'MK'
MK.Ncell = size(cellS.muTime{1},1) 

MK.masklabel = bwlabel(maskS.bwCell{1},4);
MK.celldom = unique(MK.masklabel);
[MK.nID] = getNeuronMask;

[xmicperpix ymicperpix] = getImResolution;

%Now plot color image of tuning

SIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);

SIm_hat = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;
%%
Srange = [0 100];
%sfdom = logspace(log10(sfrange(1)),log10(sfrange(end)),5);
%Sdom = DM.sfdom;
percS(find(rfit<.7)) = NaN;

for p = 2:MK.Ncell
    
    idcell = find(MK.masklabel(:) == p);
    
    SIm(idcell) = percS(p-1);
    
    %SIm_hat(idcell) = S_hat(p);
    
end

figure,imagesc(SIm,'AlphaData',sign(SIm),[0 100])
% 
% %%%Spatial frequency
% figure
% IMtens = getImTens(log2(sfprefIm),sign(sfmagIm),log2(sfrange),'jet',1);
% image(xdom,ydom,IMtens), colorbar, colormap jet
% for i = 1:length(sfdom)
%     domcell{i} = round(sfdom(i)*100)/100;
% end
% 
% sfvec = round(log2(sfdom)*100)/100;
% sfvec(end) = floor(log2(sfdom(end))*100)/100;
% sfvec(1) = ceil(log2(sfdom(1))*100)/100;
% sfvec = (sfvec - log2(sfrange(1)))/(log2(sfrange(2))-log2(sfrange(1)));
% %sfvec = round(sfvec*63+1);
% colorbar('YTick',sfvec,'YTickLabel',domcell)
% title(['color ' num2str(c)]), axis image
% plotEXcirc(idExamp,xdom,ydom)
% hold on
% hyp = 10; orig = -5;
% plot([orig hyp*cos(sfprincax*pi/180)+orig],[orig hyp*sin(sfprincax*pi/180)+orig],'k','Clipping','off')
% %set(gca,'Position',[.13 .1 xysize])
% 
% subplot(2,2,4)
% IMtens = getImTens(log2(sfprefIm_hat),sign(sfmagIm),log2([.5 8]),'jet',1);
% image(xdom,ydom,IMtens)
% %colorbar('YTick',sfvec,'YTickLabel',domcell)
% axis image
% plotEXcirc(idExamp,xdom,ydom)
% 
% %set(gca,'Position',[.53 .1 xysize])
% 
% 
% 
% %% Anatomy plot
% CH = GetTrialData([1 0],5);
% im{1} = mean(CH{1}(:,:,2:end-1),3);
% dum = im{1};
% dum = LocalZ(dum,150);
% mi = prctile(dum(:),2);
% ma = prctile(dum(:),99.9);
% dum(find(dum>ma)) = ma;
% dum(find(dum<mi)) = mi;
% 
% figure,imagesc(xdom,ydom,dum), colormap gray, axis image
% %plotEXcirc(idExamp,xdom,ydom)
% %Need circles to be red, not black
% if ~isempty(idExamp)
%     for q = 1:length(idExamp)
%         hold on
%         plot(xdom(round(MK.CoM(idExamp(q),2))),ydom(round(MK.CoM(idExamp(q),1))),'or','MarkerSize',15)
%     end
% end
%%

function dist = oridiff(angle1,angle2)

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;

function plotEXcirc(idExamp,xdom,ydom)

global MK

if ~isempty(idExamp)
    for q = 1:length(idExamp)
        hold on
        plot(xdom(round(MK.CoM(idExamp(q),2))),ydom(round(MK.CoM(idExamp(q),1))),'ok','MarkerSize',10)
    end
end

function IMtens = getImTens(pref,mag,mima,maptype,BackG)

id = find(mag>1);
mag(id) = 1;

id = find(pref<mima(1));
pref(id) = mima(1);
id = find(pref>mima(2));
pref(id) = mima(2);


prefid = (pref-mima(1))/(mima(2)-mima(1));
prefid = round(prefid*63+1);  %normalize to be colormap index

dim = size(pref);
mapvals = eval(maptype);
IMtens = zeros(dim(1),dim(2),3);
for i = 1:dim(1)
    for j = 1:dim(2)
        if mag(i,j) == 0 | isnan(mag(i,j)) | isnan(prefid(i,j));
            IMtens(i,j,:) = [BackG BackG BackG];
        else
            IMtens(i,j,:) = mag(i,j)*mapvals(prefid(i,j),:);
        end
    end
end


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