function [magodgrad dphaseROI dphaseROISEL] = plotSFODstuff(sfpref,sfpref_raw,ocdommap,ocdommap_raw,dprime2,dprimemask,trX,trY)

global ACQinfo G_handles idExamp SelectivityThresh

[xmicperpix ymicperpix] = getImResolution(1);

sfdom = getdomain('s_freq');

N_SFlines = 20;
N_ODlines = 10;

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

if ~isempty(idExamp)
    idExamp2(:,1) = idExamp(:,1)/ACQinfo.pixelsPerLine*xdom(end);
    idExamp2(:,2) = idExamp(:,2)/ACQinfo.pixelsPerLine*ydom(end);
end

figure
subplot(1,3,1)
leg = {'ow','sw','dw','pw'};

mi = prctile(ocdommap_raw(:),.5);
ma = prctile(ocdommap_raw(:),99.5);
mi = max([-1 mi]);
ma = min([1 ma]);
mima = max([-mi ma]);

plotODmap(dprime2(trY,trX),ocdommap_raw(trY,trX),0,xdom(trX),ydom(trY),[-mima mima])

%odbardom = linspace(min(ocdommap_raw(:)),max(ocdommap_raw(:)),5);
odbardom = linspace(-mima,mima,5);
clear domcell
for i = 1:length(odbardom)
    domcell{i} = round(odbardom(i)*1000)/1000;
end

iddom = linspace(1,64,length(odbardom));
colorbar('YTick',iddom,'YTickLabel',domcell)

axis([0 1000 0 1000])  %This is to make all experiments scaled appropriately, relative to each other
axis off
hold on, plot([0 500],[995 995],'k')  %Plot scale bar
text(0,1050,'500um') 

hold on
if ~isempty(idExamp)
    for i = 1:length(idExamp2(:,1))
        plot(idExamp2(i,1),idExamp2(i,2),leg{i},'markersize',10,'linewidth',2)
        hold on
    end
end

% hold on
% sfprefdum = log(sfpref);
% mi = min(sfprefdum(:)); ma = max(sfprefdum(:));
% sfprefdum = (sfprefdum-mi)/(ma-mi);
% sfdomdum = log(logspace(log10(sfdom(1)),log10(sfdom(end)),N_SFlines));
% sfdomdum = (sfdomdum-mi)/(ma-mi); 
% contour(xdom(trX),ydom(trY),sfprefdum(trY,trX),sfdomdum,'k')
% title('OD map w/ SF contour')


%plot sfmap w/ ocdom contour
sfpref_raw(1,1) = sfdom(1); sfpref_raw(1,end) = sfdom(end); 
subplot(1,3,2)
plotlogmap(dprime2(trY,trX),sfpref_raw(trY,trX),0,xdom(trX),ydom(trY))

hold on

if ~isempty(idExamp)
    for i = 1:length(idExamp2(:,1))
        plot(idExamp2(i,1),idExamp2(i,2),leg{i},'markersize',10,'linewidth',2)
        hold on
    end
end

%Make colorbar
sfbardom = sfdom;
clear domcell
for i = 1:length(sfbardom)
    domcell{i} = sfbardom(i);
end
iddom = linspace(1,64,length(sfbardom));
colorbar('YTick',iddom,'YTickLabel',domcell)

axis([0 1000 0 1000])  %This is to make all experiments scaled appropriately, relative to each other
axis off



hold on
mi = min(ocdommap(:));
ma = max(ocdommap(:));
ODdum = .20*(ocdommap-mi)/(ma-mi) + sfdom(1);
zro = -.2*mi/(ma-mi) + sfdom(1);  %keep track of the zero crossing
ocdomdom = .20*linspace(0,1,N_ODlines) + sfdom(1);
contour(xdom(trX),ydom(trY),ODdum(trY,trX),ocdomdom,'k')
hold on
contour(xdom(trX),ydom(trY),ODdum(trY,trX),[zro zro],'k','linewidth',2)
%plot both ocdom/sf contours together
logdom = logspace(log10(sfdom(1)),log10(sfdom(end)),N_SFlines);
title('SF map and OD contour')

subplot(1,3,3)
binocmap = 1-abs(ocdommap_raw(trY,trX));

% [dody dodx] = gradient(ocdommap);
% binocmap = sqrt(dody.^2 + dodx.^2);
% binocmap = medfilt2(binocmap,[5 5]);

% ma = prctile(binocmap(:),99);
% binocmap(find(binocmap>ma)) = ma;
% mi = prctile(binocmap(:),1);
% binocmap(find(binocmap<mi)) = mi;
%binocmap = binocmap-min(binocmap(:));

plotODmap(dprime2(trY,trX),binocmap,0,xdom(trX),ydom(trY),[1-mima 1])
title('Binocularity = 1-abs(OD)')
set(gca,'Xtick',[],'Ytick',[])

hold on
if ~isempty(idExamp)
    for i = 1:length(idExamp2(:,1))
        plot(idExamp2(i,1),idExamp2(i,2),leg{i},'markersize',10,'linewidth',2)
        hold on
    end
end

%Make colorbar
%odbardom = linspace(min(binocmap(:)),max(binocmap(:)),5);
odbardom = linspace(1-mima,1,5);
clear domcell
for i = 1:length(odbardom)
    domcell{i} = round(odbardom(i)*1000)/1000;
end

iddom = linspace(1,64,length(odbardom));
colorbar('YTick',iddom,'YTickLabel',domcell)

axis([0 1000 0 1000])  %This is to make all experiments scaled appropriately, relative to each other
axis off


%% Plot anatomy

CH = GetTrialData([1 0 0 0],1);
imanat = mean(CH{1},3);
h = zeros(size(imanat)); 
%h(1:3,1:3) = [.05 .1 .05; .1 1 .1; .05 .2 .05];
h(1:2,1) = [1 1]';
imanat = ifft2(abs(fft2(h)).*fft2(imanat)); %Get rid of the striping from the bidir scan
imanat = imanat.^1.8;

figure,
%subplot(2,2,4)
imagesc(xdom(trX),ydom(trY),imanat(trY,trX),[prctile(imanat(:),2) prctile(imanat(:),99.5)]), colormap gray

axis image
axis([0 1000 0 1000])  %This is to make all experiments scaled appropriately, relative to each other


hold on
if ~isempty(idExamp)
    for i = 1:length(idExamp2(:,1))
        plot(idExamp2(i,1),idExamp2(i,2),leg{i},'markersize',10,'linewidth',2)
        hold on
    end
end

hold off


%% Compute intersection

% %%
idROI = find(dprimemask>0);

dim = size(sfpref);
%dorix = sfpref(:,3:end) - sfpref(:,1:end-2);
%doriy = sfpref(3:end,:) - sfpref(1:end-2,:);

dsfx = log2(sfpref(:,3:end)./sfpref(:,1:end-2));
dsfy = log2(sfpref(3:end,:)./sfpref(1:end-2,:));
dsfx = [zeros(dim(1),1) dsfx zeros(dim(1),1)];
dsfy = [zeros(1,dim(2)); dsfy; zeros(1,dim(2))];

dodx = ocdommap(:,3:end) - ocdommap(:,1:end-2);  %Gradient needs to be of the smoothed map
dody = ocdommap(3:end,:) - ocdommap(1:end-2,:);
dodx = [zeros(dim(1),1) dodx zeros(dim(1),1)];
dody = [zeros(1,dim(2)); dody; zeros(1,dim(2))];

sfgrad = atan2(dsfy,dsfx);
angodgrad = atan2(dody,dodx);

%Find regions with reliable gradients
set(G_handles.Lwidth,'string',num2str(20/xmicperpix*3.823));
h = makeMapFilter;
sfGradSel = gradientReliability(log(sfpref),h);
odGradSel = gradientReliability(ocdommap,h);

%sfGradSel = sfGradSel-median(sfGradSel(:))+SelectivityThresh;
%odGradSel = odGradSel-median(odGradSel(:))+SelectivityThresh;

idGrad = find(sfGradSel> SelectivityThresh & odGradSel > SelectivityThresh);
idHist = intersect(idGrad,idROI);

%% Compute intersection relation betweeen maps 

dphase = angle(exp(1i*sfgrad).*exp(-1i*angodgrad));  
dphaseROI = dphase(idROI)*180/pi;
dphaseROISEL = dphase(idHist)*180/pi;
plotGradIntersection(dphaseROI,dphaseROISEL);

%% Make density plots

%OD gradient map
odZ = ocdommap-mean(ocdommap(idROI)); %Zscore od map before taking gradient
odZ = odZ/std(odZ(idROI));
dodx = odZ(:,3:end) - odZ(:,1:end-2);  %Gradient needs to be of the smoothed map
dody = odZ(3:end,:) - odZ(1:end-2,:);
dodx = [zeros(dim(1),1) dodx zeros(dim(1),1)];
dody = [zeros(1,dim(2)); dody; zeros(1,dim(2))];

dsfdx = log2(sfpref(:,3:end)) - log2(sfpref(:,1:end-2));  %Gradient needs to be of the smoothed map
dsfdy = log2(sfpref(3:end,:)) - log2(sfpref(1:end-2,:));
dsfdx = [zeros(dim(1),1) dsfdx zeros(dim(1),1)];
dsfdy = [zeros(1,dim(2)); dsfdy; zeros(1,dim(2))];
magsfgrad = sqrt(dsfdx.^2 + dsfdy.^2);

magodgrad = sqrt(dodx.^2 + dody.^2);
[mat xdom ydom] = smoothscatter(magodgrad(idROI),log2(sfpref_raw(idROI)),.0001,.01,[0 .4],[-1 3]);
figure
subplot(1,3,1)
imagesc(xdom,ydom,-(mat)), axis xy, colormap gray, axis square
[r p] = corrcoef(magodgrad(idROI),sfpref_raw(idROI));
xlabel('magnitude of OD gradient (normalized)'), ylabel('SF preference (octaves)')
title(['r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2)*10^100) 'x10^-^1^0^0'])

ocdomdum = 1-abs(ocdommap_raw(idROI));
[mat xd yd] = smoothscatter(ocdomdum,log2(sfpref_raw(idROI)),.001,.01,[0 1],[-1 3]);
subplot(1,3,2),
imagesc(xd,yd,-(mat)), axis xy, colormap gray, axis square
[r p] = corrcoef(ocdomdum,log2(sfpref_raw(idROI)));
xlabel('binocularity'), ylabel('SF preference (octaves)')
title(['r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2)*10^100) 'x10^-^1^0^0'])


[mat xdom ydom] = smoothscatter(magodgrad(idROI),magsfgrad(idROI),.0001,.001,[0 .2],[0 .2]);
subplot(1,3,3)
imagesc(xdom,ydom,-(mat)), axis xy, colormap gray, axis square
[r p] = corrcoef(magodgrad(idROI),magsfgrad(idROI));
xlabel('magnitude of OD gradient (normalized)'), ylabel('magsfgrad')
title(['r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2)*10^100) 'x10^-^1^0^0'])

function plotlogmap(mag,pref,anatflag,xdom,ydom,varargin)

global fh symbolInfo


if ~isempty(varargin)
    transparentMask = varargin{1};
    transparentID = find(transparentMask == 1);
end

mag = phi(mag-prctile(mag(:),0));
mag = mag/prctile(mag(:),98);
mag(find(mag>1)) = 1;

pref = log2(pref);
dim = size(mag);
set(gcf,'Color',[1 1 1]);

id = find(isnan(pref));
mag(id) = 0;
pref(id) = min(pref(:));

if anatflag
    
    [imanat] = getExptMean([1 0 0 0],2);
    imanat = imanat{1};
    
    mi = prctile(imanat(:),1);    
    imanat = phi(imanat-mi);
    ma = prctile(imanat(:),99);
    imanat = imanat/ma;

    imfunc = pref;
    imfunc = imfunc-min(imfunc(:));
    imfunc = imfunc/max(imfunc(:));
    imfunc = round(imfunc*63+1);

    jetid = jet;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)
            imout(i,j,:) = mag(i,j)*jetid(imfunc(i,j),:);
        end
    end
    
    imanat(:,:,2) = imanat;
    imanat(:,:,3) = imanat(:,:,1);
    
    imout = imout+sqrt(imanat);

    imout = imout/max(imout(:));
    
    x = image(xdom,ydom,imout,'CDataMapping','direct','AlphaDataMapping','none');

else      

    imfunc = pref;
    imfunc = imfunc-min(imfunc(:));
    imfunc = imfunc/max(imfunc(:));
    imfunc = round(imfunc*63+1);

    jetid = jet;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)
            imout(i,j,:) = mag(i,j)*jetid(imfunc(i,j),:);
        end
    end    

    if ~isempty(varargin)
        for i = 1:3
            imdum = imout(:,:,i);
            imdum(transparentID) = max(imout(:))*ones(size(transparentID));
            imout(:,:,i) = imdum;
        end
    end
    
    
    imout = imout/max(imout(:));
    
    image(xdom,ydom,imout,'CDataMapping','direct'); 
    
       
end
axis image


