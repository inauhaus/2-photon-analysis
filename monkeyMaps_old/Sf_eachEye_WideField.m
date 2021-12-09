function [sfpref_raw ocdommap_raw dprimemask] = Sf_eachEye_WideField(dpthresh)

global bw f0m f0m_var funcmap ACQinfo symbolInfo Analyzer G_handles idExamp

[xmicperpix ymicperpix] = getImResolution(1);

set(G_handles.HPflag,'Value',0);
set(G_handles.LPflag,'Value',1);

set(G_handles.Lwidth,'string','8');
hh = makeMapFilter;
set(G_handles.Lwidth,'string','.1');
hh2 = makeMapFilter;

anatomyflag = 0;
bwCellPlot = ones(size(funcmap));

sfdom = getdomain('s_freq');

%% get sfmap first eye (usually left eye)

for i = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp(Analyzer.loops.conds{1}.symbol{i},'s_freq');
        idsym = i;
        break
    end
end
for i = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp(Analyzer.loops.conds{1}.symbol{i},'Leye_bit');
        idsym2 = i;
        break
    end
end

symbolInfo.ID(1) = idsym;
set(G_handles.primSymbol,'value',idsym); 

symbolInfo.ID(2) = idsym2;
set(G_handles.secSymbol,'value',idsym2); 

if idsym + idsym2 == 3 
    set(G_handles.tertSymbol,'value',3);
elseif idsym + idsym2 == 4 
    set(G_handles.tertSymbol,'value',2);
elseif idsym + idsym2 == 5
    set(G_handles.tertSymbol,'value',1);
end

setsymbolstruct

funcmap = GprocessLog(f0m,bwCellPlot,hh,1);   %output is complex
sfmag1 = real(funcmap);
sfpref1 = imag(funcmap);

funcmap = GprocessLog(f0m,bwCellPlot,hh2,1);   %output is complex
sfmag_raw1 = real(funcmap);
sfpref_raw1 = imag(funcmap);

% sfpref_raw = medfilt2(sfpref_raw,[3 3]);
% sfpref = medfilt2(sfpref,[3 3]);
% 
% id = find(sfpref_raw<0.5);
% sfpref_raw(id) = 0.5;
% id = find(sfpref_raw>8);
% sfpref_raw(id) = 8;
% 
% id = find(sfpref<0.5);
% sfpref(id) = 0.5;
% id = find(sfpref>8);
% sfpref(id) = 8;


%% get sfmap second eye (usually right eye)


funcmap = GprocessLog(f0m,bwCellPlot,hh,0);   %output is complex
sfmag2 = real(funcmap);
sfpref2 = imag(funcmap);

funcmap = GprocessLog(f0m,bwCellPlot,hh2,0);   %output is complex
sfmag_raw2 = real(funcmap);
sfpref_raw2 = imag(funcmap);


%% Get mask for data selection

Nh = 9;
h = zeros(size(f0m{1}));
dum = fspecial('disk', Nh/2);
h(1:size(dum,1),1:size(dum,2)) = dum/sum(dum(:));
%h = hh;


dim = size(f0m{1});
Tens = zeros(dim(1),dim(2),length(f0m),'single'); %preallocate
Tens_var = Tens;
for k = 1:length(f0m)
    
    id = find(isnan(f0m{k}) | isnan(f0m_var{k}));
    f0m{k}(id) = min(f0m{k}(:));
    f0m_var{k}(id) = min(f0m_var{k}(:));
    
    if ~isempty(hh2)
        Tens(:,:,k) = ifft2(fft2(f0m{k}).*abs(fft2(h)));
        Tens_var(:,:,k) = ifft2(fft2(f0m_var{k}).*abs(fft2(h)));
    else
        Tens(:,:,k) = f0m{k};
        Tens_var(:,:,k) = f0m_var{k};
    end
    
end

hsum = sum(h(:)/max(h(:)));

[pk_dF idpk] = max(Tens(:,:,1:end-1),[],3);
pkSE_dF = zeros(size(idpk));
for i = 1:size(Tens_var,1)
    for j = 1:size(Tens_var,2)
        pkSE_dF(i,j) = sqrt(Tens_var(i,j,idpk(i,j)))/sqrt(getnorepeats(1))/sqrt(hsum);  %standard error of best response at each pixel
    end
end

[mi_dF idmi] = min(Tens(:,:,1:end-1),[],3);
miSE_dF = zeros(size(idmi));
for i = 1:size(Tens_var,1)
    for j = 1:size(Tens_var,2)
        miSE_dF(i,j) = sqrt(Tens_var(i,j,idmi(i,j)))/sqrt(getnorepeats(1))/sqrt(hsum);  %standard error of best response at each pixel
    end
end

base_dF = Tens(:,:,end); %response to blank
baseSE_dF = sqrt(Tens_var(:,:,end))/sqrt(getnorepeats(getnoconditions))/sqrt(hsum);  %standard error for blank


%%
%dprime = pk_dF./pkSE_dF;
%dprime = pk_dF./Tens_var(:,:,end);

dprime = (pk_dF-base_dF)./(pkSE_dF+baseSE_dF);

% dprime = (pk_dF-mi_dF)./(pkSE_dF+miSE_dF);

dprime = phi(dprime);

id = find(Tens_var(:,:,1)<=0 | imag(dprime)~=0);
dprime(id) = 0;

dprime = ifft2(fft2(dprime).*abs(fft2(hh)));  %smooth the mask
dprime(:,[1:8 end-7:end]) = 0; dprime([1:8 end-7:end],:) = 0;
% 
% dprime = dprime-median(dprime(:));
% dprime = dprime*15;
% dprime = 1./(1+exp(-dprime));

dprimemask = zeros(size(dprime));
dprimemask(find(dprime>dpthresh)) = 1;
dprimemask(:,[1:6 end-5:end]) = 0; dprimemask([1:6 end-5:end],:) = 0;

%Fill holes
%dprimemask = imfill(dprimemask);

%Get rid of outside ROIs ("small islands")
L = bwlabel(dprimemask);
Ldom = unique(L(:));
Ldom = Ldom(2:end);
clear Npix
for i = 1:length(Ldom)
    Npix(i) = length(find(L == Ldom(i)));
end
id = find(Npix<max(Npix));
for i = 1:length(id)
    idx = find(L == Ldom(id(i)));
    dprimemask(idx) = 0;
end

set(G_handles.Lwidth,'string','1');
h = makeMapFilter;
dprime2 = ifft2(fft2(double(dprimemask)).*abs(fft2(h)));
dprime2(:,1) = -.3;

%% plots

%ocdommap_raw = ocdommap_raw - median(ocdommap_raw(:));
%ocdommap = ocdommap - median(ocdommap(:));

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

logdom = logspace(log10(sfdom(1)),log10(sfdom(end)),25);

figure
subplot(1,3,1)
leg = {'ow','sw','dw','pw'};
sfdomdum = log(logspace(log10(sfdom(1)),log10(sfdom(end)),15));

sfpref_raw1(1,1) = sfdom(1); sfpref_raw1(1,1) = sfdom(end); 
plotlogmap(dprime2,sfpref_raw1,0,xdom,ydom)
hold on
sfprefdum = log(sfpref2);
sfprefdum = (sfprefdum-min(sfprefdum(:)))/(max(sfprefdum(:))-min(sfprefdum(:)));
rng = 50; sfprefdum = sfprefdum*rng;
sfdomdum = rng*(sfdomdum-log(sfdom(1)))/(log(sfdom(end))-log(sfdom(1)));
contour(xdom,ydom,sfprefdum,sfdomdum,'k')
hold on
if ~isempty(idExamp)
    for i = 1:length(idExamp2(:,1))
        plot(idExamp2(i,1),idExamp2(i,2),leg{i},'markersize',10,'linewidth',2)
        hold on
    end
end

%[idy idx] = find(1-bw);  %cover up outside the ROI
%plot(idx*xmicperpix,idy*ymicperpix,'.w','MarkerSize',15)  
%hold off
title('Sfreq map (eye 1) w/ Sfreq contour (eye 2)')

clear domcell
for i = 1:length(sfdom)
    domcell{i} = round(sfdom(i)*100)/100;
end
iddom = linspace(min(sfdomdum),max(sfdomdum),length(sfdom));
colorbar('YTick',iddom,'YTickLabel',domcell)
axis off
axis image
 
subplot(1,3,2)
leg = {'ow','sw','dw','pw'};
sfdomdum = log(logspace(log10(sfdom(1)),log10(sfdom(end)),15));

sfpref_raw2(1,1) = sfdom(1); sfpref_raw2(1,1) = sfdom(end); 
plotlogmap(dprime2,sfpref_raw2,0,xdom,ydom)
hold on
sfprefdum = log(sfpref1);
sfprefdum = (sfprefdum-min(sfprefdum(:)))/(max(sfprefdum(:))-min(sfprefdum(:)));
rng = 50; sfprefdum = sfprefdum*rng;
sfdomdum = rng*(sfdomdum-log(sfdom(1)))/(log(sfdom(end))-log(sfdom(1)));
contour(xdom,ydom,sfprefdum,sfdomdum,'k')
hold on
if ~isempty(idExamp)
    for i = 1:length(idExamp2(:,1))
        plot(idExamp2(i,1),idExamp2(i,2),leg{i},'markersize',10,'linewidth',2)
        hold on
    end
end

%[idy idx] = find(1-bw);  %cover up outside the ROI
%plot(idx*xmicperpix,idy*ymicperpix,'.w','MarkerSize',15)  
%hold off
title('Sfreq map (eye 2) w/ Sfreq contour (eye 1)')

clear domcell
for i = 1:length(sfdom)
    domcell{i} = round(sfdom(i)*100)/100;
end
iddom = linspace(min(sfdomdum),max(sfdomdum),length(sfdom));
colorbar('YTick',iddom,'YTickLabel',domcell)
axis off
axis image


%plot both ocdom/sf contours together
subplot(1,3,3)
contour(xdom,ydom,sfpref1,logdom,'k')
hold on
contour(xdom,ydom,sfpref2,logdom,'r')
axis ij
hold on
[idy idx] = find(1-dprimemask);  %cover up outside the ROI
plot(idx*xmicperpix,idy*ymicperpix,'.w','MarkerSize',15)  
axis image
title('Sfreq contour (eye 1) & Sfreq contour (eye 2)')
set(gca,'Xtick',[],'Ytick',[])
hold on,  line([20 20],[0 100])
hold off


%%
idROI = find(dprimemask>0);
idROI = idROI(1:5:end);

figure,
subplot(1,2,1)
scatter(sfpref_raw1(idROI),sfpref_raw2(idROI),'.k')
hold on, plot([0 5],[0 5],'r')
[r p] = corrcoef(sfpref_raw1(idROI),sfpref_raw2(idROI))
xlabel('SF preference (eye 1)'), ylabel('SF preference (eye 2)')
title(['r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2))])
xlim([0 5]), ylim([0 5]), axis square

subplot(1,2,2), hist(log2(sfpref_raw2(idROI)./sfpref_raw1(idROI)),20)
xlabel('log(sf2/sf1)')
ylabel('N pixels')
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



%%
function plotODmap(mag,pref,anatflag,xdom,ydom,varargin)

global fh symbolInfo

if ~isempty(varargin)
    transparentMask = varargin{1};
    transparentID = find(transparentMask == 1);
end

mag = phi(mag-prctile(mag(:),0));
mag = mag/prctile(mag(:),98);
mag(find(mag>1)) = 1;

%pref = log2(pref);
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


