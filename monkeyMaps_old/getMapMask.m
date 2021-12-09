function [dprime2 dprimemask] = getMapMask(dpthresh,hh2,hh,edgeTrunc)

%Get ROI mask based on SNR

global f0m f0m_var G_handles

[xmicperpix ymicperpix] = getImResolution(1);

Nh = 9/xmicperpix*3.823;

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

%dprime = (pk_dF-base_dF)./(pkSE_dF+baseSE_dF);
%dprime = (pk_dF-mi_dF)./(pkSE_dF+miSE_dF);

a = (pk_dF-base_dF);
b = 1./(pkSE_dF+baseSE_dF);
dprime = a.*(b-prctile(b(:),2));

dprime = phi(dprime);

id = find(Tens_var(:,:,1)<=0 | imag(dprime)~=0);
dprime(id) = 0;

dprime = ifft2(fft2(dprime).*abs(fft2(hh)));  %smooth the mask

dprime(:,[1:edgeTrunc end-edgeTrunc+1:end]) = 0; 
dprime([1:edgeTrunc end-edgeTrunc+1:end],:) = 0;

dprimemask = zeros(size(dprime));
dprimemask(find(dprime>dpthresh)) = 1;

se = strel('disk',2);
dprimemask = imopen(dprimemask,se);

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

