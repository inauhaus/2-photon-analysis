function plotODmap(mag,pref,anatflag,xdom,ydom,varargin)

global fh symbolInfo


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
    if isempty(varargin)
        imfunc = imfunc-min(imfunc(:));
        imfunc = imfunc/max(imfunc(:));
    else
        imfunc = (imfunc-varargin{1}(1))/(varargin{1}(2)-varargin{1}(1));      
        id = find(imfunc<0);
        imfunc(id) = 0;
        id = find(imfunc>1);
        imfunc(id) = 1;
    end
    
    imfunc = round(imfunc*63+1);

    jetid = jet;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)
            imout(i,j,:) = mag(i,j)*jetid(imfunc(i,j),:);
        end
    end    
    
    
    imout = imout/max(imout(:));
    
    image(xdom,ydom,imout,'CDataMapping','direct'); 
    
       
end
axis image


