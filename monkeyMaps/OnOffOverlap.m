function [D Dnorm r] = OnOffOverlap(TC_op)


%% b+w
fnames = fieldnames(TC_op{1});
for e = 1:length(TC_op)
    for i = 1:length(fnames)
        if e == 1
            eval(['TCopAll.' fnames{i} ' = [];']);
        end
        try
            dum = eval(['TC_op{' num2str(e) '}.'  fnames{i} '{1}{1};']); dim = size(dum);
            if min(dim) == 1
                eval(['TCopAll.' fnames{i} '= [TCopAll.' fnames{i} '; dum(:)];']);
                %         else
                %             eval([fnames{i} '= [' fnames{i} '; dum];']);
            end
        catch
            
        end
    end
end

%% black

fnames = fieldnames(TC_op{1});
for e = 1:length(TC_op)
    for i = 1:length(fnames)
        if e == 1
            eval(['TCopAll_black.' fnames{i} ' = [];']);
        end
        try
            dum = eval(['TC_op{' num2str(e) '}.'  fnames{i} '{1}{2};']); dim = size(dum);
            if min(dim) == 1
                eval(['TCopAll_black.' fnames{i} '= [TCopAll_black.' fnames{i} '; dum(:)];']);
                %         else
                %             eval([fnames{i} '= [' fnames{i} '; dum];']);
            end
        catch
            
        end
    end
end

%% white

fnames = fieldnames(TC_op{1});
for e = 1:length(TC_op)
    for i = 1:length(fnames)
        if e == 1
            eval(['TCopAll_white.' fnames{i} ' = [];']);
        end
        try
            dum = eval(['TC_op{' num2str(e) '}.'  fnames{i} '{1}{3};']); dim = size(dum);
            if min(dim) == 1
                eval(['TCopAll_white.' fnames{i} '= [TCopAll_white.' fnames{i} '; dum(:)];']);
                %         else
                %             eval([fnames{i} '= [' fnames{i} '; dum];']);
            end
        catch
            
        end
    end
end

%%

for p = 1:length(TCopAll.xpos)
    
    %Get a single orientation based on the b+w average
    im = TCopAll.RFraw{p};
    [bestoriid bestposid] = find(im == max(im(:)));
    
    im_b = TCopAll_black.RFraw{p};
    im_w = TCopAll_white.RFraw{p};
    degperpix = getparam('x_size')/(length(im_w(1,:))-1);  %can't use xdom because of the padding    
    
          
    rawProfile_b = im_b(bestoriid,:);
    [param profileFit_b varaccount1D_b] = Gaussfit(1:length(rawProfile_b),rawProfile_b,0);
    param(2) = abs(param(2));  %This is the same since it gets squared
    sig_b = param(2)*degperpix;
    mu_b = param(1)*degperpix;
    
    rawProfile_w = im_w(bestoriid,:);
    [param profileFit_w varaccount1D_w] = Gaussfit(1:length(rawProfile_w),rawProfile_w,0);
    param(2) = abs(param(2));  %This is the same since it gets squared
    sig_w = param(2)*degperpix;
    mu_w = param(1)*degperpix;
    
    %if ~isnan(TCopAll_black.xpos(p)) & ~isnan(TCopAll_white.xpos(p)) & abs(mu_w-mu_b)<1
    if    abs(mu_w-mu_b)<1
    %if varaccount1D_b >.8 & varaccount1D_w >.8
        D(p) = abs(mu_w-mu_b);
        Dnorm(p) = D(p)/(sig_w+sig_b);
        rdum = corrcoef(profileFit_b,profileFit_w);
        r(p) = rdum(1,2);
   else
        D(p) = NaN;
        Dnorm(p) = NaN;
        r(p) = NaN;
   end
    
end

