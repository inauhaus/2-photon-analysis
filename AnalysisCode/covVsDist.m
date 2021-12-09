function covVsDist

sig = tcourseX{1};
sig = sig-ones(size(sig,1),1)*mean(sig);
sig = sig./(ones(size(sig,1),1)*std(sig));
%sig = signew;

figure,plot(sig)

cov = sig'*sig/size(sig,1);
figure,imagesc(cov)



maskID =  bwlabel(maskS.bwCell{1},4);
id = unique(maskID(:));

xdom = 1:size(maskS.bwCell{1},2);
ydom = 1:size(maskS.bwCell{1},1);

clear CoMyx 
for i = 1:length(id)
    dum = zeros(size(maskID));
    dum(find(maskID == id(i))) = 1;
    CoMyx(i,1) = round(ydom*sum(dum')'/sum(dum(:)));
    CoMyx(i,2) = round(xdom*sum(dum)'/sum(dum(:)));
end

CoMyx(1,:) = [];

Ncell = size(sig,2);

[xmicperpix ymicperpix] = getImResolution;

        
k = 1;
clear D cov
for i = 1:Ncell
    i
    for j = i+1:Ncell
        
        D(k) = sum((CoMyx(i,:)-CoMyx(j,:)).^2).^.5;
        cov(k) = sig(:,i)'*sig(:,j)/size(sig,1);
        
        k=k+1;
    end
end

figure,scatter(D*ymicperpix,cov,'.')
xlabel('microns'), ylabel('r'), xlim([0 800]), ylim([-1 1])
