function SNR = getSNR

global f0m f0m_var Tens Tens_var

dim = size(f0m{1}(:,:,1));
f0varTens = zeros(dim(1),dim(2),length(f0m),'single'); %preallocate
f0Tens = f0varTens;
for k = 1:length(f0m)
    f0varTens(:,:,k) = f0m_var{k};
    f0Tens(:,:,k) = f0m{k};
end

[dum idma] = max(f0Tens,[],3); %id condition with best response
[dum idmi] = min(f0Tens,[],3);

M = size(f0m{1}(:,:,1),1); N = size(f0m{1}(:,:,1),2); 
maR = zeros(size(f0m{1}(:,:,1)));
miR = zeros(size(f0m{1}(:,:,1)));
maRvar = zeros(size(f0m{1}(:,:,1)));
miRvar = zeros(size(f0m{1}(:,:,1)));
for i = 1:length(f0m)

    id = find(idma == i);
    [pk pkid] = max(Tens{i},[],3);
    maR(id) = pk(id);
    
    pkidvec{i} = (pkid-1)*M*N + reshape((1:(M*N))',M,N);
    
    maRvar(id) = Tens_var{i}(pkidvec{i}(id));

end

for i = 1:length(f0m)

    id = find(idmi == i);  %get the pixels for whic this cond is its min
    miR(id) = Tens{i}(pkidvec{i}(id));
    
    miRvar(id) = Tens_var{i}(pkidvec{i}(id));

end
maRvar = sqrt(maRvar/getnorepeats(1));
miRvar = sqrt(miRvar/getnorepeats(1)); %standard Error

SNR = (maR-miR)./(maRvar+miRvar);
