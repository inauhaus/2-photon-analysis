function Sep = localSVD2(T)

%input is tensor(x,y,T), output is image(x,y)

global ACQinfo
%%
%T = T(1:2:end,1:2:end,:);
W = 11;
im = zeros(size(T,1),size(T,2));
Sep = zeros(size(T,1),size(T,2));
dim = size(T);
T = T-median(T(:));
for i =ceil(W/2):5:(dim(1)-floor(W/2))
    for j = ceil(W/2):5:(dim(2)-floor(W/2))
        
        yran = i-floor(W/2):i+floor(W/2);
        xran = j-floor(W/2):j+floor(W/2);
%         yran(find(yran<1|yran>dim(1))) = [];
%         xran(find(xran<1|xran>dim(2))) = [];
        blck = T(yran,xran,:);
        bdim = size(blck);
        bvec = reshape(blck,[bdim(1)*bdim(2) bdim(3)])';
        %bvec = bvec-median(bvec(:));
        %bvec = bvec-ones(size(bvec,1),1)*mean(bvec);
        %bvec = bvec./(ones(size(bvec,1),1)*std(bvec));
        [U,S,V] = svd(bvec','econ');
        Sep(i,j) = S(1,1)/sum(S(:));
        
        % [U,S,V] = svd(To','econ');
        % C1 = U(:,1)*V(:,1)'*S(1,1);
        % To = To-C1';
        % To = reshape(To',[dim(1) dim(2) dim(3)]);
        
        C1 = U(:,1)*S(1,1);% + U(:,2)*S(2,2)+ U(:,3)*S(3,3); %U is space;
        C1 = reshape(C1',[W W]);

        im(yran,xran) = im(yran,xran) + C1;
        
        
    end
    i
end




