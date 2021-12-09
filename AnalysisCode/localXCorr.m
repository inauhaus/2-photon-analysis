function im = localXCorr(T)

%input is tensor(x,y,T), output is image(x,y)

global ACQinfo

A = T(1:2:end,1:2:end,:);
B = T(1:2:end,2:2:end,:);
C = T(2:2:end,1:2:end,:);
D = T(2:2:end,2:2:end,:);

dim = size(A);
A = reshape(A,[dim(1)*dim(2) dim(3)])';
B = reshape(B,[dim(1)*dim(2) dim(3)])';
C = reshape(C,[dim(1)*dim(2) dim(3)])';
D = reshape(D,[dim(1)*dim(2) dim(3)])';

dimV = size(A);
A = A./(ones(dimV(1),1)*sqrt(sum(A.*A)));
B = B./(ones(dimV(1),1)*sqrt(sum(B.*B)));
C = C./(ones(dimV(1),1)*sqrt(sum(C.*C)));
D = D./(ones(dimV(1),1)*sqrt(sum(D.*D)));

im = mean(A.*B.*C.*D); %dot product
im = reshape(im',[dim(1) dim(2)]);


dim = size(im);
dimI = size(T); dimI = dimI(1:2); %Get back original dim
im = interp1(1:dim(1),im,linspace(1,dim(1),dimI(1)));
im = interp1(1:dim(2),im',linspace(1,dim(2),dimI(2)))';
% 
% im = im/std(im(:));
% im = im-median(im(:)); 




