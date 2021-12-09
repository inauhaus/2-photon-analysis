function im = localXCorr4(T,W)

%input is tensor(x,y,T), output is image(x,y)

global ACQinfo
%%

T = TensorZscore_time(T);
Gsig = 1;
%Get movies that are downsampled in space
r = floor(W/2);
rdom = -r:r;
for i = 1:W
    for j = 1:W      
        ishift = rdom(i);
        jshift = rdom(j);
        TD{i,j} = circshift(T,[ishift jshift 0]);              
        
        D = sqrt(ishift.^2 + jshift.^2);
        G = exp(-D^2/(2*Gsig^2));
        TD{i,j} = TD{i,j}*G;
    end
end

%Combine
im = ones(size(T));
for i = 1:W
    for j = 1:W      
        im = im.*TD{i,j};        
    end
end

im = mean(im,3);

figure, imagesc(im,[0 .5]), colormap gray

% 
% im = im/std(im(:));
% im = im-median(im(:)); 




