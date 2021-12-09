function [bestWindow tempIm] = findTempWindow(T,nT)

%% Get the window that is most correlated with a subset of the experiment.
%This is a modification of findTempTrial. Instead of looping through
%trials, it loops through a user defined window at nonoverlapping
%intervals. The window as set as an input. 

%T length of the window in seconds
%nT number of windows to search, starting at the beginning of the
%experiment.

global ACQinfo Analyzer G_handles

datadir = get(G_handles.datadir,'string');

display('finding best windoow for template...');

fps = 1000/(ACQinfo.msPerLine*ACQinfo.linesPerFrame);
Nframes = round(T*fps);

for t = 1:nT
    
    Nstart = (t-1)*Nframes + 1;
    imgs = single(sbxread(datadir,Nstart,Nframes));
    CHs{1} = squeeze(imgs);
    CHs{1} = CHs{1}(:,ACQinfo.unblanked,:);
     
    dim = size(CHs{1});
    
    if t == 1
        matim = zeros(dim(1)*dim(2),nT);
    end
    
    im = mean(CHs{1},3);
    %figure, imagesc(im)
    matim(:,t) = (im(:)-mean(im(:))/std(im(:)));
    
end


rmat = matim'*matim;

[dum bestWindow] = max(sum(rmat));

%tempIm = matim(:,bestWindow);
%tempIm = reshape(tempIm);

Nstart = (bestWindow-1)*Nframes + 1;
imgs = single(sbxread(datadir,Nstart,Nframes));
CHs{1} = squeeze(imgs);
CHs{1} = CHs{1}(:,ACQinfo.unblanked,:);
tempIm = mean(CHs{1},3);



