function getGaussMasks

global maskS

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
Ncell = length(celldom);

CoMy = zeros(1,Ncell);
CoMx = zeros(1,Ncell);
for p = 1:Ncell
    [idy0{p} idx0{p}] = find(masklabel == celldom(p));
    CoMy(p) = mean(idy0{p});
    CoMx(p) = mean(idx0{p});
end

%This loop takes the cell location from the previous loops to create the
%actual time course.
Wcell = 5;
dim = size(maskS.bwCell{1});
Gmask = zeros(dim);
for p = 1:Ncell
p
    %Cut out a piece around the cell, approximately centered on the cell, the dxm shift will take care of the rest:
    yran = (round(CoMy(p))-floor(Wcell/2)):(round(CoMy(p))+floor(Wcell/2));
    xran = (round(CoMx(p))-floor(Wcell/2)):(round(CoMx(p))+floor(Wcell/2));
    
    id = find(xran<1 | xran > dim(2));
    xran(id) = [];
    id = find(yran<1 | yran > dim(1));
    yran(id) = [];
    
    impiece = maskS.im{1}(yran,xran);
    
    [param ffit varaccount] = Gaussfit2Drot_const(impiece.^2,0);
    ffit = ffit-min(ffit(:));
    ffit = ffit/max(ffit(:));
    
    Gmask(yran,xran) = Gmask(yran,xran) + ffit; 

end

figure,imagesc(Gmask), colormap gray
