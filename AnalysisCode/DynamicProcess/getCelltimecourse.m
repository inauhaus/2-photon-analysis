function [cellMat] = getCelltimecourse(Tens,idcell,CoM,cellWidth)

global maskS

Ncell = length(idcell);
Nt = length(Tens(1,1,:)); 

cellMat = zeros(Ncell,Nt);

%%%%
%%%%
% Tdim = size(Tens);
% CoM = round(CoM);

% dum = cellWidth(2:end,:);
% sig = mean(dum(:));
% sig = 1.5*sig;

%smoother = fspecial('gaussian',[Tdim(1) Tdim(2)],sig);
%smoother = abs(fft2(smoother));



if isfield(maskS,'maskTens')
    
    Ncell = size(maskS.maskTens,3);
    
    %     for p = 1:Ncell
    %         dum = maskS.maskTens(:,:,p);
    %         id{p} = find(dum);
    %     end
    
    %First vectorize the tensors
    dim = size(Tens);
    Tens = reshape(Tens,[dim(1)*dim(2) dim(3)])';
    
    dim = size(maskS.maskTens);
    maskTens = reshape(maskS.maskTens,[dim(1)*dim(2) dim(3)]);
    cellMat = (Tens*maskTens)';

    
    cellMatnew = zeros([size(cellMat,1)+1 size(cellMat,2)]);
    idim = bwlabel(sign(maskS.bwCell{1}));
    cellID = unique(idim);
    %This leaves the initial entry (i.e. neuropil) of signew = 0.
    for i = 2:length(cellID)
        idx = find(idim == cellID(i));
        cid = maskS.bwCellID{1}(idx(1));
        cellMatnew(i,:) = cellMat(cid,:);
    end
    cellMat = cellMatnew;
else
    
    %%Should use this cuz it will be faster
%     dim = size(Tens(:,:,1));
%     for p = 1:Ncell
%         idlocs = idcell{p}(:)*ones(1,Nt) + ones(length(idcell{p}),1)*((0:Nt-1)*dim(1)*dim(2));
%         cellMat(p,:) = mean(Tens(idlocs))';
%     end

%I think the above was just dumb. This is better.
    dim = size(Tens);
    Tens = reshape(Tens,[dim(1)*dim(2) dim(3)])';

    for p = 1:Ncell
        cellMat(p,:) = mean(Tens(:,idcell{p}),2)';
    end

end




%%%%
%%%%

% fvec = numel(Tens(:,:,1))*(0:Nt-1);
% for p = 1:Ncell
%     tcell = zeros(length(idcell{p}),Nt);
%     for z = 1:length(idcell{p})
%         tcell(z,:) = Tens(fvec+idcell{p}(z));      
%     end    
%     cellMat(p,:) = mean(tcell);
% end

% 
% for p = 2:Ncell
%     idx = round(CoM(p,2)-3*cellWidth(p,2)):round(CoM(p,2)+3*cellWidth(p,2));
%     idy = round(CoM(p,1)-3*cellWidth(p,1)):round(CoM(p,1)+3*cellWidth(p,1));
%     
%     sig = mean(cellWidth(p,:));
%     mask = fspecial('gaussian',[length(idy) length(idx)],sig);
%     
%     tcell = zeros(length(idy),length(idx),Nt);
%     for z = 1:length(idx)
%         for k = 1:length(idy)
%             tcell(k,z,:) = Tens(k,z,1:Nt)*mask(k,z);      
%         end
%     end    
%     cellMat(p,:) = squeeze(mean(mean(tcell,1),2));
% end

