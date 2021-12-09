cellS_RemoveSingValues


SVDremovalFlag = 1;
N_PC= 5;

if SVDremovalFlag
    
    kdum = cellS.cellMat{i};
    idGoodFrames = ones(1,size(kdum,2));
    idGoodFrames(find(isnan(kdum(1,:)))) = 0;
    kern_NaNdump = kdum;
    kern_NaNdump = kern_NaNdump(:,find(idGoodFrames));
    
    kernPCdump = TensorRemoveSingularValues(kern_NaNdump',N_PC)';
    
    kdum(:,find(idGoodFrames)) = kernPCdump; %Replace with the PCA subtracted matrix. The others are NaN;
    
end

figure,imagesc(kdum - nanmean(kdum,2)*ones(1,size(kdum,2)))