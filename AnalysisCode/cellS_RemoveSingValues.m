function cellMat = cellS_RemoveSingValues(cellMat,N_PC)


SVDremovalFlag = 1;


try
    if SVDremovalFlag
        
        kdum = cellMat;
        idGoodFrames = ones(1,size(kdum,2));
        idGoodFrames(find(isnan(kdum(1,:)))) = 0;
        kern_NaNdump = kdum;
        kern_NaNdump = kern_NaNdump(:,find(idGoodFrames));
        
        kernPCdump = TensorRemoveSingularValues(kern_NaNdump',N_PC)';
        
        kdum(:,find(idGoodFrames)) = kernPCdump; %Replace with the PCA subtracted matrix. The others are NaN;
        
    end
catch
    'hi'
    kdum = cellMat;
end

cellMat = kdum;

%figure,imagesc(kdum - nanmean(kdum,2)*ones(1,size(kdum,2)))