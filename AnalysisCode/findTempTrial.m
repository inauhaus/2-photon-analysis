function bestTrial = findTempTrial

%% Get the trial average that is most correlated with the rest of the experiment.

display('finding best trial for template...');

for t = 1:getnotrials
    
    CH = GetTrialData([1 0 0 0],t);
    dim = size(CH{1});
    
    if t == 1
        matim = zeros(dim(1)*dim(2),getnotrials);
    end
    
    im = mean(CH{1},3);
    matim(:,t) = (im(:)-mean(im(:))/std(im(:)));
    
end



rmat = matim'*matim;

[dum bestTrial] = max(sum(rmat));

bestTrial
