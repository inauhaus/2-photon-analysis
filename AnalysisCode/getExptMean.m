function [CHs] = getExptMean(chvec,Ntrials)

N = 0;

CHs{1} = 0; CHs{2} = 0;
for trial = 1:Ntrials
    
    CHdum = GetTrialData(chvec,trial);
    
    N = N + length(CHdum{1}(1,1,2:end-2));
    
    for i = 1:sum(chvec)
        CHs{i} = sum(CHdum{i}(:,:,2:end-2),3) + CHs{i};
    end
    
%     if chvec(1)
%         CHs{1} = sum(CHdum{1}(:,:,2:end-2),3) + CHs{1};
%     end
%     if chvec(2)
%         CHs{2} = sum(CHdum{2}(:,:,2:end-2),3) + CHs{2};
%     end
end

for i = 1:sum(chvec)
    CHs{i} = CHs{i}/N;
end
    