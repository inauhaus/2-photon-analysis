function [H p] = getMapEntropy(orimap,sfmap,mask,shuffbit)

oriedges = 0:2:180;
sfedges = logspace(log10(.5),log10(8),5);

sfedges = linspace(-1,1,3); %Use for ocular dominance


sfpref = sfmap(find(mask));
oripref = orimap(find(mask));

if shuffbit
    [dum id] = sort(rand(length(sfpref),1));
    sfpref = sfpref(id);
end


edges{1} = oriedges;
edges{2} = sfedges;
p = hist3([oripref sfpref],'Edges',edges);
p = p(1:end-1,1:end-1);

% p = zeros(length(oriedges)-1,length(sfedges)-1);
% for i = 1:length(oriedges)-1
%     id = find(oripref>oriedges(i) & oripref<oriedges(i+1));
%     dum = histc(sfpref(id),sfedges);
%      p(i,:) = dum(1:end-1);
% end

p = p/sum(p(:));

pzero = find(p == 0);

H = -log2(p).*p;
H(pzero) = 0;
H = sum(H(:));

H = 1/std(p(:));
