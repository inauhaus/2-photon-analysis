function [U p] = getMapEntropy2(orimap,sfmap,odmap,mask,varargin)

%2 takes three maps

oriedges = 0:5:180;
sfedges = logspace(log10(.5),log10(8),30);
odedges = linspace(-1,1,20);

sfpref = sfmap(find(mask));
oripref = orimap(find(mask));
odpref = odmap(find(mask));

% if shuffbit
%     [dum id] = sort(rand(length(sfpref),1));
%     sfpref = sfpref(id);
%     
%     [dum id] = sort(rand(length(odpref),1));
%     odpref = odpref(id);
% end

edges{1} = oriedges;
edges{2} = sfedges;
p = zeros(length(oriedges)-1,length(sfedges)-1,length(odedges)-1);
for i = 1:length(odedges)-1
    id = find(odpref>odedges(i) & odpref<odedges(i+1));
    pdum = hist3([oripref(id) sfpref(id)],'Edges',edges);
    p(:,:,i) = pdum(1:end-1,1:end-1);
end

%p = p(:,3,2);

p = p/sum(p(:));

m1 = squeeze(sum(sum(p,2),3));
m2 = squeeze(sum(sum(p,1),3));
m3 = squeeze(sum(sum(p,1),2));
for i = 1:size(p,3);
    psep(:,:,i) = m1*m2*m3(i);
end
% psep = psep/norm(psep(:));
% p = p/norm(p(:));

%Get marginal
if ~isempty(varargin)
    p = squeeze(sum(p,varargin{1}));
    
    %[u s v] = svd(p);
    %psep = u(:,1)*v(:,1)'*s(1,1);
    
    
    %psep = squeeze(sum(psep,varargin{1}));
end

% pzero = find(p == 0);
% 
% H = -log2(p).*p;
% H(pzero) = 0;
% H = sum(H(:));
% 
% pzero = find(psep == 0);
% 
% Hsep = -log2(psep).*psep;
% Hsep(pzero) = 0;
% Hsep = sum(Hsep(:));

varacc = (var(p(:))-var(p(:)-psep(:)))/var(p(:));
% r = corrcoef(p(:),psep(:)); r = r(1,2);

U = varacc;

%U = s(1,1)/sum(s(:));

%U = r;

%U = H/Hsep;