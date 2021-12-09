function plotOrivSF(tcorisfraw,opref,sfpref,oridom,sfdom,varargin)

%plotOrivSF(TC.tcorisfraw{1},TC.opref{1},TC.sfpref{1},DM.oridom,DM.sfdom,[.5 1 2 8])

%length(find(sfpref>sfBands(1) & sfpref<sfBands(2)))
%length(find(sfpref>sfBands(2) & sfpref<sfBands(3)))
%length(find(sfpref>sfBands(3) & sfpref<sfBands(4)))

if isempty(varargin)
    sfBands(1) = prctile(sfpref,0)
    sfBands(2) = prctile(sfpref,33)
    sfBands(3) = prctile(sfpref,66)
    sfBands(4) = prctile(sfpref,100)
else
    sfBands = varargin{1};
end

oridom = linspace(0,180,length(tcorisfraw(1,:,1))+1);
oridom = oridom(1:end-1);
clear kernshift
for i = 1:size(tcorisfraw,1)
    
    kern = squeeze(tcorisfraw(i,:,:));
    
    [dum id] = min(abs(oridom-opref(i)));
    
    kernshift(i,:,:) = circshift(kern,[5-id 0]);
    
end

clear mukern oritcAll sftcAll
for i = 1:length(sfBands)-1
    
    id = find(sfpref>sfBands(i) & sfpref<=sfBands(i+1));
    
    mukern{i} = squeeze(nanmean(kernshift(id,:,:),1));
    
    oritcAll(:,i) = nanmean(mukern{i}');
    sftcAll(:,i) = nanmean(mukern{i});
    
    oritcAll(:,i) = oritcAll(:,i)-min(oritcAll(:,i));
    oritcAll(:,i) = oritcAll(:,i)/max(oritcAll(:,i));
end

figure
for i = 1:length(mukern)
    
   subplot(length(mukern)+1,1,i+1)
   imagesc(oridom,log2(sfdom),mukern{i}')
   set(gca,'YTick',round(log2(sfdom)*10)/10,'XTick',oridom(1:2:end),'Tickdir','out')
   colormap jet
   axis xy 
   
end

subplot(length(mukern)+1,1,1)
plot(oridom,oritcAll,'o-')
set(gca,'XTick',oridom(1:2:end))
%xlim([-20 180])

figure
subplot(1,2,1)
plot(oridom,oritcAll,'.-')
subplot(1,2,2)
plot(log2(sfdom),sftcAll,'.-')


