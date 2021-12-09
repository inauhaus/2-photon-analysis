function plotRetOriStuff(oripref,prefAxisMF,JacvecSF,SFvec,JacvecOD,ocdomvec)

oriprefAll = [];
prefAxisMFAll = [];
JacvecSFAll = [];
JacvecODAll = [];
SFvecAll = [];
ocdomvecAll = [];
for i = 1:length(oripref)    
   oriprefAll = [oriprefAll; oripref{i}(:)]; 
   prefAxisMFAll = [prefAxisMFAll; prefAxisMF{i}(:)];     
   
   JacvecSFAll = [JacvecSFAll; JacvecSF{i}(:)];  
   SFvecAll = [SFvecAll; SFvec{i}(:)];  
%    JacvecSFAll = [JacvecSFAll; JacvecSF{i}(:)];  
%    SFvecAll = [SFvecAll; SFvec{i}(:)];  
   
   JacvecODAll = [JacvecODAll; JacvecOD{i}(:)/median(JacvecOD{i}(:))];  
   ocdomvecAll = [ocdomvecAll; ocdomvec{i}(:)];  
end

figure,
subplot(1,2,1)
scatter(oriprefAll(:),prefAxisMFAll(:),'.k')
xlabel('Preferred Orientation')
ylabel('Short axis of visual space representation')
axis square
hold on
plot([0 180],[0 180])

subplot(1,2,2)
intsect = 90-abs(abs(oriprefAll(:)-prefAxisMFAll(:))-90);
%intsect = abs(oridum(:)-prefAxisMF(:));
hist(intsect,8), xlim([0 90])

[r, p] = circCorr(oriprefAll, prefAxisMFAll,180,1)

%%

figure,
scatter(log2(JacvecSFAll(:)),SFvecAll(:),'.k')
[r p] = corrcoef(JacvecSFAll(:),SFvecAll(:));
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
ylabel('SF preference (octaves)'), xlabel('Magnification factor (log)')

%%

figure,
scatter(JacvecODAll(:),ocdomvecAll(:),'.k')
[r p] = corrcoef(JacvecODAll(:),ocdomvecAll(:));
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
ylabel('monocularity'), xlabel('Magnification factor')