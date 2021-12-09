
clear opref prefAxisMF
for i = 1:length(posplanex)
   prefAxisMF(i) = getMagAxis(posplanex{i},posplaney{i}); 
   opref(i) = oprefAll{i}; 
end

axisdiff = 90-abs(abs(opref-prefAxisMF)-90);
figure,hist(axisdiff), xlim([0 90])