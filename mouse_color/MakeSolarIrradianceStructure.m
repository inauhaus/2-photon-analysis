%loaded from data

clear Rural
k = 1;
for i = 72:3192
    vname = ['VarName' num2str(i)];
    dvec = eval(vname);
    Rural.SolarElev(k) = dvec(3);
    Rural.LunarElev(k) = dvec(4);
    Rural.MoonFraction(k) = dvec(5);
    Rural.PSD(:,k) = dvec(6:566);
    
    k = k+1;
end

Rural.nm = 280:840;

%save solarIrradianceData Rural

%%

clear PSDmu elevmu
elevationBinEdges = [-30:30:90];
for i = 1:length(elevationBinEdges)-1
    
    id = find(Rural.SolarElev > elevationBinEdges(i) & Rural.SolarElev < elevationBinEdges(i+1));
    
    PSDmu(:,i) = mean(Rural.PSD(:,id),2);
    PSDmu_norm(:,i) = PSDmu(:,i)/max(PSDmu(:,i));
    
    elevmu(i) = mean(Rural.SolarElev(id));
    
end

figure,semilogy(Rural.nm,PSDmu)

figure,plot(elevmu,max(PSDmu))