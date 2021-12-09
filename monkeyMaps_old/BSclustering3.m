function BSclustering2(oripref,sfpref,tcori,tcsf,dori,dsf,doriEuc,dsfEuc,dist,animID,animID_d)

%bootstrap analysis for ori/sf clustering

%trims the stuff I don't use anyway (the correlation analysis).  Also, it
%includes the expt ID to do the local clustering

%%

doriAll = []; dsfAll = []; doriEucAll = []; dsfEucAll = []; distAll = []; animIDdAll = [];
for i = 1:length(dori)
    
    doriAll = [doriAll; dori{i}(:)];
    dsfAll = [dsfAll; dsf{i}(:)];
    doriEucAll = [doriEucAll; doriEuc{i}(:)];
    dsfEucAll = [dsfEucAll; dsfEuc{i}(:)];
    distAll = [distAll; dist{i}(:)];
    animIDdAll = [animIDdAll; animID_d{i}(:)];
    
end

prcdom = 0:20:100;
PW.Ddom = [];
for i = 1:length(prcdom)
    PW.Ddom = [PW.Ddom prctile(distAll,prcdom(i))];
end

for i = 1:length(PW.Ddom)-1    
   labs{i} = [num2str(round(PW.Ddom(i))) ' to '  num2str(round(PW.Ddom(i+1)))];
end

clear dori dsf doripair dsfpair dtcoripair dtcsfpair doriEuc dsfEuc animID_d
for i = 1:length(PW.Ddom)-1
    
    id = find(distAll>PW.Ddom(i) & distAll<=PW.Ddom(i+1));
    distdom(i) = mean(distAll(id));
    dori{i} = doriAll(id);
    dsf{i} = dsfAll(id);
    doriEuc{i} = doriEucAll(id);
    dsfEuc{i} = dsfEucAll(id);
    animID_d{i} = animIDdAll(id);
    
end

%%

Nsim = 10;

for i = 1:length(PW.Ddom)
    
    
    
    
end






