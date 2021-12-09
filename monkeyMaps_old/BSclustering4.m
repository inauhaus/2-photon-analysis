function BSclustering4(oripref,sfpref,phase,tcori,tcsf,tccolor,doriAll,dsfAll,dphaseAll,doriEucAll,dsfEucAll,dcolorEucAll,distAll,animID,animIDdAll)

%4 is the same as two, but with the addition of color

%bootstrap analysis for ori/sf clustering

%trims the stuff I don't use anyway (the correlation analysis).  Also, it
%includes the expt ID to do the local clustering

%phase is treated differently throughout the code

%%

%Get rid of NaN from phase pairs
id = find(~isnan(dphaseAll));
dphaseAll = dphaseAll(id);
animIDdAll_Phase = animIDdAll(id); distAll_Phase = distAll(id);
%id = find(~isnan(phase));
phase = phase(id,:);  %This is actually from all the pairs, so its huge.  i.e. can't assign a single phase for each cell

%Get rid of NaN from ori, sf, and color pairs
id = find(isnan(doriAll.*dsfAll));
doriAll(id) = []; dsfAll(id) = []; doriEucAll(id) = [];  dsfEucAll(id) = [];  distAll(id) = []; animIDdAll(id) = [];
dcolorEucAll(id) = [];

%Get rid of NaN from ori and sf cells
id = find(~isnan(oripref.*sfpref));
oripref = oripref(id);
sfpref = sfpref(id);
tcori = tcori(id,:);
tcsf = tcsf(id,:);
tccolor = tccolor(id,:);
animID = animID(id);


prcdom = 0:20:100;
Dco = [];
for i = 1:length(prcdom)
    Dco = [Dco prctile(distAll,prcdom(i))];  %distanace domain cutoff points
end

for i = 1:length(Dco)-1    
   labs{i} = [num2str(round(Dco(i))) ' to '  num2str(round(Dco(i+1)))];
end

clear distdom dori dsf doripair dsfpair dtcoripair dtcsfpair doriEuc dsfEuc dcolorEuc animID_d
for i = 1:length(Dco)-1
    
    id = find(distAll>Dco(i) & distAll<=Dco(i+1));
    distdom(i) = mean(distAll(id));
    dori{i} = doriAll(id);
    dsf{i} = dsfAll(id);    
    doriEuc{i} = doriEucAll(id);
    dsfEuc{i} = dsfEucAll(id);
    dcolorEuc{i} = dcolorEucAll(id);
    animID_d{i} = animIDdAll(id);
    
    id = find(distAll_Phase>Dco(i) & distAll_Phase<=Dco(i+1));
    distdom_Phase(i) = mean(distAll_Phase(id));
    dphase{i} = dphaseAll(id);
    dphaseEuc{i} = (1-cos(dphaseAll(id)*pi/180))/2+.5;
    animID_d_Phase{i} = animIDdAll_Phase(id);
    
end

%%

Ncell = length(oripref);

NcellPhase = length(phase(:)); %this is a huge number 2 x number of pairs

Nsim = 10;

for d = 1:length(distdom)
    d
    Ngrabs = length(dori{d});
    NgrabsPhase = length(dphase{d});  %this will be smaller
    
    doriBS{d} = zeros(Ngrabs,Nsim);
    dsfBS{d} = zeros(Ngrabs,Nsim);

    doriEucBS{d} = zeros(Ngrabs,Nsim);
    dsfEucBS{d} = zeros(Ngrabs,Nsim);   
    dcolorEucBS{d} = zeros(Ngrabs,Nsim); 
    
    dphaseBS{d} = zeros(NgrabsPhase,Nsim);
    
    for i = 1:Nsim
        for j = 1:Ngrabs

            rpick1 = round(1+rand(1)*(Ncell-1));  %randomly select a cell         
            rpick2 = rpick1;
            while rpick2 == rpick1
                rpick2 = round(1+rand(1)*(Ncell-1));  %randomly select a cell
            end

            doridum = abs(oridiff(oripref(rpick1)*pi/180,oripref(rpick2)*pi/180)*180/pi); %degrees
            dsfdum = abs(log2(sfpref(rpick1)/sfpref(rpick2)));  %octaves
            
            doriBS{d}(j,i) = doridum;
            dsfBS{d}(j,i) = dsfdum;
            
            v1 = tcori(rpick1,:); v2 = tcori(rpick2,:);            
            doriEucdum = corrcoef(v1,v2); doriEucdum = -doriEucdum(1,2); doriEucdum = (doriEucdum+1)/2;
            doriEucBS{d}(j,i) = doriEucdum;
            
            v1 = tcsf(rpick1,:); v2 = tcsf(rpick2,:);            
            dsfEucdum = corrcoef(v1,v2); dsfEucdum = -dsfEucdum(1,2); dsfEucdum = (dsfEucdum+1)/2;
            dsfEucBS{d}(j,i) = dsfEucdum;           
            
            v1 = tccolor(rpick1,:); v2 = tccolor(rpick2,:);            
            dcolorEucdum = corrcoef(v1,v2); dcolorEucdum = -dcolorEucdum(1,2); dcolorEucdum = (dcolorEucdum+1)/2;
            dcolorEucBS{d}(j,i) = dcolorEucdum;  

        end
        
        for j = 1:NgrabsPhase
           
            rpick1 = round(1+rand(1)*(NcellPhase-1));  %randomly select a cell         
            rpick2 = rpick1;
            while rpick2 == rpick1
                rpick2 = round(1+rand(1)*(NcellPhase-1));  %randomly select a cell
            end
            
            %dphasedum = abs(angle(exp(1i*pi/180*(rand(1)*360-rand(1)*360)))*180/pi); %degrees  (uniform dist)

            dphasedum = abs(angle(exp(1i*pi/180*(phase(rpick1)-phase(rpick2))))*180/pi); %degrees          
            
            dphaseBS{d}(j,i) = dphasedum;
            
            dphaseEucBS{d}(j,i) = (1-cos(dphasedum*pi/180))/2+0.5;  %normalized from 0 to 1
            
        end

    end
     
end

%%

for d = 1:length(distdom)
    Ngrabs = length(dori{d});

    %Mean
    medoriBS_trial{d} = mean(doriBS{d}); %mean of each trial
    medsfBS_trial{d} = mean(dsfBS{d});
    medphaseBS_trial{d} = mean(dphaseBS{d});
    medoriEucBS_trial{d} = mean(doriEucBS{d}); %mean of each trial
    medsfEucBS_trial{d} = mean(dsfEucBS{d});
    medcolorEucBS_trial{d} = mean(dcolorEucBS{d});
    medphaseEucBS_trial{d} = mean(dphaseEucBS{d});

    medoriBS(d) = mean(medoriBS_trial{d});  %total mean for each distance
    medsfBS(d) = mean(medsfBS_trial{d});
    medphaseBS(d) = mean(medphaseBS_trial{d});
    medoriEucBS(d) = mean(medoriEucBS_trial{d});  %total mean for each distance
    medsfEucBS(d) = mean(medsfEucBS_trial{d});
    medcolorEucBS(d) = mean(medcolorEucBS_trial{d});
    medphaseEucBS(d) = mean(medphaseEucBS_trial{d});

    %Stand Dev
    stdoriBS_trial{d} = std(doriBS{d})/sqrt(Ngrabs); %sig of each trial
    stdsfBS_trial{d} = std(dsfBS{d})/sqrt(Ngrabs);
    stdphaseBS_trial{d} = std(dphaseBS{d})/sqrt(NgrabsPhase);
    stdoriEucBS_trial{d} = std(doriEucBS{d})/sqrt(Ngrabs); %sig of each trial
    stdsfEucBS_trial{d} = std(dsfEucBS{d})/sqrt(Ngrabs);
    stdcolorEucBS_trial{d} = std(dcolorEucBS{d})/sqrt(Ngrabs);
    stdphaseEucBS_trial{d} = std(dphaseEucBS{d})/sqrt(NgrabsPhase);

    stdoriBS(d) = mean(stdoriBS_trial{d});  %mean sig across trials, for each distance
    stdsfBS(d) = mean(stdsfBS_trial{d});
    stdphaseBS(d) = mean(stdphaseBS_trial{d});
    stdoriEucBS(d) = mean(stdoriEucBS_trial{d});  %mean sig across trials, for each distance
    stdsfEucBS(d) = mean(stdsfEucBS_trial{d});
    stdcolorEucBS(d) = mean(stdcolorEucBS_trial{d});
    stdphaseEucBS(d) = mean(stdphaseEucBS_trial{d});
end


%%

for d = 1:length(dori)
    
    %preference difference
    medratio_ori(d) = medoriBS(d)/mean(dori{d});
    medratio_sf(d) = medsfBS(d)/mean(dsf{d});
    medratio_phase(d) = medphaseBS(d)/mean(dphase{d});
    
    stdratio_ori(d) = stdoriBS(d)/mean(dori{d});
    stdratio_sf(d) = stdsfBS(d)/mean(dsf{d});
    stdratio_phase(d) = stdphaseBS(d)/mean(dphase{d});
    
    id = find(medoriBS_trial{d} > mean(dori{d}));
    percOri(d) = length(id)/Nsim;
    id = find(medsfBS_trial{d} > mean(dsf{d}));
    percSf(d) = length(id)/Nsim;
    id = find(medphaseBS_trial{d} > mean(dphase{d}));
    percPhase(d) = length(id)/Nsim;
    
    %Euclidean difference
    medratio_oriEuc(d) = medoriEucBS(d)/mean(doriEuc{d});
    medratio_sfEuc(d) = medsfEucBS(d)/mean(dsfEuc{d});
    medratio_colorEuc(d) = medcolorEucBS(d)/mean(dcolorEuc{d});
    medratio_phaseEuc(d) = medphaseEucBS(d)/mean(dphaseEuc{d});
    
    stdratio_oriEuc(d) = stdoriEucBS(d)/mean(doriEuc{d});
    stdratio_sfEuc(d) = stdsfEucBS(d)/mean(dsfEuc{d});
    stdratio_colorEuc(d) = stdcolorEucBS(d)/mean(dcolorEuc{d});
    stdratio_phaseEuc(d) = stdphaseEucBS(d)/mean(dphaseEuc{d});
    
    id = find(medoriEucBS_trial{d} > mean(doriEuc{d}));
    percOriEuc(d) = length(id)/Nsim;
    id = find(medsfEucBS_trial{d} > mean(dsfEuc{d}));
    percSfEuc(d) = length(id)/Nsim;
    id = find(medcolorEucBS_trial{d} > mean(dcolorEuc{d}));
    percColorEuc(d) = length(id)/Nsim;
    id = find(medphaseEucBS_trial{d} > mean(dphaseEuc{d}));
    percPhaseEuc(d) = length(id)/Nsim;
end


percPhase
percOri
percSf

%%
%Plot diff preference
figure, errorbar([distdom' distdom' distdom'],([medratio_ori' medratio_sf'  medratio_phase']),([stdratio_ori' stdratio_sf'  stdratio_phase']))
hold on, plot([distdom(1) distdom(end)],[1 1],'--k')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)    

ylabel('mean(ShuffledDist)/mean(ActualDist)')
xlabel('Distance')
title('Difference in preferred')
legend('ori','spfreq','phase')

hold on
for i = 1:length(dori)
   if percOri(i)<.01 | percOri(i)>.99 
       plot(distdom(i)+2,(medratio_ori(i))+.1,'*b')
   end
   if percSf(i)<.01  | percSf(i)>.99 
       plot(distdom(i)+2,(medratio_sf(i)+.1),'*g')
   end
end
ylim([.5 medratio_ori(1)+.5])

%Plot Euclidian (I put the dphase in here too (same one above))
figure, errorbar([distdom' distdom' distdom' distdom'],([medratio_oriEuc' medratio_sfEuc' medratio_colorEuc' medratio_phaseEuc']),([stdratio_oriEuc' stdratio_sfEuc' stdratio_colorEuc' stdratio_phaseEuc']))
hold on, plot([distdom(1) distdom(end)],[1 1],'--k')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)    

ylabel('mean(ShuffledDist)/mean(ActualDist)')
xlabel('Distance')
title('Euclidean')
legend('ori','spfreq','color','phase')

hold on
for i = 1:length(dori)
   if percOriEuc(i)<.01 | percOriEuc(i)>.99 
       plot(distdom(i)+2,(medratio_oriEuc(i))+.1,'*b')
   end
   if percSfEuc(i)<.01  | percSfEuc(i)>.99 
       plot(distdom(i)+2,(medratio_sfEuc(i))+.1,'*g')
   end
end
ylim([.5 medratio_oriEuc(1)+.5])
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now do bootstrap by shuffling individual ROIs

doriBS = []; dsfBS = [];  dphaseBS = []; 

adom = unique(animID);

animIDdAll_Phasedum = [animIDdAll_Phase animIDdAll_Phase];  %this is just used here to access both sides of each pair
for R = 1:length(adom)  %parse the values according to the ROI

    id = find(adom(R) == animID);
    oripref_R{R} = oripref(id);
    sfpref_R{R} = sfpref(id);    
    tcori_R{R} = tcori(id,:);
    tcsf_R{R} = tcsf(id,:);
    tccolor_R{R} = tccolor(id,:);
    
    id = find(adom(R) == animIDdAll_Phasedum);
    phase_R{R} = phase(id);
end


for R = 1:length(adom)  %parse the pairings according to the ROI

    id = find(adom(R) == animIDdAll ); %id pairs in this ROI
    dori_R{R} = doriAll(id);
    dsf_R{R} = dsfAll(id);
    
    id = find(adom(R) == animIDdAll_Phase);
    dphase_R{R} = dphaseAll(id);
end

Nsim = 2;

for d = 1:length(distdom)  %loop each distance
    d
    for R = 1:length(adom)  %loop each ROI
        
        Ncell = length(oripref_R{R}); %number of cells in this ROI
        
        Npairs = length(dori_R{R}); %number of pairs in this ROI

        NpairsPhase = length(dphase_R{R}); %number of pairs for which I could compute phase difference, in this ROI
    
        ida = find(animID_d{d} == adom(R)); %id the pairs for this ROI and distance
        Ngrabs = length(ida);  %number of pairs in ROI at this distance
        
        ida = find(animID_d_Phase{d} == adom(R)); %id the (available phase) pairs for this ROI and distance
        NgrabsPhase = length(ida);  %number of available phase pairs in ROI at this distance
        
        doriBS{d,R} = zeros(Ngrabs,Nsim);
        dsfBS{d,R} = zeros(Ngrabs,Nsim);
        doriEucBS{d,R} = zeros(Ngrabs,Nsim);
        dsfEucBS{d,R} = zeros(Ngrabs,Nsim);
        dcolorEucBS{d,R} = zeros(Ngrabs,Nsim);
        
        dphaseBS{d,R} = zeros(NgrabsPhase,Nsim);
        
        for i = 1:Nsim
            
            for j = 1:Ngrabs


                %pick two cells and compute the difference
                rpick1 = round(1+rand(1)*(Ncell-1));  %randomly select a cell in ROI
                rpick2 = rpick1;
                while rpick2 == rpick1
                    rpick2 = round(1+rand(1)*(Ncell-1));  %randomly select a cell in ROI
                end
% 
%                 doridum = abs(oridiff(oripref_R{R}(rpick1)*pi/180,oripref_R{R}(rpick2)*pi/180)*180/pi); %degrees
%                 dsfdum = abs(log2(sfpref_R{R}(rpick1)/sfpref_R{R}(rpick2)));  %octaves
%                 dphasedum = abs(angle(exp(1i*pi/180*(phase_R{R}(rpick1)-phase_R{R}(rpick2))))*180/pi); %degrees
                
                %Pick a pair where difference was already computed (necessary for phase)
                rpick = round(1+rand(1)*(Npairs-1));                  
                
                doridum = dori_R{R}(rpick);
                dsfdum = dsf_R{R}(rpick);                                 

                doriBS{d,R}(j,i) = (doridum/mean(dori{d}));
                dsfBS{d,R}(j,i) = (dsfdum/mean(dsf{d}));                
                
                %doriBS{d,R}(j,i) = (doridum/mean(dori{d}(ida))); %normalize by the mean within this distance
                %dsfBS{d,R}(j,i) = (dsfdum/mean(dsf{d}(ida)));

                v1 = tcori_R{R}(rpick1,:); v2 = tcori_R{R}(rpick2,:);
                doriEucdum = corrcoef(v1,v2); doriEucdum = -doriEucdum(1,2); doriEucdum = (doriEucdum+1)/2;

                v1 = tcsf_R{R}(rpick1,:); v2 = tcsf_R{R}(rpick2,:);
                dsfEucdum = corrcoef(v1,v2); dsfEucdum = -dsfEucdum(1,2); dsfEucdum = (dsfEucdum+1)/2;
                
                v1 = tccolor_R{R}(rpick1,:); v2 = tccolor_R{R}(rpick2,:);
                dcolorEucdum = corrcoef(v1,v2); dcolorEucdum = -dcolorEucdum(1,2); dcolorEucdum = (dcolorEucdum+1)/2;
                
                doriEucBS{d,R}(j,i) = (doriEucdum/mean(doriEuc{d}));
                dsfEucBS{d,R}(j,i) = (dsfEucdum/mean(dsfEuc{d}));
                dcolorEucBS{d,R}(j,i) = (dcolorEucdum/mean(dcolorEuc{d}));
                
                %doriEucBS{d,R}(j,i) = (doriEucdum/mean(doriEuc{d}(ida))); %normalize by the mean within this distance
                %dsfEucBS{d,R}(j,i) = (dsfEucdum/mean(dsfEuc{d}(ida)));

            end
            
            
            for j = 1:NgrabsPhase
                
                rpickPhase = round(1+rand(1)*(NpairsPhase-1));
                dphasedum = dphase_R{R}(rpickPhase);
                dphaseBS{d,R}(j,i) = (dphasedum/mean(dphase{d}));
                
            end

        end
    end
    
end

%%
ND = length(distdom); 
doriBS2 = cell(1,ND); dsfBS2 = cell(1,ND); dphaseBS2 = cell(1,ND); doriEucBS2 = cell(1,ND); dsfEucBS2 = cell(1,ND); dcolorEucBS2 = cell(1,ND);
for d = 1:ND
    for R = 1:length(adom)
        doriBS2{d} = [doriBS2{d}; doriBS{d,R}];  %append simulation across ROIs
        dsfBS2{d} = [dsfBS2{d}; dsfBS{d,R}];
        dphaseBS2{d} = [dphaseBS2{d}; dphaseBS{d,R}];
        doriEucBS2{d} = [doriEucBS2{d}; doriEucBS{d,R}];
        dsfEucBS2{d} = [dsfEucBS2{d}; dsfEucBS{d,R}];
        dcolorEucBS2{d} = [dcolorEucBS2{d}; dcolorEucBS{d,R}];
    end
end


for d = 1:ND
    Ngrabs = length(doriBS2{d}(:,1));

    %Mean
    medoriBS_trial{d} = mean(doriBS2{d}); %mean of each trial
    medsfBS_trial{d} = mean(dsfBS2{d});
    medphaseBS_trial{d} = mean(dphaseBS2{d});
    medoriEucBS_trial{d} = mean(doriEucBS2{d}); %mean of each trial
    medsfEucBS_trial{d} = mean(dsfEucBS2{d});
    medcolorEucBS_trial{d} = mean(dcolorEucBS2{d});

    medoriBS(d) = mean(medoriBS_trial{d});  %total mean for each distance
    medsfBS(d) = mean(medsfBS_trial{d});
    medphaseBS(d) = mean(medphaseBS_trial{d});
    medoriEucBS(d) = mean(medoriEucBS_trial{d});  %total mean for each distance
    medsfEucBS(d) = mean(medsfEucBS_trial{d});
    medcolorEucBS(d) = mean(medcolorEucBS_trial{d});

    %Stand Dev
    stdoriBS_trial{d} = std(doriBS2{d})/sqrt(Ngrabs); %sig of each trial
    stdsfBS_trial{d} = std(dsfBS2{d})/sqrt(Ngrabs);
    stdphaseBS_trial{d} = std(dphaseBS2{d})/sqrt(Ngrabs);
    stdoriEucBS_trial{d} = std(doriEucBS2{d})/sqrt(Ngrabs); %sig of each trial
    stdsfEucBS_trial{d} = std(dsfEucBS2{d})/sqrt(Ngrabs);
    stdcolorEucBS_trial{d} = std(dcolorEucBS2{d})/sqrt(Ngrabs);

    stdoriBS(d) = mean(stdoriBS_trial{d});  %mean sig across trials, for each distance
    stdsfBS(d) = mean(stdsfBS_trial{d});
    stdphaseBS(d) = mean(stdphaseBS_trial{d});
    stdoriEucBS(d) = mean(stdoriEucBS_trial{d});  %mean sig across trials, for each distance
    stdsfEucBS(d) = mean(stdsfEucBS_trial{d});
    stdcolorEucBS(d) = mean(stdcolorEucBS_trial{d});

end


for d = 1:length(dori)
    
    %preference difference
%     medratio_ori(d) = medoriBS(d)/mean(dori{d});
%     medratio_sf(d) = medsfBS(d)/mean(dsf{d});
%     
%     stdratio_ori(d) = stdoriBS(d)/mean(dori{d});
%     stdratio_sf(d) = stdsfBS(d)/mean(dsf{d});
    
    medratio_ori(d) = medoriBS(d);  %already divided by the mean above
    medratio_sf(d) = medsfBS(d);
    medratio_phase(d) = medphaseBS(d);
    
    stdratio_ori(d) = stdoriBS(d);
    stdratio_sf(d) = stdsfBS(d);
    stdratio_phase(d) = stdphaseBS(d);
    
    id = find(medoriBS_trial{d} > 1);
    percOri(d) = length(id)/Nsim;
    id = find(medsfBS_trial{d} > 1);
    percSf(d) = length(id)/Nsim;
    id = find(medphaseBS_trial{d} > 1);
    percPhase(d) = length(id)/Nsim;
    
    %Euclidean difference
%     medratio_oriEuc(d) = medoriEucBS(d)/mean(doriEuc{d});
%     medratio_sfEuc(d) = medsfEucBS(d)/mean(dsfEuc{d});
%     
%     stdratio_oriEuc(d) = stdoriEucBS(d)/mean(doriEuc{d});
%     stdratio_sfEuc(d) = stdsfEucBS(d)/mean(dsfEuc{d});
    
    medratio_oriEuc(d) = medoriEucBS(d);
    medratio_sfEuc(d) = medsfEucBS(d);
    medratio_colorEuc(d) = medcolorEucBS(d);
    
    stdratio_oriEuc(d) = stdoriEucBS(d);
    stdratio_sfEuc(d) = stdsfEucBS(d);
    stdratio_colorEuc(d) = stdcolorEucBS(d);
    
    id = find(medoriEucBS_trial{d} > mean(doriEuc{d}));
    percOriEuc(d) = length(id)/Nsim;
    id = find(medsfEucBS_trial{d} > mean(dsfEuc{d}));
    percSfEuc(d) = length(id)/Nsim;
    id = find(medcolorEucBS_trial{d} > mean(dcolorEuc{d}));
    percColorEuc(d) = length(id)/Nsim;
    
end


%%
%Plot diff preference
figure, errorbar([distdom' distdom' distdom'],([medratio_ori' medratio_sf'  medratio_phase']),([stdratio_ori' stdratio_sf'  stdratio_sf']))
hold on, plot([distdom(1) distdom(end)],[1 1],'--k')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)    

ylabel('mean(ShuffledDist)/mean(ActualDist)')
xlabel('Distance')
title('Difference in preferred')
legend('ori','spfreq','phase')

hold on
for i = 1:length(dori)
   if percOri(i)<.01 | percOri(i)>.99 
       plot(distdom(i)+2,(medratio_ori(i))+.1,'*b')
   end
   if percSf(i)<.01  | percSf(i)>.99 
       plot(distdom(i)+2,(medratio_sf(i)+.1),'*g')
   end
end
ylim([.5 medratio_ori(1)+.5])

%Plot Euclidian
figure, errorbar([distdom' distdom' distdom' distdom'],([medratio_oriEuc' medratio_sfEuc' medratio_colorEuc' medratio_phase']),([stdratio_oriEuc' stdratio_sfEuc' stdratio_colorEuc' stdratio_sf']))
hold on, plot([distdom(1) distdom(end)],[1 1],'--k')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)    

ylabel('mean(ShuffledDist)/mean(ActualDist)')
xlabel('Distance')
title('Euclidean')
legend('ori','spfreq','color','phase')

hold on
for i = 1:length(dori)
   if percOriEuc(i)<.01 | percOriEuc(i)>.99 
       plot(distdom(i)+2,(medratio_oriEuc(i))+.1,'*b')
   end
   if percSfEuc(i)<.01  | percSfEuc(i)>.99 
       plot(distdom(i)+2,(medratio_sfEuc(i))+.1,'*g')
   end
end
ylim([.5 medratio_oriEuc(1)+.5])
    
%% Plot un-normalized histograms
ND = length(distdom);

figure
for i = 1:ND
    
    subplot(3,ND,i)
    dom = 0:10:90;
    [h] = histc(dori{i},dom);
    bar(dom,h/sum(h))
    xlabel('dori')
    xlim([-10 100]), ylim([0 .8])
    set(gca,'Xtick',dom)
    title(['mu = ' num2str(45/mean(dori{i}))])
    
    subplot(3,ND,i+ND)    
    dom = [0:.5:5];
    [h] = histc(dsf{i},dom);
    bar(dom,h/sum(h))
    %plotBarGraph(dsf{i},0,[0:.5:4])
    xlabel('dsf')
    xlim([-.5 4.5]), ylim([0 .5])
    title(['mu = ' num2str(mean(dsf{i}))])
    
    subplot(3,ND,i+2*ND)    
    dom = [0:20:180];
    [h] = histc(dphase{i},dom);
    bar(dom,h/sum(h))
    %plotBarGraph(dphase{i},0,[10:20:170])
    xlabel('dphase')
    xlim([-20 200]), ylim([0 .3])
    title(['mu = ' num2str(90/mean(dphase{i}))])
       
    muIm(1,i) = mean(dori{i});
    muIm(2,i) = mean(dsf{i});
    muIm(3,i) = mean(dphase{i});
        
end

    

for i = 1:ND
    dc(i) = mean(dcolorEuc{i});
    do(i) = mean(doriEuc{i});
    ds(i) = mean(dsfEuc{i});
    
end
figure,plot(dc/nanmean(dcolorEucAll))
hold on,plot(do/nanmean(doriEucAll))
hold on,plot(ds/nanmean(dsfEucAll))