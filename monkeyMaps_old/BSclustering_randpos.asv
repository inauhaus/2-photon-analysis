function BSclustering_randpos(oripref,size,pos,BWdiff,RF,doriAll,dsizeAll,dposAll,dBWdiffAll,dRFEucAll,distAll,animID,animIDdAll)

%modification of BSclustering5 (the one for randori)

%%

%Get rid of NaN from ori and size pairs
id = find(isnan(doriAll.*dsizeAll.*dposAll.*dBWdiffAll));
doriAll(id) = []; dsizeAll(id) = []; dposAll(id) = [];  dBWdiffAll(id) = [];
dRFEucAll(id) = []; 
distAll(id) = []; animIDdAll(id) = [];


%Get rid of NaN from ori and size cells
id = find(~isnan(oripref.*size.*pos.*BWdiff));
oripref = oripref(id);
size = size(id);
pos = pos(id);
BWdiff = BWdiff(id);

RF = RF(id);
animID = animID(id);


prcdom = 0:25:100;
Dco = [];
for i = 1:length(prcdom)
    Dco = [Dco prctile(distAll,prcdom(i))];  %distanace domain cutoff points
end

for i = 1:length(Dco)-1    
   labs{i} = [num2str(round(Dco(i))) ' to '  num2str(round(Dco(i+1)))];
end

clear distdom dori dsize dpos dBWdiff dRFEuc animID_d
for i = 1:length(Dco)-1
    
    id = find(distAll>Dco(i) & distAll<=Dco(i+1));
    distdom(i) = mean(distAll(id));
    dori{i} = doriAll(id);
    dsize{i} = dsizeAll(id);    
    dpos{i} = dposAll(id);
    dBWdiff{i} = dBWdiffAll(id);
    dRFEuc{i} = dRFEucAll(id);
    animID_d{i} = animIDdAll(id);
   
end

%%

doriBS = []; dsizeBS = []; dBWdiffBS = [];

Ncell = length(oripref);

Nsim = 10;

for d = 1:length(distdom)
    d
    Ngrabs = length(dori{d});
    
    doriBS{d} = zeros(Ngrabs,Nsim);
    dsizeBS{d} = zeros(Ngrabs,Nsim);
    dBWdiffBS{d} = zeros(Ngrabs,Nsim);

    dRFEucBS{d} = zeros(Ngrabs,Nsim);
    
    for i = 1:Nsim
        for j = 1:Ngrabs

            rpick1 = round(1+rand(1)*(Ncell-1));  %randomly select a cell         
            rpick2 = rpick1;
            while rpick2 == rpick1
                rpick2 = round(1+rand(1)*(Ncell-1));  %randomly select a cell
            end

            doridum = abs(oridiff(oripref(rpick1)*pi/180,oripref(rpick2)*pi/180)*180/pi); %degrees
            dsizedum = abs(log2(size(rpick1)/size(rpick2)));  %octaves
            dBWdiffdum = abs(BWdiff(rpick1)-BWdiff(rpick2));  %degrees
            
            doriBS{d}(j,i) = doridum;
            dsizeBS{d}(j,i) = dsizedum;
            dBWdiffBS{d}(j,i) = dBWdiffdum;
            
        end

    end
     
end

%%

for d = 1:length(distdom)
    Ngrabs = length(dori{d});

    %Mean
    medoriBS_trial{d} = mean(doriBS{d}); %mean of each trial
    medsizeBS_trial{d} = mean(dsizeBS{d});
    medBWdiffBS_trial{d} = mean(dBWdiffBS{d});

    medoriBS(d) = mean(medoriBS_trial{d});  %total mean for each distance
    medsizeBS(d) = mean(medsizeBS_trial{d});
    medBWdiffBS(d) = mean(medBWdiffBS_trial{d});

    %Stand Dev
    stdoriBS_trial{d} = std(doriBS{d})/sqrt(Ngrabs); %sig of each trial
    stdsizeBS_trial{d} = std(dsizeBS{d})/sqrt(Ngrabs);
    stdBWdiffBS_trial{d} = std(dBWdiffBS{d})/sqrt(Ngrabs);

    stdoriBS(d) = mean(stdoriBS_trial{d});  %mean sig across trials, for each distance
    stdsizeBS(d) = mean(stdsizeBS_trial{d});
    stdBWdiffBS(d) = mean(stdBWdiffBS_trial{d});
end


%%

for d = 1:length(dori)
    
    %preference difference
    medratio_ori(d) = medoriBS(d)/mean(dori{d});
    medratio_size(d) = medsizeBS(d)/mean(dsize{d});
    medratio_BWdiff(d) = medBWdiffBS(d)/mean(dBWdiff{d});
    
    stdratio_ori(d) = stdoriBS(d)/mean(dori{d});
    stdratio_size(d) = stdsizeBS(d)/mean(dsize{d});
    stdratio_BWdiff(d) = stdBWdiffBS(d)/mean(dBWdiff{d});
    
    id = find(medoriBS_trial{d} > mean(dori{d}));
    percOri(d) = length(id)/Nsim;
    id = find(medsizeBS_trial{d} > mean(dsize{d}));
    percsize(d) = length(id)/Nsim;
    id = find(medBWdiffBS_trial{d} > mean(dBWdiff{d}));
    percBWdiff(d) = length(id)/Nsim;
    
end


percOri
percsize
percBWdiff
%%
%Plot diff preference (n.b.  position bootstrap doesn't make sense here, so I got rid of it)
figure, errorbar([distdom' distdom' distdom'],([medratio_ori' medratio_size'  medratio_BWdiff']),([stdratio_ori' stdratio_size' stdratio_BWdiff']))
hold on, plot([distdom(1) distdom(end)],[1 1],'--k')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)    

ylabel('mean(ShuffledDist)/mean(ActualDist)')
xlabel('Distance')
title('Difference in preferred')
legend('ori','size','Off - On')

hold on
for i = 1:length(dori)
   if percOri(i)<.01 | percOri(i)>.99 
       plot(distdom(i)+2,(medratio_ori(i))+.1,'*b')
   end
   if percsize(i)<.01  | percsize(i)>.99 
       plot(distdom(i)+2,(medratio_size(i)+.1),'*g')
   end
   if percBWdiff(i)<.01  | percBWdiff(i)>.99 
       plot(distdom(i)+2,(medratio_BWdiff(i)+.1),'*r')
   end
end
ylim([.5 medratio_ori(1)+.5])

    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now do bootstrap by shuffling individual ROIs

doriBS = []; dsizeBS = []; dposBS = [];  dBWdiffBS = [];
adom = unique(animID);

for R = 1:length(adom)  %parse the values according to the ROI

    id = find(adom(R) == animID);
    oripref_R{R} = oripref(id);
    size_R{R} = size(id);    
    pos_R{R} = pos(id);
    BWdiff_R{R} = BWdiff(id);
    
    RF_R{R} = RF(id);
    
end


for R = 1:length(adom)  %parse the pairings according to the ROI

    id = find(adom(R) == animIDdAll ); %id pairs in this ROI
    dori_R{R} = doriAll(id);
    dsize_R{R} = dsizeAll(id);
    dpos_R{R} = dposAll(id);
    dBWdiff_R{R} = dBWdiffAll(id);
    
end

Nsim = 2;

for d = 1:length(distdom)  %loop each distance
    mudo = mean(dori{d});
    mudsize = mean(dsize{d});
    mudpos = mean(dpos{d});
    mudBWdiff = mean(dBWdiff{d});
    mudoE = mean(dRFEuc{d});

    d
    for R = 1:length(adom)  %loop each ROI
        
        Ncell = length(oripref_R{R}); %number of cells in this ROI
        
        Npairs = length(dori_R{R}); %number of pairs in this ROI
    
        ida = find(animID_d{d} == adom(R)); %id the pairs for this ROI and distance
        Ngrabs = length(ida);  %number of pairs in ROI at this distance
        
        doriBS{d,R} = zeros(Ngrabs,Nsim);
        dsizeBS{d,R} = zeros(Ngrabs,Nsim);
        dposBS{d,R} = zeros(Ngrabs,Nsim);
        dBWdiffBS{d,R} = zeros(Ngrabs,Nsim);
        dRFEucBS{d,R} = zeros(Ngrabs,Nsim);
        
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
%                 dsizedum = abs(log2(size_R{R}(rpick1)/size_R{R}(rpick2)));  %octaves
                
                %Pick a pair where difference was already computed (necessary for phase)
                rpick = round(1+rand(1)*(Npairs-1));                  
                
                doridum = dori_R{R}(rpick);
                dsizedum = dsize_R{R}(rpick);   
                dposdum = dpos_R{R}(rpick);
                dBWdiffdum = dBWdiff_R{R}(rpick);

                doriBS{d,R}(j,i) = (doridum/mudo);
                dsizeBS{d,R}(j,i) = (dsizedum/mudsize);                
                dposBS{d,R}(j,i) = (dposdum/mudpos);
                dBWdiffBS{d,R}(j,i) = (dBWdiffdum/mudBWdiff);
                
                %doriBS{d,R}(j,i) = (doridum/mean(dori{d}(ida))); %normalize by the mean within this distance
                %dsizeBS{d,R}(j,i) = (dsizedum/mean(dsize{d}(ida)));

                v1 = RF_R{R}{rpick1}(:); v2 = RF_R{R}{rpick2}(:);
                dRFEucdum = corrcoef(v1,v2); dRFEucdum = -dRFEucdum(1,2); dRFEucdum = (dRFEucdum+1)/2;

                
%                 dRFEucBS{d,R}(j,i) = dRFEucdum/mudoE;
                
                dRFEucBS{d,R}(j,i) = dRFEucdum;
                
                %dRFEucBS{d,R}(j,i) = (dRFEucdum/mean(dRFEuc{d}(ida))); %normalize by the mean within this distance

            end
            
        end
    end
    
end

%%
ND = length(distdom); 
doriBS2 = cell(1,ND); dsizeBS2 = cell(1,ND); dposBS2 = cell(1,ND); dBWdiffBS2 = cell(1,ND); dRFEucBS2 = cell(1,ND);
for d = 1:ND
    for R = 1:length(adom)
        doriBS2{d} = [doriBS2{d}; doriBS{d,R}];  %append simulation across ROIs
        dsizeBS2{d} = [dsizeBS2{d}; dsizeBS{d,R}];
        dposBS2{d} = [dposBS2{d}; dposBS{d,R}];
        dBWdiffBS2{d} = [dBWdiffBS2{d}; dBWdiffBS{d,R}];
        dRFEucBS2{d} = [dRFEucBS2{d}; dRFEucBS{d,R}];
    end
end


for d = 1:ND
    Ngrabs = length(doriBS2{d}(:,1));

    %Mean
    medoriBS_trial{d} = mean(doriBS2{d}); %mean of each trial
    medsizeBS_trial{d} = mean(dsizeBS2{d});
    medposBS_trial{d} = mean(dposBS2{d});
    medBWdiffBS_trial{d} = mean(dBWdiffBS2{d});
    medRFEucBS_trial{d} = mean(dRFEucBS2{d}); %mean of each trial

    medoriBS(d) = mean(medoriBS_trial{d});  %total mean for each distance
    medsizeBS(d) = mean(medsizeBS_trial{d});
    medposBS(d) = mean(medposBS_trial{d});
    medBWdiffBS(d) = mean(medBWdiffBS_trial{d});
    medRFEucBS(d) = mean(medRFEucBS_trial{d});  %total mean for each distance

    %Stand Dev
    stdoriBS_trial{d} = std(doriBS2{d})/sqrt(Ngrabs); %sig of each trial
    stdsizeBS_trial{d} = std(dsizeBS2{d})/sqrt(Ngrabs);
    stdposBS_trial{d} = std(dposBS2{d})/sqrt(Ngrabs);
    stdBWdiffBS_trial{d} = std(dBWdiffBS2{d})/sqrt(Ngrabs);
    stdRFEucBS_trial{d} = std(dRFEucBS2{d})/sqrt(Ngrabs); %sig of each trial

    stdoriBS(d) = mean(stdoriBS_trial{d});  %mean sig across trials, for each distance
    stdsizeBS(d) = mean(stdsizeBS_trial{d});
    stdposBS(d) = mean(stdposBS_trial{d});
    stdBWdiffBS(d) = mean(stdBWdiffBS_trial{d});
    stdRFEucBS(d) = mean(stdRFEucBS_trial{d});  %mean sig across trials, for each distance

end


for d = 1:length(dori)
    
    %preference difference
%     medratio_ori(d) = medoriBS(d)/mean(dori{d});
%     medratio_size(d) = medsizeBS(d)/mean(dsize{d});
%     
%     stdratio_ori(d) = stdoriBS(d)/mean(dori{d});
%     stdratio_size(d) = stdsizeBS(d)/mean(dsize{d});
    
    medratio_ori(d) = medoriBS(d);  %already divided by the mean above
    medratio_size(d) = medsizeBS(d);
    medratio_pos(d) = medposBS(d);
    medratio_BWdiff(d) = medBWdiffBS(d);
    
    stdratio_ori(d) = stdoriBS(d);
    stdratio_size(d) = stdsizeBS(d);
    stdratio_pos(d) = stdposBS(d);
    stdratio_BWdiff(d) = stdBWdiffBS(d);
    
    id = find(medoriBS_trial{d} > 1);
    percOri(d) = length(id)/Nsim;
    id = find(medsizeBS_trial{d} > 1);
    percsize(d) = length(id)/Nsim;
    id = find(medposBS_trial{d} > 1);
    percpos(d) = length(id)/Nsim;
    id = find(medBWdiffBS_trial{d} > 1);
    percBWdiff(d) = length(id)/Nsim;
    
    %Euclidean difference
    medratio_RFEuc(d) = medRFEucBS(d)/mean(dRFEuc{d});
    
    stdratio_RFEuc(d) = stdRFEucBS(d)/mean(dRFEuc{d});
    
    id = find(medRFEucBS_trial{d} > mean(dRFEuc{d}));
    percRFEuc(d) = length(id)/Nsim;
    
end


%%
%Plot diff preference
figure, errorbar([distdom' distdom' distdom' distdom'],([medratio_ori' medratio_size' medratio_pos' medratio_BWdiff']),([stdratio_ori' stdratio_size' stdratio_pos' stdratio_BWdiff']))
hold on, plot([distdom(1) distdom(end)],[1 1],'--k')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)    

ylabel('mean(ShuffledDist)/mean(ActualDist)')
xlabel('Distance')
title('Difference in preferred')
legend('ori','size','position','Off-On')

hold on
for i = 1:length(dori)
   if percOri(i)<.01 | percOri(i)>.99 
       plot(distdom(i)+2,(medratio_ori(i))+.1,'*b')
   end
   if percsize(i)<.01  | percsize(i)>.99 
       plot(distdom(i)+2,(medratio_size(i)+.1),'*g')
   end
   if percpos(i)<.01  | percpos(i)>.99 
       plot(distdom(i)+2,(medratio_pos(i)+.1),'*r')
   end
   if percBWdiff(i)<.01  | percBWdiff(i)>.99 
       plot(distdom(i)+2,(medratio_BWdiff(i)+.1),'*c')
   end
end
ylim([.5 medratio_ori(1)+.5])

%Plot Euclidian
figure, errorbar(distdom',medratio_RFEuc',stdratio_RFEuc')
hold on, plot([distdom(1) distdom(end)],[1 1],'--k')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)    

ylabel('mean(ShuffledDist)/mean(ActualDist)')
xlabel('Distance')
title('Euclidean')
legend('RF overlap')

hold on
for i = 1:length(dori)
   if percRFEuc(i)<.01 | percRFEuc(i)>.99 
       plot(distdom(i)+2,(medratio_RFEuc(i))+.1,'*b')
   end
end
ylim([.5 medratio_RFEuc(1)+.5])
    
%% Plot un-normalized histograms
ND = length(distdom);

figure
for i = 1:ND
    
    subplot(4,ND,i)
    dom = 0:10:90;
    [h] = histc(dori{i},dom);
    bar(dom,h/sum(h))
    xlabel('dori')
    xlim([-10 100]), ylim([0 .8])
    set(gca,'Xtick',dom)
    title(['mu = ' num2str(45/mean(dori{i}))])
    
    subplot(4,ND,i+ND)    
    dom = [0:.05:75];
    [h] = histc(dsize{i},dom);
    bar(dom,h/sum(h))
    %plotBarGraph(dsize{i},0,[0:.5:4])
    xlabel('dsize')
    xlim([-.1 .75]), ylim([0 .5])
    title(['mu = ' num2str(mean(dsize{i}))])
    
    subplot(4,ND,i+2*ND)    
    dom = [0:.05:75];
    [h] = histc(dpos{i},dom);
    bar(dom,h/sum(h))
    %plotBarGraph(dsize{i},0,[0:.5:4])
    xlabel('dpos')
    xlim([-.1 .75]), ylim([0 .5])
    title(['mu = ' num2str(mean(dpos{i}))])
       
    subplot(4,ND,i+3*ND)    
    dom = [0:.05:75];
    [h] = histc(dBWdiff{i},dom);
    bar(dom,h/sum(h))
    %plotBarGraph(dsize{i},0,[0:.5:4])
    xlabel('dBWdiff')
    xlim([-.1 .75]), ylim([0 .5])
    title(['mu = ' num2str(mean(dBWdiff{i}))])
    
    muIm(1,i) = mean(dori{i});
    muIm(2,i) = mean(dsize{i});
    muIm(3,i) = mean(dpos{i});
    muIm(4,i) = mean(dBWdiff{i});
        
end
%%  Expected correlation coefficient and dpos/dori at each distance

for i = 1:length(dori)

    domu(i) = mean(dori{i});
    dosig(i) = std(dori{i})/sqrt(length(dori{i}));

    dsizemu(i) = mean(dsize{i});
    dsizesig(i) = std(dsize{i})/sqrt(length(dsize{i}));

    dum = dpos{i};
    dposmu(i) = mean(dum);
    dpossig(i) = std(dum)/sqrt(length(dum));

    dum = dBWdiff{i};
    dBWdiffmu(i) = mean(dum);
    dBWdiffsig(i) = std(dum)/sqrt(length(dum));

  
end


figure,
subplot(4,1,1)
errorbar(distdom,domu,dosig,'k'), ylabel('|dori| (deg)')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)  
subplot(4,1,2)
errorbar(distdom,dsizemu,dsizesig,'k'), ylabel('|dsize| (octaves)')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)  
subplot(4,1,3)
errorbar(distdom,dposmu,dpossig,'k'), ylabel('|dpos|  (deg)')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)  
subplot(4,1,4)
errorbar(distdom,dBWdiffmu,dBWdiffsig,'k'), ylabel('|dBWdiff|')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)  

xlabel('Distance')

