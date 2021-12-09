function BSclustering_color(lumpref,LMphase,dlumAll,dLMphaseAll,distAll,animID,animIDdAll)

%modification of BSclustering_randpos (the one for randlum)

%%

%Get rid of NaN from lum and LMphase pairs
id = find(isnan(dlumAll.*dLMphaseAll));
dlumAll(id) = []; dLMphaseAll(id) = [];

distAll(id) = []; animIDdAll(id) = [];


%Get rid of NaN from lum and LMphase cells
id = find(~isnan(lumpref.*LMphase));
lumpref = lumpref(id);
LMphase = LMphase(id);

animID = animID(id);


prcdom = 0:20:100;
Dco = [];
for i = 1:length(prcdom)
    Dco = [Dco prctile(distAll,prcdom(i))];  %distanace domain cutoff points
end

for i = 1:length(Dco)-1    
   labs{i} = [num2str(round(Dco(i))) ' to '  num2str(round(Dco(i+1)))];
end

clear distdom dlum dLMphase animID_d
for i = 1:length(Dco)-1
    
    id = find(distAll>Dco(i) & distAll<=Dco(i+1));
    distdom(i) = mean(distAll(id));
    dlum{i} = dlumAll(id);
    dLMphase{i} = dLMphaseAll(id);    

    animID_d{i} = animIDdAll(id);
   
end

%%

dlumBS = []; dLMphaseBS = []; 

Ncell = length(lumpref);

Nsim = 10;

for d = 1:length(distdom)
    d
    Ngrabs = length(dlum{d});
    
    dlumBS{d} = zeros(Ngrabs,Nsim);
    dLMphaseBS{d} = zeros(Ngrabs,Nsim);
    
    for i = 1:Nsim
        for j = 1:Ngrabs

            rpick1 = round(1+rand(1)*(Ncell-1));  %randomly select a cell         
            rpick2 = rpick1;
            while rpick2 == rpick1
                rpick2 = round(1+rand(1)*(Ncell-1));  %randomly select a cell
            end
            
            dlumdum = abs(lumpref(rpick1) - lumpref(rpick2));

            dLMphasedum = abs(LMphase(rpick1) - LMphase(rpick2));
            
            dlumBS{d}(j,i) = dlumdum;
            dLMphaseBS{d}(j,i) = dLMphasedum;
            
        end

    end
     
end

%%

for d = 1:length(distdom)
    Ngrabs = length(dlum{d});

    %Mean
    medlumBS_trial{d} = mean(dlumBS{d}); %mean of each trial
    medLMphaseBS_trial{d} = mean(dLMphaseBS{d});

    medlumBS(d) = mean(medlumBS_trial{d});  %total mean for each distance
    medLMphaseBS(d) = mean(medLMphaseBS_trial{d});

    %Stand Dev
    stdlumBS_trial{d} = std(dlumBS{d})/sqrt(Ngrabs); %sig of each trial
    stdLMphaseBS_trial{d} = std(dLMphaseBS{d})/sqrt(Ngrabs);

    stdlumBS(d) = mean(stdlumBS_trial{d});  %mean sig across trials, for each distance
    stdLMphaseBS(d) = mean(stdLMphaseBS_trial{d});

end


%%

for d = 1:length(dlum)
    
    %preference difference
    medratio_lum(d) = medlumBS(d)/mean(dlum{d});
    medratio_LMphase(d) = medLMphaseBS(d)/mean(dLMphase{d});
    
    stdratio_lum(d) = stdlumBS(d)/mean(dlum{d});
    stdratio_LMphase(d) = stdLMphaseBS(d)/mean(dLMphase{d});
    
    id = find(medlumBS_trial{d} > mean(dlum{d}));
    perclum(d) = length(id)/Nsim;
    id = find(medLMphaseBS_trial{d} > mean(dLMphase{d}));
    percLMphase(d) = length(id)/Nsim;
    
end


perclum
percLMphase

%%
%Plot diff preference (n.b.  position bootstrap doesn't make sense here, so I got rid of it)
figure, errorbar([distdom' distdom'],([medratio_lum' medratio_LMphase']),([stdratio_lum' stdratio_LMphase']))
hold on, plot([distdom(1) distdom(end)],[1 1],'--k')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)    

ylabel('mean(ShuffledDist)/mean(ActualDist)')
xlabel('Distance')
title('Difference in preferred')
legend('L-M/L+M','LMphaseDiff')

hold on
for i = 1:length(dlum)
   if perclum(i)<.01 | perclum(i)>.99 
       plot(distdom(i)+2,(medratio_lum(i))+.1,'*b')
   end
   if percLMphase(i)<.01  | percLMphase(i)>.99 
       plot(distdom(i)+2,(medratio_LMphase(i)+.1),'*g')
   end
end
ylim([.5 medratio_lum(1)+.5])

    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now do bootstrap by shuffling individual ROIs

dlumBS = []; dLMphaseBS = [];
adom = unique(animID);

for R = 1:length(adom)  %parse the values according to the ROI

    id = find(adom(R) == animID);
    lumpref_R{R} = lumpref(id);
    LMphase_R{R} = LMphase(id);    
    
end


for R = 1:length(adom)  %parse the pairings according to the ROI

    id = find(adom(R) == animIDdAll ); %id pairs in this ROI
    dlum_R{R} = dlumAll(id);
    dLMphase_R{R} = dLMphaseAll(id);
    
end

Nsim = 2;

for d = 1:length(distdom)  %loop each distance
    mudo = mean(dlum{d});
    mudLMphase = mean(dLMphase{d});

    d
    for R = 1:length(adom)  %loop each ROI
        
        Ncell = length(lumpref_R{R}); %number of cells in this ROI
        
        Npairs = length(dlum_R{R}); %number of pairs in this ROI
    
        ida = find(animID_d{d} == adom(R)); %id the pairs for this ROI and distance
        Ngrabs = length(ida);  %number of pairs in ROI at this distance
        
        dlumBS{d,R} = zeros(Ngrabs,Nsim);
        dLMphaseBS{d,R} = zeros(Ngrabs,Nsim);
        
        for i = 1:Nsim
            
            for j = 1:Ngrabs


                %pick two cells and compute the difference
                rpick1 = round(1+rand(1)*(Ncell-1));  %randomly select a cell in ROI
                rpick2 = rpick1;
                while rpick2 == rpick1
                    rpick2 = round(1+rand(1)*(Ncell-1));  %randomly select a cell in ROI
                end
% 
%                 dlumdum = abs(lumdiff(lumpref_R{R}(rpick1)*pi/180,lumpref_R{R}(rpick2)*pi/180)*180/pi); %degrees
%                 dLMphasedum = abs(log2(LMphase_R{R}(rpick1)/LMphase_R{R}(rpick2)));  %octaves
                
                %Pick a pair where difference was already computed (necessary for phase)
                rpick = round(1+rand(1)*(Npairs-1));                  
                
                dlumdum = dlum_R{R}(rpick);
                dLMphasedum = dLMphase_R{R}(rpick);   

                dlumBS{d,R}(j,i) = (dlumdum/mudo);
                dLMphaseBS{d,R}(j,i) = (dLMphasedum/mudLMphase);                
                
                %dlumBS{d,R}(j,i) = (dlumdum/mean(dlum{d}(ida))); %normalize by the mean within this distance
                %dLMphaseBS{d,R}(j,i) = (dLMphasedum/mean(dLMphase{d}(ida)));

            end
            
        end
    end
    
end

%%
ND = length(distdom); 
dlumBS2 = cell(1,ND); dLMphaseBS2 = cell(1,ND);
for d = 1:ND
    for R = 1:length(adom)
        dlumBS2{d} = [dlumBS2{d}; dlumBS{d,R}];  %append simulation across ROIs
        dLMphaseBS2{d} = [dLMphaseBS2{d}; dLMphaseBS{d,R}];
    end
end


for d = 1:ND
    Ngrabs = length(dlumBS2{d}(:,1));

    %Mean
    medlumBS_trial{d} = mean(dlumBS2{d}); %mean of each trial
    medLMphaseBS_trial{d} = mean(dLMphaseBS2{d});

    medlumBS(d) = mean(medlumBS_trial{d});  %total mean for each distance
    medLMphaseBS(d) = mean(medLMphaseBS_trial{d});

    %Stand Dev
    stdlumBS_trial{d} = std(dlumBS2{d})/sqrt(Ngrabs); %sig of each trial
    stdLMphaseBS_trial{d} = std(dLMphaseBS2{d})/sqrt(Ngrabs);

    stdlumBS(d) = mean(stdlumBS_trial{d});  %mean sig across trials, for each distance
    stdLMphaseBS(d) = mean(stdLMphaseBS_trial{d});

end


for d = 1:length(dlum)
    
    %preference difference
%     medratio_lum(d) = medlumBS(d)/mean(dlum{d});
%     medratio_LMphase(d) = medLMphaseBS(d)/mean(dLMphase{d});
%     
%     stdratio_lum(d) = stdlumBS(d)/mean(dlum{d});
%     stdratio_LMphase(d) = stdLMphaseBS(d)/mean(dLMphase{d});
    
    medratio_lum(d) = medlumBS(d);  %already divided by the mean above
    medratio_LMphase(d) = medLMphaseBS(d);
    
    stdratio_lum(d) = stdlumBS(d);
    stdratio_LMphase(d) = stdLMphaseBS(d);
    
    id = find(medlumBS_trial{d} > 1);
    perclum(d) = length(id)/Nsim;
    id = find(medLMphaseBS_trial{d} > 1);
    percLMphase(d) = length(id)/Nsim;
    
end


%%
%Plot diff preference
figure, errorbar([distdom' distdom'],([medratio_lum' medratio_LMphase']),([stdratio_lum' stdratio_LMphase']))
hold on, plot([distdom(1) distdom(end)],[1 1],'--k')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)    

ylabel('mean(ShuffledDist)/mean(ActualDist)')
xlabel('Distance')
title('Difference in preferred')
legend('L-M/L+M','LMphaseDiff')

hold on
for i = 1:length(dlum)
   if perclum(i)<.01 | perclum(i)>.99 
       plot(distdom(i)+2,(medratio_lum(i))+.1,'*b')
   end
   if percLMphase(i)<.01  | percLMphase(i)>.99 
       plot(distdom(i)+2,(medratio_LMphase(i)+.1),'*g')
   end
end
ylim([.5 medratio_lum(1)+.5])

    
%% Plot un-normalized histograms
ND = length(distdom);

figure
for i = 1:ND
    
    subplot(2,ND,i)
    dom = 0:.01:.1;
    [h] = histc(dlum{i},dom);
    bar(dom,h/sum(h))
    xlabel('d L-M/L+M')
    xlim([0 .1]), ylim([0 .8])
    set(gca,'Xtick',dom)
    title(['mu = ' num2str(45/mean(dlum{i}))])
    
    subplot(2,ND,i+ND)    
    dom = [0:2:90];
    [h] = histc(dLMphase{i}*180/pi,dom);
    bar(dom,h/sum(h))
    %plotBarGraph(dLMphase{i},0,[0:.5:4])
    xlabel('dLMphase')
    xlim([-.1 90]), ylim([0 .5])
    title(['mu = ' num2str(mean(dLMphase{i}))])

    
    muIm(1,i) = mean(dlum{i});
    muIm(2,i) = mean(dLMphase{i});

        
end


%%  Expected correlation coefficient and dlum/dLMphase at each distance

for i = 1:length(dlum)

    domu(i) = mean(dlum{i});
    dosig(i) = std(dlum{i})/sqrt(length(dlum{i}));

    dLMphasemu(i) = mean(dLMphase{i});
    dLMphasesig(i) = std(dLMphase{i})/sqrt(length(dLMphase{i}));


end


figure,
subplot(2,1,1)
errorbar(distdom,domu,dosig,'k'), ylabel('|dlum| (deg)')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)  
subplot(2,1,2)
errorbar(distdom,dLMphasemu,dLMphasesig,'k'), ylabel('|dLMphase| (octaves)')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)  

xlabel('Distance')


