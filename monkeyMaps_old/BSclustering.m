function BSclustering(oripref,sfpref,tcori,tcsf,dori,dsf,doripair,dsfpair,dtcoripair,dtcsfpair,doriEuc,dsfEuc,dist)

%bootstrap analysis for ori/sf clustering

%%

doriAll = []; dsfAll = []; doripairAll = []; dsfpairAll = []; dtcoripairAll = []; dtcsfpairAll = []; doriEucAll = []; dsfEucAll = []; distAll = [];
for i = 1:length(dori)
    
    doriAll = [doriAll; dori{i}(:)];
    dsfAll = [dsfAll; dsf{i}(:)];
    doripairAll = [doripairAll; doripair{i}];
    dsfpairAll = [dsfpairAll; dsfpair{i}];
    dtcoripairAll = [dtcoripairAll; dtcoripair{i}];
    dtcsfpairAll = [dtcsfpairAll; dtcsfpair{i}];
    doriEucAll = [doriEucAll; doriEuc{i}(:)];
    dsfEucAll = [dsfEucAll; dsfEuc{i}(:)];
    distAll = [distAll; dist{i}(:)];
    
end

prcdom = 0:20:100;
PW.Ddom = [];
for i = 1:length(prcdom)
    PW.Ddom = [PW.Ddom prctile(distAll,prcdom(i))];
end

for i = 1:length(PW.Ddom)-1    
   labs{i} = [num2str(round(PW.Ddom(i))) ' to '  num2str(round(PW.Ddom(i+1)))];
end

clear dori dsf doripair dsfpair dtcoripair dtcsfpair doriEuc dsfEuc 
for i = 1:length(PW.Ddom)-1
    
    id = find(distAll>PW.Ddom(i) & distAll<=PW.Ddom(i+1));
    distdom(i) = mean(distAll(id));
    dori{i} = doriAll(id);
    dsf{i} = dsfAll(id);
    doripair{i} = doripairAll(id,:);
    dsfpair{i} = dsfpairAll(id,:);
    dtcoripair{i} = dtcoripairAll(id,:);
    dtcsfpair{i} = dtcsfpairAll(id,:);
    doriEuc{i} = doriEucAll(id);
    dsfEuc{i} = dsfEucAll(id);
    
end
%%
% clear rori por rsf psf
% for d = 1:length(dori) 
%     d
% %     r = cxcorr(doripair{d}(:,1),doripair{d}(:,2),180);
% %     if r<0
% %         [r p] = circCorr(180-doripair{d}(:,1),doripair{d}(:,2),180,1);
% %         r = -r;
% %     else
% %         [r p] = circCorr(doripair{d}(:,1),doripair{d}(:,2),180,1);
% %     end
%     [r p] = circCorr(doripair{d}(:,1),doripair{d}(:,2),180,1);
%     rori(d) = r;  pori(d) = p;
%     
%     r = corrcoef(log(dsfpair{d}(:,1)),log(dsfpair{d}(:,2)));     r = r(1,2);
%     if r<0
%         [r p] = corrcoef_BS(-log(dsfpair{d}(:,1)),log(dsfpair{d}(:,2)),1,1);
%         r = -r;
%     else
%         [r p] = corrcoef_BS(log(dsfpair{d}(:,1)),log(dsfpair{d}(:,2)),1,1);
%     end
%     rsf(d) = r; psf(d) = p;
%     %rsf(d) = r(1,2); psf(d) = p(1,2);
% end
% 
% figure, plot(distdom,rori,'-ok')
% hold on, plot(distdom,rsf,'-or')
% ylabel('r'),xlabel('distance')
% legend('ori','sfreq')
% hold on
% for i = 1:length(dori)
%    if pori(i)<.01 
%        plot(distdom(i),rori(i),'*k')
%    end
%    if psf(i)<.01 
%        plot(distdom(i),rsf(i),'*r')
%    end
% end
% 
% hold off
%%
% Nori = length(dtcoripair{1}(1,:))/2;
% Nsf = length(dtcsfpair{1}(1,:))/2;
% clear rori pori rsf psf
% for d = 1:length(dori)
%     
%     oridum1 = dtcoripair{d}(:,1:Nori);  oridum2 = dtcoripair{d}(:,Nori+1:end);
%     
%     mi = min(oridum1,[],2); oridum1 = oridum1-mi*ones(1,Nori);
%     ma = sum(oridum1,2); oridum1 = oridum1./(ma*ones(1,Nori));
%     
%     [r p] = corrcoef(oridum1(:),oridum2(:));
%     rori(d) = r(1,2); pori(d) = p(1,2);
%     
%     sfdum1 = dtcsfpair{d}(:,1:Nsf);  sfdum2 = dtcsfpair{d}(:,Nsf+1:end);
%     mi = min(sfdum1,[],2); sfdum1 = sfdum1-mi*ones(1,Nsf);
%     ma = sum(sfdum1,2); sfdum1 = sfdum1./(ma*ones(1,Nsf));
%     
%     [r p] = corrcoef(sfdum1(:),sfdum2(:));
%     rsf(d) = r(1,2); psf(d) = p(1,2);
% end
% 
% figure, plot(distdom,rori,'-ok')
% hold on, plot(distdom,rsf,'-or')
% ylabel('r'),xlabel('distance')
% legend('ori','sfreq')
%%

doriBS = []; dsfBS = []; 

id = find(~isnan(oripref.*sfpref));
oripref = oripref(id);
sfpref = sfpref(id);
tcori = tcori(id,:);
tcsf = tcsf(id,:);

Ncell = length(oripref);

Nsim = 10;

for d = 1:length(PW.Ddom)-1
    Ngrabs = length(dori{d});
    doriBS{d} = zeros(Ngrabs,Nsim);
    dsfBS{d} = zeros(Ngrabs,Nsim);
    doriEucBS{d} = zeros(Ngrabs,Nsim);
    dsfEucBS{d} = zeros(Ngrabs,Nsim);
    for i = 1:Nsim
        i
        for j = 1:Ngrabs

            rpick1 = round(1+rand(1)*(Ncell-1));  %randomly select a cell
            rpick2 = round(1+rand(1)*(Ncell-1));  %randomly select a cell

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

        end

    end
    
    medoriBS_trial{d} = mean(doriBS{d}); %mean of each trial
    medsfBS_trial{d} = mean(dsfBS{d});
    medoriEucBS_trial{d} = mean(doriEucBS{d}); %mean of each trial
    medsfEucBS_trial{d} = mean(dsfEucBS{d});
    
    medoriBS(d) = mean(medoriBS_trial{d});  %total mean for each distance
    medsfBS(d) = mean(medsfBS_trial{d});
    medoriEucBS(d) = mean(medoriEucBS_trial{d});  %total mean for each distance
    medsfEucBS(d) = mean(medsfEucBS_trial{d});

    
end


%%

for d = 1:length(dori)
    
    %preference difference
    medratio_ori(d) = medoriBS(d)/mean(dori{d});
    medratio_sf(d) = medsfBS(d)/mean(dsf{d});
    
    id = find(medoriBS_trial{d} > mean(dori{d}));
    percOri(d) = length(id)/Nsim;
    id = find(medsfBS_trial{d} > mean(dsf{d}));
    percSf(d) = length(id)/Nsim;
    
    %Euclidean difference
    medratio_oriEuc(d) = medoriEucBS(d)/mean(doriEuc{d});
    medratio_sfEuc(d) = medsfEucBS(d)/mean(dsfEuc{d});
    
    id = find(medoriEucBS_trial{d} > mean(doriEuc{d}));
    percOriEuc(d) = length(id)/Nsim;
    id = find(medsfEucBS_trial{d} > mean(dsfEuc{d}));
    percSfEuc(d) = length(id)/Nsim;
end



percOri
percSf
%%
%Plot diff preference
figure, plot(distdom,log2([medratio_ori' medratio_sf']),'o-')
hold on, plot([distdom(1) distdom(end)],[0 0],'--k')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)    

ylabel('mean(ShuffledDist)/mean(ActualDist)')
xlabel('Distance')
title('Difference in preferred')
legend('ori','spfreq')

hold on
for i = 1:length(dori)
   if percOri(i)<.01 | percOri(i)>.99 
       plot(distdom(i),log2(medratio_ori(i)),'*b')
   end
   if percSf(i)<.01  | percSf(i)>.99 
       plot(distdom(i),log2(medratio_sf(i)),'*g')
   end
end

%Plot Euclidian
figure, plot(distdom,log2([medratio_oriEuc' medratio_sfEuc']),'o-')
hold on, plot([distdom(1) distdom(end)],[0 0],'--k')
set(gca,'XTick',distdom)
set(gca,'XTickLabel',labs)    

ylabel('mean(ShuffledDist)/mean(ActualDist)')
xlabel('Distance')
title('Euclidean')
legend('ori','spfreq')

hold on
for i = 1:length(dori)
   if percOriEuc(i)<.01 | percOriEuc(i)>.99 
       plot(distdom(i),log2(medratio_oriEuc(i)),'*b')
   end
   if percSfEuc(i)<.01  | percSfEuc(i)>.99 
       plot(distdom(i),log2(medratio_sfEuc(i)),'*g')
   end
end
    
    
