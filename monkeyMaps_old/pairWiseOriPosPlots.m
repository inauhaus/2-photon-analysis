function pairWiseOriPosPlots(dori,dpos,doriNorm,dposNorm,dist)

global PW

doriAll = []; dposAll = []; doriNormAll = []; dposNormAll = []; distAll = [];
for i = 1:length(dori)
    
    doriAll = [doriAll; dori{i}(:)];
    dposAll = [dposAll; dpos{i}(:)];
    doriNormAll = [doriNormAll; doriNorm{i}(:)];
    dposNormAll = [dposNormAll; dposNorm{i}(:)];
    distAll = [distAll; dist{i}(:)];
    
end

% H = [distAll ones(length(distAll),1)];
% p = inv(H'*H)*H'*doriAll;
% doriAll = doriAll - H*p;
% 
% p = inv(H'*H)*H'*dposAll;
% dposAll = dposAll - H*p;

%PW.Ddom = [0 prctile(distAll,33) prctile(distAll,66) max(distAll)];
%PW.Ddom = [0 prctile(distAll,20) prctile(distAll,40) prctile(distAll,60) prctile(distAll,80) max(distAll)];
PW.Ddom = [0 prctile(distAll,25) prctile(distAll,50) prctile(distAll,75) max(distAll)];
%PW.Ddom = [0 max(distAll)];

% Ddom = 0:40:max(PW.Ddom);
% for i = 1:length(Ddom)-1
%     
%     id = find(distAll>Ddom(i) & distAll<=Ddom(i+1));
%     doriAll(id) = doriAll(id) - mean(doriAll(id));
%     dposAll(id) = dposAll(id) - mean(dposAll(id));
%     doriNormAll(id) = doriNormAll(id) - mean(doriNormAll(id));
%     dposNormAll(id) = dposNormAll(id) - mean(dposNormAll(id));
%     
% end

clear dori dpos doriNorm dposNorm 
for i = 1:length(PW.Ddom)-1
    
    id = find(distAll>PW.Ddom(i) & distAll<=PW.Ddom(i+1));
    dori{i} = doriAll(id);
    dpos{i} = dposAll(id);
    doriNorm{i} = doriNormAll(id);
    dposNorm{i} = dposNormAll(id);
    
end

%%  Plot 2D histogram for each distance
ND = length(dori);
figure
for i = 1:ND
    [mat xdom ydom] = smoothscatter(dori{i},dpos{i},.01,.01);

    matnorm = ones(length(mat(:,1)),1)*sum(mat);
    %matnorm = sum(mat,2)*ones(1,length(mat(1,:)));
    %mat = mat./matnorm;
    
    %mat = sum(mat,2)*sum(mat,1);

    subplot(2,ND,i),scatter(dori{i},dpos{i},'.k')
    %xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
    xlabel('dori'), ylabel('dpos')
    
    subplot(2,ND,i+ND),imagesc(xdom,ydom,(mat)), axis xy
    xlabel('dori (Norm)'), ylabel('dpos (Norm)'), 
    colormap pink
    
    [r p] = corrcoef(dori{i},dpos{i}); r = r(1,2); p = p(1,2);
    
    title(['P(dori,dpos|'  num2str(PW.Ddom(i))  ');   r = ' num2str(r) ';  p = ' num2str(p)])
end

%% dpos distribution given dori

yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];
%doridom = [0 15 30 45 90];
doridom = [prctile(doriAll,0) prctile(doriAll,25) prctile(doriAll,50) prctile(doriAll,75) prctile(doriAll,100)];
ns = length(doridom)-1;
clear prc mu_dpos
for Did = 1:length(PW.Ddom)-1  %loop through pairwise cortical distance

    posdom = linspace(0,max(dpos{Did})+.0001,10);
    
    xaxdom = doridom(1:end-1);
    
    clear Ppos
    %Ppos = zeros(length(posdom),Nori);
    for i = 1:length(doridom)-1  %loop through pairwise orientation difference
        
        id = find(dori{Did} >= doridom(i) & dori{Did} < doridom(i+1));
        dum = histc(dpos{Did}(id),posdom)'; dum = dum(1:end-1);
        Ppos(:,i) = dum;
        Ppos(:,i) = Ppos(:,i)/sum(Ppos(:,i));
        
        mu_dpos(Did,i) = nanmean(dpos{Did}(id));
        sig_dpos(Did,i) = nanstd(dpos{Did}(id))/sqrt(length(find(~isnan(dpos{Did}(id)))));

        for q = 1:length(prcdom)
            prc{q}(Did,i) = prctile(dpos{Did}(id),prcdom(q)*100);
        end
        
    end
    
    mu_dpos = mu_dpos(:,1:ns);
    cumPpos = [zeros(1,ns); cumsum(Ppos(:,1:ns))];
    dom = 0:length(cumPpos(:,1))-1;
    domI = linspace(dom(1),dom(end),length(dom)*10);
    cumPposI = interp1(dom,cumPpos,domI);
    
    figure(5)
    subplot(2,3,Did)
    plot(posdom,cumPpos,'-o'), ylim([0 1.2]), xlabel('dpos'), 
    title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])
    
    figure(6)    
    subplot(1,length(prcdom)+2,1)
    imagesc(xaxdom,yaxdom,mu_dpos)
    %set(gca,'Xtick',[]), set(gca,'Ytick',[])
    ylabel('Cortical Distance'), xlabel('dori')
    title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori
% 
%             [dum id] = min(abs(cumPposI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %posreq percentile location, as a function of dOri and cortical distance
% 
%         end
        subplot(1,length(prcdom)+2,q+1)
        imagesc(xaxdom,yaxdom,prc{q})
        %set(gca,'Xtick',[]), set(gca,'Ytick',[])
        %plot(prc{q}')
        ylabel('Cortical Distance'), xlabel('dori')
        title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(doridom)-1    
   labs{i} = [num2str(round(doridom(i))) ' to '  num2str(round(doridom(i+1)))];
end

figure(6)
subplot(1,length(prcdom)+2,length(prcdom)+2)
errorbar(mu_dpos',sig_dpos')
set(gca,'XTick',1:(length(doridom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dori range'), ylabel('E(dpos)')
%% dpos distribution given dori (Norm)


yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];
%doridom = [0 .1 .2 .3 1];
doridom = [prctile(doriNormAll,0) prctile(doriNormAll,25) prctile(doriNormAll,50) prctile(doriNormAll,75) prctile(doriNormAll,100)];
ns = length(doridom)-1;
clear prc mu_dpos
for Did = 1:length(PW.Ddom)-1  %loop through pairwise distance
    
    %doridom = [prctile(doriNorm{Did},0) prctile(doriNorm{Did},25) prctile(doriNorm{Did},50) prctile(doriNorm{Did},75) prctile(doriNorm{Did},100)];

    %posdom = linspace(0,1.0001,10);
    posdom = linspace(0,.3,10);

    xaxdom = doridom(1:end-1);
    
    clear Ppos
    %Ppos = zeros(length(posdom),Nori);
    for i = 1:length(doridom)-1

        id = find(doriNorm{Did} >= doridom(i) & doriNorm{Did} < doridom(i+1));
        dum = histc(dposNorm{Did}(id),posdom); dum = dum(1:end-1);
        Ppos(:,i) = dum;
        Ppos(:,i) = Ppos(:,i)/sum(Ppos(:,i));
        mu_dpos(Did,i) = nanmean(dposNorm{Did}(id));
        sig_dpos(Did,i) = nanstd(dposNorm{Did}(id))/sqrt(length(find(~isnan(dposNorm{Did}(id)))));  
        
        for q = 1:length(prcdom)
            prc{q}(Did,i) = prctile(dposNorm{Did}(id),prcdom(q)*100);
        end
    end
    
    mu_dpos = mu_dpos(:,1:ns);
    cumPpos = [zeros(1,ns); cumsum(Ppos(:,1:ns))];
    %cumPpos = Ppos;
    dom = 0:length(cumPpos(:,1))-1;
    domI = linspace(dom(1),dom(end),length(dom)*10);
    cumPposI = interp1(dom,cumPpos,domI);
    
    figure(3)
    subplot(2,3,Did)
    plot(posdom,cumPpos,'-o'), ylim([0 1.2]), xlabel('dposNorm') 
    title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])
    
    figure(4)    
    subplot(1,length(prcdom)+2,1)
    imagesc(xaxdom,yaxdom,mu_dpos)
    %set(gca,'Xtick',[]), set(gca,'Ytick',[])
    ylabel('Cortical Distance'), xlabel('doriNorm')
    title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori
% 
%             [dum id] = min(abs(cumPposI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %posreq percentile location, as a function of dOri and cortical distance
% 
%         end
        subplot(1,length(prcdom)+2,q+1)
        imagesc(xaxdom,yaxdom,prc{q})
        %set(gca,'Xtick',[]), set(gca,'Ytick',[])
        %plot(prc{q}')
        ylabel('Cortical Distance'), xlabel('doriNorm')
        title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(doridom)-1    
   labs{i} = [num2str(doridom(i)) ' to '  num2str(doridom(i+1))];
end
figure(4)
subplot(1,length(prcdom)+2,length(prcdom)+2)
errorbar(mu_dpos',sig_dpos')
set(gca,'XTick',1:(length(doridom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dori range (norm)'), ylabel('E(dpos) (norm)')
%% Now dori distribution given dpos (Norm)

yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];
%dposdom = [0 .2 .4 .6 1];
dposdom = [prctile(dposNormAll,0) prctile(dposNormAll,25) prctile(dposNormAll,50) prctile(dposNormAll,75) prctile(dposNormAll,100)];

ns = length(dposdom)-1;
clear prc mu_dori N
for Did = 1:length(PW.Ddom)-1  %loop through pairwise distance

    oridom = linspace(0,1.0001,10);

    xaxdom = dposdom(1:end-1);
    
    clear Pori
    %Pori = zeros(length(posdom),Nori);
    for i = 1:length(dposdom)-1

        id = find(dposNorm{Did} >= dposdom(i) & dposNorm{Did} < dposdom(i+1));
        oriDist{Did,i} = doriNorm{Did}(id);
        
        dum = histc(oriDist{Did,i},oridom);
        Pori(:,i) = dum(1:end-1);  %last element of histc is outside of range
        N(Did,i) = sum(Pori(:,i));
        Pori(:,i) = Pori(:,i)/N(Did,i);

        mu_dori(Did,i) = nanmean(doriNorm{Did}(id));
        sig_dori(Did,i) = nanmean(doriNorm{Did}(id))/sqrt(length(find(~isnan(doriNorm{Did}(id)))));
        
        for q = 1:length(prcdom)
            prc{q}(Did,i) = prctile(oriDist{Did,i},prcdom(q)*100);
        end
    end
    
    mu_dori = mu_dori(:,1:ns);
    cumPori = [zeros(1,ns); cumsum(Pori(:,1:ns))];
    %cumPpos = Ppos;
    dom = 0:length(cumPori(:,1))-1;
    domI = linspace(dom(1),dom(end),length(dom)*10);
    cumPoriI = interp1(dom,cumPori,domI);
    
    figure(34)
    subplot(2,3,Did)
    plot(cumPori,'-o'), ylim([0 1.2]), xlabel('doriNorm') 
    title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])

    figure(45)    
    subplot(1,length(prcdom)+2,1)
    imagesc(xaxdom,yaxdom,mu_dori)
    %set(gca,'Xtick',[]), set(gca,'Ytick',[])
    ylabel('Cortical Distance'), xlabel('dposNorm')
    title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori 
%             [dum id] = min(abs(cumPoriI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %posreq percentile location, as a function of dOri and cortical distance 
%         end
        subplot(1,length(prcdom)+2,q+1)

        imagesc(xaxdom,yaxdom,prc{q})
        %set(gca,'Xtick',[]), set(gca,'Ytick',[])
        %plot(prc{q}')
        ylabel('Cortical Distance'), xlabel('dposNorm')
        title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(dposdom)-1    
   labs{i} = [num2str(dposdom(i)) ' to '  num2str(dposdom(i+1))];
end

figure(45)
subplot(1,length(prcdom)+2,length(prcdom)+2)
errorbar(mu_dori',sig_dori')
set(gca,'XTick',1:(length(dposdom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dpos range (norm)'), ylabel('E(dori) (norm)')
%% dori distribution given dpos


yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];
%dposdom = [0 .7 1.4 2.1 10];
dposdom = [0 prctile(dposAll,25) prctile(dposAll,50) prctile(dposAll,75) 10];
ns = length(dposdom)-1;
clear prc mu_dori
for Did = 1:length(PW.Ddom)-1  %loop through pairwise distance

    oridom = linspace(0,90,10);

    xaxdom = doridom(1:end-1);
    
    clear Pori
    %Ppos = zeros(length(posdom),Nori);
    for i = 1:length(dposdom)-1

        id = find(dpos{Did} >= dposdom(i) & dpos{Did} < dposdom(i+1));
        
        dum = histc(dori{Did}(id),oridom); 
        Pori(:,i) = dum(1:end-1); %last element of histc is outside of range
        Pori(:,i) = Pori(:,i)/sum(Pori(:,i));
        
        mu_dori(Did,i) = nanmean(dori{Did}(id));
        mu_dsig(Did,i) = nanstd(dori{Did}(id))/sqrt(length(find(~isnan(dori{Did}(id)))));
        
        for q = 1:length(prcdom)
            prc{q}(Did,i) = prctile(dori{Did}(id),prcdom(q)*100);
        end
    end
    
    mu_dori = mu_dori(:,1:ns);
    cumPori = [zeros(1,ns); cumsum(Pori(:,1:ns))];
    %cumPpos = Ppos;
    dom = 0:length(cumPori(:,1))-1;
    domI = linspace(dom(1),dom(end),length(dom)*10);
    cumPoriI = interp1(dom,cumPori,domI);
    
    figure(38)
    subplot(2,3,Did)
    plot(cumPori,'-o'), ylim([0 1.2]), xlabel('dori') 
    title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])
    
    figure(39)    
    subplot(1,length(prcdom)+2,1)
    imagesc(xaxdom,yaxdom,mu_dori)
    %set(gca,'Xtick',[]), set(gca,'Ytick',[])
    ylabel('Cortical Distance'), xlabel('dpos')
    title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori
% 
%             [dum id] = min(abs(cumPposI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %posreq percentile location, as a function of dOri and cortical distance
% 
%         end
        subplot(1,length(prcdom)+2,q+1)
        imagesc(xaxdom,yaxdom,prc{q})
        %set(gca,'Xtick',[]), set(gca,'Ytick',[])
        %plot(prc{q}')
        ylabel('Cortical Distance'), xlabel('dpos')
        title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(dposdom)-1    
   labs{i} = [num2str(dposdom(i)) ' to '  num2str(dposdom(i+1))];
end

figure(39)
subplot(1,length(prcdom)+2,length(prcdom)+2)
errorbar(mu_dori',sig_dori')
set(gca,'XTick',1:(length(dposdom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dpos range'), ylabel('E(dori)')
%% anova
% dim = size(oriDist);
% vecAll = []; Distgroup = [];  posgroup = [];
% k = 1;
% for i = 1:dim(1)
%     for j = 1:dim(2)
%         vecAll = [vecAll; (oriDist{i,j}(:))];
%         Distgroup = [Distgroup i*ones(1,length(oriDist{i,j}(:)))];
%         posgroup = [posgroup j*ones(1,length(oriDist{i,j}(:)))];
%         k = k+1;
%     end
% end
% 
% [p,tbl,stats] = anovan(vecAll,{Distgroup posgroup})