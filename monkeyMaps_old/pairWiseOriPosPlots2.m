function pairWiseOriPosPlots2(doriAll,dposAll,doriNormAll,dposNormAll,distAll,animID)

global PW

%normalize distribution for each animal?
% animdom = unique(animID);
% for i = 1:length(animdom)
%     id = find(animID == animdom(i));
%     dposAll(id) = (dposAll(id)-nanmean(dposAll(id)))/nanstd(dposAll(id));
%     dposNormAll(id) = (dposNormAll(id)-nanmean(dposNormAll(id)))/nanstd(dposNormAll(id));
% end

PW.Ddom = [0 prctile(distAll,33) prctile(distAll,66) max(distAll)];
%PW.Ddom = [0 prctile(distAll,50) max(distAll)];
PW.Ddom = round(PW.Ddom);

clear dori dpos doriNorm dposNorm 
for i = 1:length(PW.Ddom)-1
    
    id = find(distAll>PW.Ddom(i) & distAll<=PW.Ddom(i+1) & ~isnan(doriAll) & ~isnan(dposAll) & ~isnan(doriNormAll) & ~isnan(dposNormAll));
    dori{i} = doriAll(id);
    dpos{i} = dposAll(id);
    doriNorm{i} = doriNormAll(id);
    dposNorm{i} = dposNormAll(id);
    
end

%%  Plot 2D histogram for each distance
ND = length(doriNorm);
figure(99)
for i = 1:ND
    [mat xdom ydom] = smoothscatter(doriNorm{i},dposNorm{i},.08,.004,[0 90],[0 .7]);
    mat = mat/sum(mat(:));

    %matnorm = ones(length(mat(:,1)),1)*sum(mat);
    %matnorm = sum(mat,2)*ones(1,length(mat(1,:)));
    %mat = mat./matnorm;
    
    %mat = sum(mat,2)*sum(mat,1);

    %subplot(2,ND,i),scatter(doriNorm{i},dposNorm{i},'.k')
    %xlabel('dori'), ylabel('dpos')
    
    subplot(1,ND,i),imagesc(xdom,ydom,-(mat).^.35), axis xy
    xlabel('dori'), ylabel('dpos (normalized)')
    colormap gray
    
    set(gca,'Xtick',[0 45 90])
    
    [r p] = corrcoef(doriNorm{i},dposNorm{i});
    title(['P(dori,dpos|'  num2str(PW.Ddom(i))  ');   r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2))])
    
    axis square
end


%% dpos distribution given dori (Norm)

yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];

clear prc mu_dpos sig_dpos xdom
for Did = 1:length(PW.Ddom)-1  %loop through pairwise distance
    
    %doridom = [-1 prctile(doriNorm{Did},25) prctile(doriNorm{Did},50) prctile(doriNorm{Did},75) 1];
    doridom = [prctile(doriNorm{Did},0) prctile(doriNorm{Did},33) prctile(doriNorm{Did},66) prctile(doriNorm{Did},100)];
    ns = length(doridom)-1;
    
    posdom = linspace(-1,1.0001,10);

    xaxdom = doridom(1:end-1);
    
    clear Ppos
    %Ppos = zeros(length(posdom),Nori);
    for i = 1:length(doridom)-1

        id = find(doriNorm{Did} >= doridom(i) & doriNorm{Did} < doridom(i+1));
        
        dum = histc(dposNorm{Did}(id),posdom); dum = dum(1:end-1);
        Ppos(:,i) = dum;
        Ppos(:,i) = Ppos(:,i)/sum(Ppos(:,i));
        
        idno = find(~isnan(dposNorm{Did}(id)));
        mu_dpos(Did,i) = nanmean(dposNorm{Did}(id));
        sig_dpos(Did,i) = nanstd(dposNorm{Did}(id))/sqrt(length(idno));
        
        xdom{Did}(i) = nanmean(doriNorm{Did}(id));
        
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
    
%     figure(3)
%     subplot(2,3,Did)
%     plot(cumPpos,'-o'), ylim([0 1.2]), xlabel('dposNorm') 
%     title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])
%     
%     figure(4)    
%     subplot(1,length(prcdom)+1,1)
%     imagesc(xaxdom,yaxdom,mu_dpos)
%     %set(gca,'Xtick',[]), set(gca,'Ytick',[])
%     ylabel('Cortical Distance'), xlabel('doriNorm')
%     title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori
% 
%             [dum id] = min(abs(cumPposI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %posreq percentile location, as a function of dOri and cortical distance
% 
%         end

%         subplot(1,length(prcdom)+1,q+1)
%         imagesc(xaxdom,yaxdom,prc{q})
%         ylabel('Cortical Distance'), xlabel('doriNorm')
%         title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(doridom)-1    
   labs{i} = [num2str(doridom(i)) ' to '  num2str(doridom(i+1))];
end

figure, bar(mu_dpos')
hold on,errorbar(mu_dpos',sig_dpos','.k')
set(gca,'XTick',1:(length(doridom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dori range (norm)'), ylabel('E(dpos) (norm)')

NDis = length(PW.Ddom)-1;
figure(99)
for i = 1:NDis
    subplot(1,3,i)
    hold on, plot(xdom{i},mu_dpos(i,:),'.-b')
    for j = 1:length(sig_dpos)
        hold on
        plot([xdom{i}(j) xdom{i}(j)],[mu_dpos(i,j)-sig_dpos(i,j) ,mu_dpos(i,j)+sig_dpos(i,j)],'b')
    end
end


%% Now dori distribution given dpos (Norm)

yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];

clear prc mu_dori N xdom
for Did = 1:length(PW.Ddom)-1  %loop through pairwise distance

    %dposdom = [-1 prctile(dposNorm{Did},25) prctile(dposNorm{Did},50) prctile(dposNorm{Did},75) 1];
    dposdom = [prctile(dposNorm{Did},0) prctile(dposNorm{Did},33) prctile(dposNorm{Did},66) prctile(doriNorm{Did},100)];
    ns = length(dposdom)-1;
    
    oridom = linspace(0,90,5);

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

        idno = find(~isnan(doriNorm{Did}(id)));
        mu_dori(Did,i) = nanmean(doriNorm{Did}(id));
        sig_dori(Did,i) = nanstd(doriNorm{Did}(id))/sqrt(length(idno));
        
        xdom{Did}(i) = nanmean(dposNorm{Did}(id));
        
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
    
%     figure(34)
%     subplot(2,3,Did)
%     plot(cumPori,'-o'), ylim([0 1.2]), xlabel('doriNorm') 
%     title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])
% 
%     figure(45)    
%     subplot(1,length(prcdom)+1,1)
%     imagesc(xaxdom,yaxdom,mu_dori)
%     %set(gca,'Xtick',[]), set(gca,'Ytick',[])
%     ylabel('Cortical Distance'), xlabel('dposNorm')
%     title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori 
%             [dum id] = min(abs(cumPoriI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %posreq percentile location, as a function of dOri and cortical distance 
%         end
%         subplot(1,length(prcdom)+1,q+1)
% 
%         imagesc(xaxdom,yaxdom,prc{q})
%         ylabel('Cortical Distance'), xlabel('dposNorm')
%         title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(dposdom)-1    
   labs{i} = [num2str(dposdom(i)) ' to '  num2str(dposdom(i+1))];
end
figure, bar(mu_dori')
hold on,errorbar(mu_dori',sig_dori','.k')
set(gca,'XTick',1:(length(dposdom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dpos range (norm)'), ylabel('E(dori) (norm)')

NDis = length(PW.Ddom)-1;

figure(99)
for i = 1:NDis
    subplot(1,3,i)
    hold on, plot(mu_dori(i,:),xdom{i},'.-r')
    for j = 1:length(sig_dpos)
        hold on
        plot([mu_dori(i,j)-sig_dori(i,j) ,mu_dori(i,j)+sig_dori(i,j)],[xdom{i}(j) xdom{i}(j)],'r')
    end
end


%%%%%%%%%%%%%%%Now regular dori and dpos%%%%%%%%%%%%%


%%  Plot 2D histogram for each distance
ND = length(dori);
figure(98)
for i = 1:ND
    [mat xdom ydom] = smoothscatter(dori{i},dpos{i},.05,.005);

    %matnorm = ones(length(mat(:,1)),1)*sum(mat);
    %matnorm = sum(mat,2)*ones(1,length(mat(1,:)));
    %mat = mat./matnorm;
    
    %mat = sum(mat,2)*sum(mat,1);

    %subplot(2,ND,i),scatter(dori{i},dpos{i},'.k')
    %xlabel('dori'), ylabel('dpos')
    
    subplot(1,ND,i),imagesc(xdom,ydom,-(mat).^.35), axis xy
    xlabel('dori'), ylabel('dpos'),
    colormap gray
    
    set(gca,'Xtick',[0 45 90])
    
    [r p] = corrcoef(dori{i},dpos{i});
    title(['P(dori,dpos|'  num2str(PW.Ddom(i))  ');   r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2))])
    
    axis square

end




%% dpos distribution given dori

yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];

clear prc mu_dpos sig_dpos xdom
for Did = 1:length(PW.Ddom)-1  %loop through pairwise cortical distance

    %doridom = [0 prctile(dori{Did},25) prctile(dori{Did},50) prctile(dori{Did},75) 90];
    doridom = [0 prctile(dori{Did},33) prctile(dori{Did},66) 90];
    ns = length(doridom)-1;
    
    posdom = linspace(0,max(dpos{Did})+.0001,10);
    
    xaxdom = doridom(1:end-1);
    
    clear Ppos
    %Ppos = zeros(length(posdom),Nori);
    for i = 1:length(doridom)-1  %loop through pairwise orientation difference

        id = find(dori{Did} >= doridom(i) & dori{Did} < doridom(i+1));
        dum = histc(dpos{Did}(id),posdom)'; dum = dum(1:end-1);
        Ppos(:,i) = dum;
        Ppos(:,i) = Ppos(:,i)/sum(Ppos(:,i));
        
        idno = find(~isnan(dpos{Did}(id)));
        mu_dpos(Did,i) = nanmean(dpos{Did}(id));
        sig_dpos(Did,i) = nanstd(dpos{Did}(id))/sqrt(length(idno));
        
        xdom{Did}(i) = nanmean(dori{Did}(id));

        for q = 1:length(prcdom)
            prc{q}(Did,i) = prctile(dpos{Did}(id),prcdom(q)*100);
        end
        
    end
    
    mu_dpos = mu_dpos(:,1:ns);
    cumPpos = [zeros(1,ns); cumsum(Ppos(:,1:ns))];
    dom = 0:length(cumPpos(:,1))-1;
    domI = linspace(dom(1),dom(end),length(dom)*10);
    cumPposI = interp1(dom,cumPpos,domI);
    
%     figure(5)
%     subplot(2,3,Did)
%     plot([posdom],cumPpos,'-o'), ylim([0 1.2]), xlabel('dpos'), 
%     title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])
%     
%     figure(6)    
%     subplot(1,length(prcdom)+1,1)
%     imagesc(xaxdom,yaxdom,mu_dpos)
%     ylabel('Cortical Distance'), xlabel('dori')
%     title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori
% 
%             [dum id] = min(abs(cumPposI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %posreq percentile location, as a function of dOri and cortical distance
% 
%         end

%         subplot(1,length(prcdom)+1,q+1)
%         imagesc(xaxdom,yaxdom,prc{q})
%         ylabel('Cortical Distance'), xlabel('dori')
%         title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(doridom)-1    
   labs{i} = [num2str(round(doridom(i))) ' to '  num2str(round(doridom(i+1)))];
end
figure,bar(mu_dpos')
hold on,errorbar(mu_dpos',sig_dpos','.k')
set(gca,'XTick',1:(length(doridom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dori range'), ylabel('E(dpos)')

NDis = length(PW.Ddom)-1;
figure(98)
for i = 1:NDis
    subplot(1,3,i)
    hold on, plot(xdom{i},mu_dpos(i,:),'.-b')
    for j = 1:length(sig_dpos)
        hold on
        plot([xdom{i}(j) xdom{i}(j)],[mu_dpos(i,j)-sig_dpos(i,j) ,mu_dpos(i,j)+sig_dpos(i,j)],'b')
    end
end


%% dori distribution given dpos


yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];
%dposdom = [0 .7 1.4 2.1 10];

clear prc mu_dori sig_dori xdom
for Did = 1:length(PW.Ddom)-1  %loop through pairwise distance

    %dposdom = [prctile(dpos{Did},0) prctile(dpos{Did},25) prctile(dpos{Did},50) prctile(dpos{Did},75) prctile(dpos{Did},6100)];
    dposdom = [prctile(dpos{Did},0) prctile(dpos{Did},33) prctile(dpos{Did},66) prctile(dpos{Did},100)];
    ns = length(dposdom)-1;
    
    oridom = linspace(0,90,10);

    xaxdom = doridom(1:end-1);
    
    clear Pori mu
    %Ppos = zeros(length(posdom),Nori);
    for i = 1:length(dposdom)-1

        id = find(dpos{Did} >= dposdom(i) & dpos{Did} < dposdom(i+1));
        
        dum = histc(dori{Did}(id),oridom); 
        Pori(:,i) = dum(1:end-1); %last element of histc is outside of range
        Pori(:,i) = Pori(:,i)/sum(Pori(:,i));
        
        idno = find(~isnan(dori{Did}(id)));
        mu_dori(Did,i) = nanmean(dori{Did}(id));
        sig_dori(Did,i) = nanstd(dori{Did}(id))/sqrt(length(idno));
        
        xdom{Did}(i) = nanmean(dpos{Did}(id));
        
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
    
%     figure(38)
%     subplot(2,3,Did)
%     plot(cumPori,'-o'), ylim([0 1.2]), xlabel('dori') 
%     title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])
%     
%     figure(39)    
%     subplot(1,length(prcdom)+1,1)
%     imagesc(xaxdom,yaxdom,mu_dori)
%     ylabel('Cortical Distance'), xlabel('dpos')
%     title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori
% 
%             [dum id] = min(abs(cumPposI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %posreq percentile location, as a function of dOri and cortical distance
% 
%         end

%         subplot(1,length(prcdom)+1,q+1)
%         imagesc(xaxdom,yaxdom,prc{q})
%         ylabel('Cortical Distance'), xlabel('dpos')
%         title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(dposdom)-1    
   labs{i} = [num2str(dposdom(i)) ' to '  num2str(dposdom(i+1))];
end
figure, bar(mu_dori')
hold on,errorbar(mu_dori',sig_dori','.k')
set(gca,'XTick',1:(length(dposdom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dpos range'), ylabel('E(dori)')

NDis = length(PW.Ddom)-1;

figure(98)
for i = 1:NDis
    subplot(1,3,i)
    hold on, plot(mu_dori(i,:),xdom{i},'.-r')
    for j = 1:length(sig_dpos)
        hold on
        plot([mu_dori(i,j)-sig_dori(i,j) ,mu_dori(i,j)+sig_dori(i,j)],[xdom{i}(j) xdom{i}(j)],'r')
    end
end
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