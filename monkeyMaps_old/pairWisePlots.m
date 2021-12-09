function pairWisePlots(doriAll,dsfAll,doriEucAll,dsfEucAll,distAll)

global PW

PW.Ddom = [0 prctile(distAll,33) prctile(distAll,66) max(distAll)];
%PW.Ddom = [0 prctile(distAll,50) max(distAll)];
PW.Ddom = round(PW.Ddom);

clear dori dsf doriEuc dsfEuc 
for i = 1:length(PW.Ddom)-1
    
    id = find(distAll>PW.Ddom(i) & distAll<=PW.Ddom(i+1) & ~isnan(doriAll) & ~isnan(dsfAll) & ~isnan(doriEucAll) & ~isnan(dsfEucAll));
    dori{i} = doriAll(id);
    dsf{i} = dsfAll(id);
    doriEuc{i} = doriEucAll(id);
    dsfEuc{i} = dsfEucAll(id);
    
end

%%  Plot 2D histogram for each distance
ND = length(doriEuc);
figure(99)
for i = 1:ND
    [mat xdom ydom] = smoothscatter(doriEuc{i},dsfEuc{i},.01,.01,[-1 1], [-1 1]);
    mat = mat/sum(mat(:));

    %matnorm = ones(length(mat(:,1)),1)*sum(mat);
    %matnorm = sum(mat,2)*ones(1,length(mat(1,:)));
    %mat = mat./matnorm;
    
    %mat = sum(mat,2)*sum(mat,1);

    %subplot(2,ND,i),scatter(doriEuc{i},dsfEuc{i},'.k')
    %xlabel('dori'), ylabel('dsf')
    
    subplot(1,ND,i),imagesc(xdom,ydom,-(mat).^.35), axis xy
    xlabel('dori'), ylabel('dsf'),
    colormap gray
    
    [r p] = corrcoef(doriEuc{i},dsfEuc{i});
    title(['P(dori,dsf|'  num2str(PW.Ddom(i))  ');   r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2))])
    
    axis square
end


%% dsf distribution given dori (Euc)

yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];

clear prc mu_dsf sig_dsf xdom
for Did = 1:length(PW.Ddom)-1  %loop through pairwise distance
    
    %doridom = [-1 prctile(doriEuc{Did},25) prctile(doriEuc{Did},50) prctile(doriEuc{Did},75) 1];
    doridom = [-1 prctile(doriEuc{Did},33) prctile(doriEuc{Did},66) 1];
    ns = length(doridom)-1;
    
    sfdom = linspace(-1,1.0001,10);

    xaxdom = doridom(1:end-1);
    
    clear Psf
    %Psf = zeros(length(sfdom),Nori);
    for i = 1:length(doridom)-1

        id = find(doriEuc{Did} >= doridom(i) & doriEuc{Did} < doridom(i+1));
        
        dum = histc(dsfEuc{Did}(id),sfdom); dum = dum(1:end-1);
        Psf(:,i) = dum;
        Psf(:,i) = Psf(:,i)/sum(Psf(:,i));
        
        idno = find(~isnan(dsfEuc{Did}(id)));
        mu_dsf(Did,i) = nanmean(dsfEuc{Did}(id));
        sig_dsf(Did,i) = nanstd(dsfEuc{Did}(id))/sqrt(length(idno));
        
        xdom{Did}(i) = nanmean(doriEuc{Did}(id));
        
        for q = 1:length(prcdom)
            prc{q}(Did,i) = prctile(dsfEuc{Did}(id),prcdom(q)*100);
        end
    end
    
    mu_dsf = mu_dsf(:,1:ns);
    cumPsf = [zeros(1,ns); cumsum(Psf(:,1:ns))];
    %cumPsf = Psf;
    dom = 0:length(cumPsf(:,1))-1;
    domI = linspace(dom(1),dom(end),length(dom)*10);
    cumPsfI = interp1(dom,cumPsf,domI);
    
%     figure(3)
%     subplot(2,3,Did)
%     plot(cumPsf,'-o'), ylim([0 1.2]), xlabel('dsfEuc') 
%     title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])
%     
%     figure(4)    
%     subplot(1,length(prcdom)+1,1)
%     imagesc(xaxdom,yaxdom,mu_dsf)
%     %set(gca,'Xtick',[]), set(gca,'Ytick',[])
%     ylabel('Cortical Distance'), xlabel('doriEuc')
%     title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori
% 
%             [dum id] = min(abs(cumPsfI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %sfreq percentile location, as a function of dOri and cortical distance
% 
%         end

%         subplot(1,length(prcdom)+1,q+1)
%         imagesc(xaxdom,yaxdom,prc{q})
%         ylabel('Cortical Distance'), xlabel('doriEuc')
%         title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(doridom)-1    
   labs{i} = [num2str(doridom(i)) ' to '  num2str(doridom(i+1))];
end

figure, bar(mu_dsf')
hold on,errorbar(mu_dsf',sig_dsf','.k')
set(gca,'XTick',1:(length(doridom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dori range (norm)'), ylabel('E(dsf) (norm)')

NDis = length(PW.Ddom)-1;
figure(99)
for i = 1:NDis
    subplot(1,3,i)
    hold on, plot(xdom{i},mu_dsf(i,:),'.-b')
    for j = 1:length(sig_dsf)
        hold on
        plot([xdom{i}(j) xdom{i}(j)],[mu_dsf(i,j)-sig_dsf(i,j) ,mu_dsf(i,j)+sig_dsf(i,j)],'b')
    end
end


%% Now dori distribution given dsf (Euc)

yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];

clear prc mu_dori N xdom
for Did = 1:length(PW.Ddom)-1  %loop through pairwise distance

    %dsfdom = [-1 prctile(dsfEuc{Did},25) prctile(dsfEuc{Did},50) prctile(dsfEuc{Did},75) 1];
    dsfdom = [-1 prctile(dsfEuc{Did},33) prctile(dsfEuc{Did},66) 1];
    ns = length(dsfdom)-1;
    
    oridom = linspace(-1,1.0001,5);

    xaxdom = dsfdom(1:end-1);
    
    clear Pori
    %Pori = zeros(length(sfdom),Nori);
    for i = 1:length(dsfdom)-1

        id = find(dsfEuc{Did} >= dsfdom(i) & dsfEuc{Did} < dsfdom(i+1));
        oriDist{Did,i} = doriEuc{Did}(id);
        
        dum = histc(oriDist{Did,i},oridom);
        Pori(:,i) = dum(1:end-1);  %last element of histc is outside of range
        N(Did,i) = sum(Pori(:,i));
        Pori(:,i) = Pori(:,i)/N(Did,i);

        idno = find(~isnan(doriEuc{Did}(id)));
        mu_dori(Did,i) = nanmean(doriEuc{Did}(id));
        sig_dori(Did,i) = nanstd(doriEuc{Did}(id))/sqrt(length(idno));
        
        xdom{Did}(i) = nanmean(dsfEuc{Did}(id));
        
        for q = 1:length(prcdom)
            prc{q}(Did,i) = prctile(oriDist{Did,i},prcdom(q)*100);
        end
    end
    
    mu_dori = mu_dori(:,1:ns);
    cumPori = [zeros(1,ns); cumsum(Pori(:,1:ns))];
    %cumPsf = Psf;
    dom = 0:length(cumPori(:,1))-1;
    domI = linspace(dom(1),dom(end),length(dom)*10);
    cumPoriI = interp1(dom,cumPori,domI);
    
%     figure(34)
%     subplot(2,3,Did)
%     plot(cumPori,'-o'), ylim([0 1.2]), xlabel('doriEuc') 
%     title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])
% 
%     figure(45)    
%     subplot(1,length(prcdom)+1,1)
%     imagesc(xaxdom,yaxdom,mu_dori)
%     %set(gca,'Xtick',[]), set(gca,'Ytick',[])
%     ylabel('Cortical Distance'), xlabel('dsfEuc')
%     title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori 
%             [dum id] = min(abs(cumPoriI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %sfreq percentile location, as a function of dOri and cortical distance 
%         end
%         subplot(1,length(prcdom)+1,q+1)
% 
%         imagesc(xaxdom,yaxdom,prc{q})
%         ylabel('Cortical Distance'), xlabel('dsfEuc')
%         title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(dsfdom)-1    
   labs{i} = [num2str(dsfdom(i)) ' to '  num2str(dsfdom(i+1))];
end
figure, bar(mu_dori')
hold on,errorbar(mu_dori',sig_dori','.k')
set(gca,'XTick',1:(length(dsfdom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dsf range (norm)'), ylabel('E(dori) (norm)')

NDis = length(PW.Ddom)-1;

figure(99)
for i = 1:NDis
    subplot(1,3,i)
    hold on, plot(mu_dori(i,:),xdom{i},'.-r')
    for j = 1:length(sig_dsf)
        hold on
        plot([mu_dori(i,j)-sig_dori(i,j) ,mu_dori(i,j)+sig_dori(i,j)],[xdom{i}(j) xdom{i}(j)],'r')
    end
end


%%%%%%%%%%%%%%%Now regular dori and dsf%%%%%%%%%%%%%


%%  Plot 2D histogram for each distance
ND = length(dori);
figure(98)
for i = 1:ND
    [mat xdom ydom] = smoothscatter(dori{i},dsf{i},.6,.03,[0 90],[0 3]);

    %matnorm = ones(length(mat(:,1)),1)*sum(mat);
    %matnorm = sum(mat,2)*ones(1,length(mat(1,:)));
    %mat = mat./matnorm;
    
    %mat = sum(mat,2)*sum(mat,1);

    %subplot(2,ND,i),scatter(dori{i},dsf{i},'.k')
    %xlabel('dori'), ylabel('dsf')
    
    subplot(1,ND,i),imagesc(xdom,ydom,-(mat).^.35), axis xy
    xlabel('dori'), ylabel('dsf'),
    set(gca,'Xtick',[0 45 90])
    colormap gray
    
    [r p] = corrcoef(dori{i},dsf{i});
    title(['P(dori,dsf|'  num2str(PW.Ddom(i))  ');   r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2))])
    
    axis square

end




%% dsf distribution given dori

yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];

clear prc mu_dsf sig_dsf xdom
for Did = 1:length(PW.Ddom)-1  %loop through pairwise cortical distance

    %doridom = [0 prctile(dori{Did},25) prctile(dori{Did},50) prctile(dori{Did},75) 90];
    doridom = [0 prctile(dori{Did},33) prctile(dori{Did},66) 90];
    ns = length(doridom)-1;
    
    sfdom = linspace(0,max(dsf{Did})+.0001,10);
    
    xaxdom = doridom(1:end-1);
    
    clear Psf
    %Psf = zeros(length(sfdom),Nori);
    for i = 1:length(doridom)-1  %loop through pairwise orientation difference

        id = find(dori{Did} >= doridom(i) & dori{Did} < doridom(i+1));
        dum = histc(dsf{Did}(id),sfdom)'; dum = dum(1:end-1);
        Psf(:,i) = dum;
        Psf(:,i) = Psf(:,i)/sum(Psf(:,i));
        
        idno = find(~isnan(dsf{Did}(id)));
        mu_dsf(Did,i) = nanmean(dsf{Did}(id));
        sig_dsf(Did,i) = nanstd(dsf{Did}(id))/sqrt(length(idno));
        
        xdom{Did}(i) = nanmean(dori{Did}(id));

        for q = 1:length(prcdom)
            prc{q}(Did,i) = prctile(dsf{Did}(id),prcdom(q)*100);
        end
        
    end
    
    mu_dsf = mu_dsf(:,1:ns);
    cumPsf = [zeros(1,ns); cumsum(Psf(:,1:ns))];
    dom = 0:length(cumPsf(:,1))-1;
    domI = linspace(dom(1),dom(end),length(dom)*10);
    cumPsfI = interp1(dom,cumPsf,domI);
    
%     figure(5)
%     subplot(2,3,Did)
%     plot([sfdom],cumPsf,'-o'), ylim([0 1.2]), xlabel('dsf'), 
%     title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])
%     
%     figure(6)    
%     subplot(1,length(prcdom)+1,1)
%     imagesc(xaxdom,yaxdom,mu_dsf)
%     ylabel('Cortical Distance'), xlabel('dori')
%     title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori
% 
%             [dum id] = min(abs(cumPsfI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %sfreq percentile location, as a function of dOri and cortical distance
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
figure,bar(mu_dsf')
hold on,errorbar(mu_dsf',sig_dsf','.k')
set(gca,'XTick',1:(length(doridom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dori range'), ylabel('E(dsf)')

NDis = length(PW.Ddom)-1;
figure(98)
for i = 1:NDis
    subplot(1,3,i)
    hold on, plot(xdom{i},mu_dsf(i,:),'.-b')
    for j = 1:length(sig_dsf)
        hold on
        plot([xdom{i}(j) xdom{i}(j)],[mu_dsf(i,j)-sig_dsf(i,j) ,mu_dsf(i,j)+sig_dsf(i,j)],'b')
    end
end


%% dori distribution given dsf


yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];
%dsfdom = [0 .7 1.4 2.1 10];

clear prc mu_dori sig_dori xdom
for Did = 1:length(PW.Ddom)-1  %loop through pairwise distance

    dsfdom = [0 prctile(dsf{Did},25) prctile(dsf{Did},50) prctile(dsf{Did},75) 10];
    dsfdom = [0 prctile(dsf{Did},33) prctile(dsf{Did},66) 10];
    ns = length(dsfdom)-1;
    
    oridom = linspace(0,90,10);

    xaxdom = doridom(1:end-1);
    
    clear Pori mu
    %Psf = zeros(length(sfdom),Nori);
    for i = 1:length(dsfdom)-1

        id = find(dsf{Did} >= dsfdom(i) & dsf{Did} < dsfdom(i+1));
        
        dum = histc(dori{Did}(id),oridom); 
        Pori(:,i) = dum(1:end-1); %last element of histc is outside of range
        Pori(:,i) = Pori(:,i)/sum(Pori(:,i));
        
        idno = find(~isnan(dori{Did}(id)));
        mu_dori(Did,i) = nanmean(dori{Did}(id));
        sig_dori(Did,i) = nanstd(dori{Did}(id))/sqrt(length(idno));
        
        xdom{Did}(i) = nanmean(dsf{Did}(id));
        
        for q = 1:length(prcdom)
            prc{q}(Did,i) = prctile(dori{Did}(id),prcdom(q)*100);
        end
    end
    
    mu_dori = mu_dori(:,1:ns);
    cumPori = [zeros(1,ns); cumsum(Pori(:,1:ns))];
    %cumPsf = Psf;
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
%     ylabel('Cortical Distance'), xlabel('dsf')
%     title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dori
% 
%             [dum id] = min(abs(cumPsfI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %sfreq percentile location, as a function of dOri and cortical distance
% 
%         end

%         subplot(1,length(prcdom)+1,q+1)
%         imagesc(xaxdom,yaxdom,prc{q})
%         ylabel('Cortical Distance'), xlabel('dsf')
%         title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(dsfdom)-1    
   labs{i} = [num2str(dsfdom(i)) ' to '  num2str(dsfdom(i+1))];
end
figure, bar(mu_dori')
hold on,errorbar(mu_dori',sig_dori','.k')
set(gca,'XTick',1:(length(dsfdom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dsf range'), ylabel('E(dori)')

NDis = length(PW.Ddom)-1;

figure(98)
for i = 1:NDis
    subplot(1,3,i)
    hold on, plot(mu_dori(i,:),xdom{i},'.-r')
    for j = 1:length(sig_dsf)
        hold on
        plot([mu_dori(i,j)-sig_dori(i,j) ,mu_dori(i,j)+sig_dori(i,j)],[xdom{i}(j) xdom{i}(j)],'r')
    end
end
%% anova
% dim = size(oriDist);
% vecAll = []; Distgroup = [];  sfgroup = [];
% k = 1;
% for i = 1:dim(1)
%     for j = 1:dim(2)
%         vecAll = [vecAll; (oriDist{i,j}(:))];
%         Distgroup = [Distgroup i*ones(1,length(oriDist{i,j}(:)))];
%         sfgroup = [sfgroup j*ones(1,length(oriDist{i,j}(:)))];
%         k = k+1;
%     end
% end
% 
% [p,tbl,stats] = anovan(vecAll,{Distgroup sfgroup})