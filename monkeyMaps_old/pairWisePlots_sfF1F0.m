function pairWisePlots_sfF1F0(dsfAll,dF1F0All,distAll)

global PW

PW.Ddom = [0 prctile(distAll,33) prctile(distAll,66) max(distAll)];
%PW.Ddom = [0 prctile(distAll,50) max(distAll)];
PW.Ddom = round(PW.Ddom);

clear dF1F0 dsf
for i = 1:length(PW.Ddom)-1
    
    id = find(distAll>PW.Ddom(i) & distAll<=PW.Ddom(i+1) & ~isnan(dF1F0All) & ~isnan(dsfAll));
    dF1F0{i} = dF1F0All(id);
    dsf{i} = dsfAll(id);
    
end



%%%%%%%%%%%%%%%Now regular dF1F0 and dsf%%%%%%%%%%%%%

%%  Plot 2D histogram for each distance
ND = length(dF1F0);
figure(98)
for i = 1:ND
    [mat xdom ydom] = smoothscatter(dF1F0{i},dsf{i},.04,.04,[-3 3],[-3 3]);

    %matnorm = ones(length(mat(:,1)),1)*sum(mat);
    %matnorm = sum(mat,2)*ones(1,length(mat(1,:)));
    %mat = mat./matnorm;
    
    %mat = sum(mat,2)*sum(mat,1);

    %subplot(2,ND,i),scatter(dF1F0{i},dsf{i},'.k')
    %xlabel('dF1F0'), ylabel('dsf')
    
    subplot(1,ND,i),imagesc(xdom,ydom,-(mat).^.35), axis xy
    xlabel('dF1F0'), ylabel('dsf'),
    colormap gray
    
    [r p] = corrcoef(dF1F0{i},dsf{i});
    title(['P(dF1F0,dsf|'  num2str(PW.Ddom(i))  ');   r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2))])
    
    axis square

end




%% dsf distribution given dF1F0

yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];

clear prc mu_dsf sig_dsf xdom
for Did = 1:length(PW.Ddom)-1  %loop through pairwise cortical distance

    %dF1F0dom = [0 prctile(dF1F0{Did},25) prctile(dF1F0{Did},50) prctile(dF1F0{Did},75) 90];
    dF1F0dom = [prctile(dF1F0{Did},2) prctile(dF1F0{Did},33) prctile(dF1F0{Did},66) prctile(dF1F0{Did},98)];
    ns = length(dF1F0dom)-1;
    
    sfdom = linspace(0,max(dsf{Did})+.0001,10);
    
    xaxdom = dF1F0dom(1:end-1);
    
    clear Psf
    %Psf = zeros(length(sfdom),NF1F0);
    for i = 1:length(dF1F0dom)-1  %loop through pairwise F1F0entation difference

        id = find(dF1F0{Did} >= dF1F0dom(i) & dF1F0{Did} < dF1F0dom(i+1));
        dum = histc(dsf{Did}(id),sfdom)'; dum = dum(1:end-1);
        Psf(:,i) = dum;
        Psf(:,i) = Psf(:,i)/sum(Psf(:,i));
        
        idno = find(~isnan(dsf{Did}(id)));
        mu_dsf(Did,i) = nanmean(dsf{Did}(id));
        sig_dsf(Did,i) = nanstd(dsf{Did}(id))/sqrt(length(idno));
        
        xdom{Did}(i) = nanmean(dF1F0{Did}(id));

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
%     ylabel('Cortical Distance'), xlabel('dF1F0')
%     title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dF1F0
% 
%             [dum id] = min(abs(cumPsfI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %sfreq percentile location, as a function of dF1F0 and cortical distance
% 
%         end

%         subplot(1,length(prcdom)+1,q+1)
%         imagesc(xaxdom,yaxdom,prc{q})
%         ylabel('Cortical Distance'), xlabel('dF1F0')
%         title([num2str(prcdom(q)*100) ' percentile'])
    end

end

for i = 1:length(dF1F0dom)-1    
   labs{i} = [num2str(round(dF1F0dom(i))) ' to '  num2str(round(dF1F0dom(i+1)))];
end
figure,bar(mu_dsf')
hold on,errorbar(mu_dsf',sig_dsf','.k')
set(gca,'XTick',1:(length(dF1F0dom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dF1F0 range'), ylabel('E(dsf)')

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


%% dF1F0 distribution given dsf


yaxdom = (PW.Ddom(1:end-1) + PW.Ddom(2:end))/2;

prcdom = [.15 .5 .85];
%dsfdom = [0 .7 1.4 2.1 10];

clear prc mu_dF1F0 sig_dF1F0 xdom
for Did = 1:length(PW.Ddom)-1  %loop through pairwise distance

    %dsfdom = [0 prctile(dsf{Did},25) prctile(dsf{Did},50) prctile(dsf{Did},75) 10];
    dsfdom = [prctile(dsf{Did},2) prctile(dsf{Did},33) prctile(dsf{Did},66) prctile(dsf{Did},98)];
    ns = length(dsfdom)-1;
    
    F1F0dom = linspace(0,90,10);

    xaxdom = dF1F0dom(1:end-1);
    
    clear PF1F0 mu
    %Psf = zeros(length(sfdom),NF1F0);
    for i = 1:length(dsfdom)-1

        id = find(dsf{Did} >= dsfdom(i) & dsf{Did} < dsfdom(i+1));
        
        dum = histc(dF1F0{Did}(id),F1F0dom); 
        PF1F0(:,i) = dum(1:end-1); %last element of histc is outside of range
        PF1F0(:,i) = PF1F0(:,i)/sum(PF1F0(:,i));
        
        idno = find(~isnan(dF1F0{Did}(id)));
        mu_dF1F0(Did,i) = nanmean(dF1F0{Did}(id));
        sig_dF1F0(Did,i) = nanstd(dF1F0{Did}(id))/sqrt(length(idno));
        
        xdom{Did}(i) = nanmean(dsf{Did}(id));
        
        for q = 1:length(prcdom)
            prc{q}(Did,i) = prctile(dF1F0{Did}(id),prcdom(q)*100);
        end
    end
    
    mu_dF1F0 = mu_dF1F0(:,1:ns);
    cumPF1F0 = [zeros(1,ns); cumsum(PF1F0(:,1:ns))];
    %cumPsf = Psf;
    dom = 0:length(cumPF1F0(:,1))-1;
    domI = linspace(dom(1),dom(end),length(dom)*10);
    cumPF1F0I = interp1(dom,cumPF1F0,domI);
    
%     figure(38)
%     subplot(2,3,Did)
%     plot(cumPF1F0,'-o'), ylim([0 1.2]), xlabel('dF1F0') 
%     title([num2str(PW.Ddom(Did)) ' to ' num2str(PW.Ddom(Did+1)) ' microns'])
%     
%     figure(39)    
%     subplot(1,length(prcdom)+1,1)
%     imagesc(xaxdom,yaxdom,mu_dF1F0)
%     ylabel('Cortical Distance'), xlabel('dsf')
%     title(['mean'])
    
    for q = 1:length(prcdom)
%         for k = 1:ns %loop through dF1F0
% 
%             [dum id] = min(abs(cumPsfI(:,k)-prcdom(q)));
%             prc{q}(Did,k) = id; %sfreq percentile location, as a function of dF1F0 and cortical distance
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
figure, bar(mu_dF1F0')
hold on,errorbar(mu_dF1F0',sig_dF1F0','.k')
set(gca,'XTick',1:(length(dsfdom)-1) )
set(gca,'XTickLabel',labs)
xlabel('dsf range'), ylabel('E(dF1F0)')

NDis = length(PW.Ddom)-1;

figure(98)
for i = 1:NDis
    subplot(1,3,i)
    hold on, plot(mu_dF1F0(i,:),xdom{i},'.-r')
    for j = 1:length(sig_dsf)
        hold on
        plot([mu_dF1F0(i,j)-sig_dF1F0(i,j) ,mu_dF1F0(i,j)+sig_dF1F0(i,j)],[xdom{i}(j) xdom{i}(j)],'r')
    end
end
%% anova
% dim = size(F1F0Dist);
% vecAll = []; Distgroup = [];  sfgroup = [];
% k = 1;
% for i = 1:dim(1)
%     for j = 1:dim(2)
%         vecAll = [vecAll; (F1F0Dist{i,j}(:))];
%         Distgroup = [Distgroup i*ones(1,length(F1F0Dist{i,j}(:)))];
%         sfgroup = [sfgroup j*ones(1,length(F1F0Dist{i,j}(:)))];
%         k = k+1;
%     end
% end
% 
% [p,tbl,stats] = anovan(vecAll,{Distgroup sfgroup})