function sfOriaxes2(dori,dsf,ax,dist,animID)

%Ian Nauhaus

animdom = unique(animID)

for j = 1:length(animdom)  %loop through each ROI

    id = find(animID == animdom(j));
    
    doriAll{j} = dori(id);
    dsfAll{j} = dsf(id);
    axAll{j} = ax(id);
    distAll{j} = dist(id);
   
    
end

%%


axedges = 0:10:180; dax = axedges(2)-axedges(1);
axW = 2*dax;
axdom = axedges(1:end-1)+dax/2;
figure,
for i = 1:length(doriAll) %loop through each ROI
    axAll{i} = angle(exp(1i*axAll{i}*2*pi/180))*180/pi/2;
    axAll{i} = axAll{i} + 90;    
    
    doriAlln = doriAll{i}./distAll{i}*1000; %degrees/mm
    dsfAlln = dsfAll{i}./distAll{i}*1000; %octaves/mm

%     doriAlln = doriAll{i}; %degrees/mm
%     dsfAlln = dsfAll{i}; %octaves/mm
    
    varflag = 1;
    f = 1;
    if varflag
        doriAlln = abs(doriAlln);
        dsfAlln = abs(dsfAlln);
        f = 2;
    end
    
    id = find(distAll{i}>10);
    distAll{i} = distAll{i}(id); axAll{i} = axAll{i}(id); doriAlln = doriAlln(id); dsfAlln = dsfAlln(id);    
    
    for j = 1:length(axdom)  %loop through axis domain
    
        id = find(axAll{i} > axdom(j)-axW/2 & axAll{i} < axdom(j)+axW/2);
        
        dori_mu(i,j) = nanmedian(doriAlln(id));
        dori_SE(i,j) = std(doriAlln(id))/sqrt(length(id));
        
        dsf_mu(i,j) = nanmedian(dsfAlln(id));
        dsf_SE(i,j) = std(dsfAlln(id))/sqrt(length(id));
        
    end
    
    subplot(2,length(doriAll),i)
    scatter(axAll{i},doriAlln,'.')
    hold on
    %errorbar(axdom,dori_mu(i,:),dori_SE(i,:),'k','LineWidth',3), ylim([prctile(doriAlln,2) prctile(doriAlln,98)])
    plot(axdom,dori_mu(i,:),'k','LineWidth',3), ylim([prctile(doriAlln,2) prctile(doriAlln,98)])
    Resori(i) = sum(dori_mu(i,:).*exp(1i*axdom*pi/180*f));  %resultant
    Magori(i) = 2*abs(Resori(i)/(length(dori_mu(i,:)))); %estimate of dori/mm (peak amplitude)
    Parori(i) = abs(Resori(i)/sum((dori_mu(i,:))));  %normalized magnitude (parallism)
    oriprincax(i) = angle(Resori(i))*180/pi/f;
    title(['ph = ' num2str(round(oriprincax(i))) ' amp = ' num2str(round(Magori(i)*10)/10)])
    ylabel('degrees/mm')
    hold off
    
    subplot(2,length(doriAll),i+length(doriAll))
    scatter(axAll{i},dsfAlln,'.')
    hold on    
%     errorbar(axdom,dsf_mu(i,:),dsf_SE(i,:),'k','LineWidth',3), ylim([prctile(dsfAlln,2) prctile(dsfAlln,98)])
    plot(axdom,dsf_mu(i,:),'k','LineWidth',3), ylim([prctile(dsfAlln,2) prctile(dsfAlln,98)])    
    Ressf(i) = sum(dsf_mu(i,:).*exp(1i*axdom*pi/180*f));
    Magsf(i) = 2*abs(Ressf(i)/(length(dsf_mu(i,:)))); %estimate of dori/mm (peak amplitude)
    Parsf(i) = abs(Ressf(i)/sum((dsf_mu(i,:))));  %normalized magnitude (parallism)
    sfprincax(i) = angle(Ressf(i))*180/pi/f;
    
    ylabel('octaves/mm')
    
    M(i) = sqrt(abs(Magori(i).*Magsf(i)));

    %daxis(i) = abs(oridiff(sfprincax(i)*pi/180,oriprincax(i)*pi/180)*180/pi);
    daxis(i) = 90-abs(abs(sfprincax(i) - oriprincax(i)) - 90);  %how it should be defined in the paper
    
    title(['ph=' num2str(round(sfprincax(i))) ' amp=' num2str(round(Magsf(i)*10)/10) ' diff=' num2str(round(daxis(i)))])
    hold off
end

figure,compass(M.*cos(daxis*pi/180),M.*sin(daxis*pi/180))

hcorn = [0:15:90]; ddom = hcorn(2)-hcorn(1); hdom = hcorn(1:end-1)+ddom/2;
dum = histc(daxis,hcorn);
figure,bar(hdom,dum(1:end-1))
set(gca,'XTick',hcorn)
figure,scatter(abs(Resori),abs(Ressf)), xlabel('deg/mm'), ylabel('octaves/mm')
function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;

