function sfOriaxes(dori,dsf,ax,dist)

%Ian Nauhaus
%%
for j = 1:length(dori)  %loop through each ROI
    doriAll{j} = [];
    dsfAll{j} = [];
    axAll{j} = [];
    distAll{j} = [];
    for i = 1:length(dori{j}) %loop through each distance
        doriAll{j} = [doriAll{j}(:); dori{j}{i}(:)];
        dsfAll{j} = [dsfAll{j}(:); dsf{j}{i}(:)];
        axAll{j} = [axAll{j}(:); ax{j}{i}(:)];
        distAll{j} = [distAll{j}(:); dist{j}{i}(:)];
    end
    
    for i = 1:length(dori{j}) %loop through each distance
        doriAll{j} = [doriAll{j}(:); -dori{j}{i}(:)];
        dsfAll{j} = [dsfAll{j}(:); -dsf{j}{i}(:)];
        axAll{j} = [axAll{j}(:); ax{j}{i}(:)+180];
        distAll{j} = [distAll{j}(:); dist{j}{i}(:)];
    end
end


varflag = 1;  %if we want it to be non-directional (axis), set to 1
f = 1+varflag;

axmax = 90*(2-varflag);
axedges = -axmax:10:axmax; dax = axedges(2)-axedges(1);
axW = 4*dax;
axdom = axedges(1:end-1)+dax/2;
daxis = zeros(1,length(doriAll))*NaN;
figure,
for i = 1:length(doriAll) %loop through each ROI
%for i = 1:1
    axAll{i} = angle(exp(1i*axAll{i}*pi/180*f))*180/pi/f;
    
    doriAlln = doriAll{i}./distAll{i}*1000; %degrees/mm
    dsfAlln = dsfAll{i}./distAll{i}*1000; %octaves/mm
    
    if varflag
        doriAlln = abs(doriAlln);
        dsfAlln = abs(dsfAlln);
    end
    
    id = find(distAll{i}>0 & distAll{i}<150);
    distAll{i} = distAll{i}(id); axAll{i} = axAll{i}(id); doriAlln = doriAlln(id); dsfAlln = dsfAlln(id);    
    
    for j = 1:length(axdom)  %loop through axis domain
    
        id = find(axAll{i} > axdom(j)-axW/2 & axAll{i} < axdom(j)+axW/2);
        
        dori_mu(i,j) = trimmean(doriAlln(id),20);
        dori_SE(i,j) = std(doriAlln(id))/sqrt(length(id));
        
        dsf_mu(i,j) = trimmean(dsfAlln(id),20);
        dsf_SE(i,j) = std(dsfAlln(id))/sqrt(length(id));
        
    end
    
    subplot(2,length(doriAll),i)
    scatter(axAll{i},doriAlln,'.')
    hold on
    %errorbar(axdom,dori_mu(i,:),dori_SE(i,:),'k','LineWidth',3), ylim([prctile(doriAlln,2) prctile(doriAlln,98)])
    plot(axdom,dori_mu(i,:),'k','LineWidth',3), ylim([prctile(doriAlln,2) prctile(doriAlln,98)])
    Resori(i) = sum(dori_mu(i,:).*exp(1i*axdom*pi/180*f));  %resultant
    %Parori(i) = abs(Resori(i)/(.5*length(dori_mu(i,:)))); %estimate of dori/mm
    Parori(i) = abs(Resori(i)/sum(abs(dori_mu(i,:))));  %normalized magnitude (parallism)
    oriprincax(i) = angle(Resori(i))*180/pi/f;
    title(['ph = ' num2str(round(oriprincax(i))) ' amp = ' num2str(round(Parori(i)*100)/100)])
    ylabel('degrees/mm')
    %xlim([-90 90])
    hold off
    
    subplot(2,length(doriAll),i+length(doriAll))
    scatter(axAll{i},dsfAlln,'.')
    hold on    
%     errorbar(axdom,dsf_mu(i,:),dsf_SE(i,:),'k','LineWidth',3), ylim([prctile(dsfAlln,2) prctile(dsfAlln,98)])
    plot(axdom,dsf_mu(i,:),'k','LineWidth',3), ylim([prctile(dsfAlln,2) prctile(dsfAlln,98)])    
    Ressf(i) = sum(dsf_mu(i,:).*exp(1i*axdom*pi/180*f));
    %Parsf(i) = abs(Ressf(i)/(.5*length(dsf_mu(i,:)))); %estimate of dori/mm
    Parsf(i) = abs(Ressf(i)/sum(abs(dsf_mu(i,:))));  %normalized magnitude (parallism)
    sfprincax(i) = angle(Ressf(i))*180/pi/f;
    %xlim([-90 90])
    ylabel('octaves/mm')
    
    M(i) = sqrt(abs(Parori(i).*Parsf(i)));
    
    if varflag
        daxis(i) = abs(oridiff(sfprincax(i)*pi/180,oriprincax(i)*pi/180)*180/pi);
    else
        daxis(i) = angle(exp(1i*oriprincax(i)*pi/180) .* exp(-1i*sfprincax(i)*pi/180))*180/pi;
    end
    
    title(['ph=' num2str(round(sfprincax(i))) ' amp=' num2str(round(Parsf(i)*100)/100) ' diff=' num2str(round(daxis(i)))])
    hold off
end

figure,compass(M.*cos(daxis*pi/180),M.*sin(daxis*pi/180))

hcorn = [0:15:90]; ddom = hcorn(2)-hcorn(1); hdom = hcorn(1:end-1)+ddom/2;
dum = histc(daxis,hcorn);
figure,bar(hdom,dum(1:end-1))

figure,scatter(abs(Resori),abs(Ressf)), xlabel('deg/mm'), ylabel('octaves/mm')
function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;

