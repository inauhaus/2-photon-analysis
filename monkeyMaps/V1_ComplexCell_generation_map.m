function [XsigVsf BWLinVsf BWVsf sfpref F1vsBand] = V1_ComplexCell_generation_map(sfBands,xsigDom,sfsigDom,sf2sigModel)

%sfBands: Vector; Center of pooling functions in sfreq dimension: cyc/deg
%xsigDom: Vector; Pooling width (sigma) in spatial dimension; degrees
%sfsigDom: Vector; Pooling width (sigma) in sfreq dimension: cyc/degS

%Outputs all have the same dimensions;
%Dims 1,2,3 are sfBands, xsigDom, sfsigDom, respectively

%sfBands = logspace(log10(.5),log10(8),5);

% sfBands = linspace(.5,8,5);
% 
% xsigDom = [.25]+.001;
% sfsigDom = [3]+.001;

%%
%First create sampling distribution of simple cell inputs
%spatial frequency and position are the two domains

%Create T(x-xp,sf): the population of V1 simple receptive fields
%xp is position shift.

xN = 300;
sfN = 300;
MF = 1;  %mm/deg
cortex_mm = linspace(-3,3,xN); %location of each neuron in mm
xCenterdom = cortex_mm/MF; %RF location
%sfprefdom = log2(logspace(log10(.1),log10(20),sfN)); 
mapPeriod = 0.7;  %mm
sfprefdom = sin(2*pi*(cortex_mm/mapPeriod))*7.5/2 + 7.5/2 +.5;

%sfprefdom = (linspace((.1),(20),sfN));
xdom = linspace(-3,3,xN);



[xCentermesh xmesh] = meshgrid(xCenterdom,xdom); %x y

[sfprefmesh xmesh] = meshgrid(sfprefdom,xdom); %x y


sigmesh = getSizefromSFpref(sf2sigModel,sfprefmesh);

%sigmesh = -sfprefmesh - log2(4);

%Derivative of Gaussian model
% SimplePop = exp((-(xmesh-xCentermesh).^2)./(2*sigmesh.^2));
% SimplePop = diff(SimplePop,[],3);
% SimplePop2 = diff(SimplePop,[],3);
% 
% SimplePop(:,:,end+1) = SimplePop(:,:,end);
% SimplePop2(:,:,end+1) = SimplePop2(:,:,end);
% SimplePop2(:,:,end+1) = SimplePop2(:,:,end);

%Gabor model
Env = exp((-(xmesh-xCentermesh).^2)./(2*sigmesh.^2));
%SimplePop = Env.*cos(xmesh.*sfprefmesh*2*pi);
%SimplePop2 = Env.*sin(xmesh.*sfprefmesh*2*pi);

SimplePop = Env.*cos((xmesh-xCentermesh).*sfprefmesh*2*pi);
SimplePop2 = Env.*sin((xmesh-xCentermesh).*sfprefmesh*2*pi);

ma = max(abs(SimplePop),[],1);
ma2 = max(abs(SimplePop2),[],1);
ma = sqrt(sum(SimplePop.*SimplePop));
ma2 = sqrt(sum(SimplePop2.*SimplePop2));
for i = 1:size(SimplePop,2)
    
    SimplePop(:,i) = SimplePop(:,i)./ma(i);  
    SimplePop2(:,i) = SimplePop2(:,i)./ma(i);
    
end

figure
for i = 1:50:xN
   
    %dum = SimplePop(i,100,:);
    %dum2 = SimplePop2(i,100,:);
    %Env1 = squeeze(sqrt(dum.^2 + dum2.^2));
    
    Env1 = squeeze(Env(:,i));
    plot(xdom,Env1,'k')
    hold on
end
    
    

%%
%xsigDom = .2; %Gaussian pooling function in degrees of visual space
%sfsigDom = .01;  %Gaussian pooling octave

clear BWVsf BWLinVsf XsigVsf sfpref F1vsBand

xPosDom = linspace(-mapPeriod/2,mapPeriod/2,10);

cdom = jet;
id = round(linspace(1,length(cdom(:,1)),length(xPosDom)));
cdom = cdom(id,:);

figure(88)
clf

figure(87)
clf

for i = 1:length(xPosDom)
    
    %octaves
    %Gweight = exp(-xCentermesh.^2/(2*xsig^2)) .* exp(-(log2(sfprefmesh)-log2(sfBands(i))).^2/(2*sfsig^2));
    
    %cyc/deg
    xsig = xsigDom(1);
    Gweight = exp(-(xCentermesh-xPosDom(i)).^2/(2*xsig^2));
    
    %imagesc(xCenterdom,sfprefdom,Gweight(:,:,50)),
    figure(87)
    plot(xCenterdom,sfprefdom,'Color',cdom(i,:)),
    xlim([-.5 .5]), ylim([.1 7]), xlabel('degrees'), ylabel('cyc/deg')
    hold on
    %dum = get(gca,'YTick');
    %set(gca,'YTickLabel',2.^dum);
    
%     [dum sfexamp(1)] = min(abs(sfprefdom-sfBands(2)));
%     [dum sfexamp(2)] = min(abs(sfprefdom-sfBands(end-1)));
%     [dum xexamp(1)] = min(abs(xCenterdom-(-xsigDom)));
%     [dum xexamp(2)] = min(abs(xCenterdom-(xsigDom)));
    
%     plot(xCenterdom(xexamp),sfprefdom(sfexamp),'.k','MarkerSize',15)
%     hold on
%     plot(xCenterdom(fliplr(xexamp)),sfprefdom(sfexamp),'.k','MarkerSize',15)
%     
%     figure(86)
%     subplot(1,2,1)
%     plot(xCenterdom,squeeze(SimplePop2(sfexamp(1),xexamp(1),:)),'k')
%     xlim([-1 1])
%     hold on
%     plot(xCenterdom,squeeze(SimplePop2(sfexamp(1),xexamp(2),:)),'Color',[.5 .5 .5])
%     xlim([-1 1])
%     
%     subplot(1,2,2)
%     plot(xCenterdom,squeeze(SimplePop2(sfexamp(2),xexamp(1),:)),'k')
%     hold on
%     plot(xCenterdom,squeeze(SimplePop2(sfexamp(2),xexamp(2),:)),'Color',[.5 .5 .5])
    
    In = ((Gweight.*SimplePop).^2 + (Gweight.*SimplePop2).^2);
    %In = (Gweight.*SimplePop).^2;
    %%%%
    
    %     dum = Gweight(:,:,50);
    %     [idy idx] = find(dum == max(dum(:)));
    %     RFout = squeeze(Env(idy(1),idx(1),:));
    %     figure, plot(RFout)
    %%%%
    
    
    RFout = sqrt(squeeze(sum(In,2)));
    RFout = RFout/max(RFout);
    
    %Get sf tuning curve for each cell (F0)
    TCin = abs(fft(Gweight.*SimplePop2,[],1)); %Mag of response of each neuron to each spatial frequency. e.g. Number of spikes it passes downstream.
    
    F0out = squeeze(sum(TCin,2)); %Now add up the spikes for each spatial frequency
    F0out = F0out/pi; %half wave rectification
    
    %Get F1 modulation for each cell
    F1in = fft(Gweight.*SimplePop2,[],1);
    F1out = abs(squeeze(sum(F1in,2))); %add up the vectors. High sfs will interfere deconstructively... low sfs will be construct.
    F1out = F1out/2; %half wave rectification
    %modDepth = F1out./F0out;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%In a half-wave rectification of the sine wave, A*sin(x), the F1 harmonic has an
    %%amplitude of A/2, and the DC is A/pi. So, max "F1/F0" is pi.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    F0dom = linspace(0,length(xdom)/range(xdom),length(F0out));
    F0out = F0out(1:end/2);
    F0dom = F0dom(1:end/2);
    
    %Normalize
    [F0ma idpref] = max(F0out);
    F0out = F0out/F0ma;
    
    %F1vsBand(i) = F1out(idpref)*sfBands(i); %Multiply by the sf to capture carrier amplitude. Otherwise, peak energy falls from fewer sample points in the envelope.
    F1vsBand(i) = F1out(idpref)/F0ma*2; %Mod ratio [0 2]
    %F1vsBand(i) = max(F1out);
    
    figure(88)
    subplot(3,2,1), plot(xdom,RFout,'Color',cdom(i,:)), hold on
    plot([xdom(1) xdom(end)],[.61 .61],'--','Color',[.5 .5 .5])
    xlim([-2 2])
    xlabel('deg')
    subplot(3,2,2), plot(F0dom,F0out,'Color',cdom(i,:)), hold on
    plot([F0dom(1) F0dom(end)],[.61 .61],'--','Color',[.5 .5 .5])
    xlim([0 15])
    xlabel('c/deg')
    
    %Get RF size (sig)
    xdomI = linspace(xdom(1),xdom(end),length(xdom)*30000);
    RFoutI = interp1(xdom,RFout,xdomI,'spline');
    [dum idma] = max(RFoutI);
    [dum idsig] = min(abs(RFoutI(idma:end)-0.6065));
    XsigVsf(i) = idsig*(xdomI(end)-xdomI(end-1));
    
    %Get Sfreq bandwidth (sig)
    TCdomI = linspace(F0dom(1),F0dom(end),length(F0dom)*30000);
    TCoutI = interp1(F0dom,F0out,TCdomI,'spline');
    
    [dum idma] = max(TCoutI);
    idma = idma(1);
    sfpref(i) = TCdomI(idma);
    
    [sflo sfhi] = gethcolco(TCdomI,TCoutI,0.61);
    
    
    %             [dum idsigHi] = min(abs(TCoutI(idma:end)-.6065));
    %             idsigHi = idsigHi + idma - 1;
    %             [dum idsigLow] = min(abs(TCoutI(1:idma)-.6065));
    %             dsf = TCdomI(end)-TCdomI(end-1);
    BWLinVsf(i) = (sfhi-sflo)/2; %sigma
    BWLinVsf(i) = (sfhi-sfpref(i)); %sigma
    BWVsf(i) = (sfhi/sflo);
    
end


figure, plot(sfBands,F1vsBand,'-o')
xlabel('Center of SF Band')
ylabel('F1 amplitude')

cdom = [1 0 0; 0 1 0; 0 0 1];
%BW (2sig) from poster is ~1.4 c/deg
figure(88)
for k = 1:length(sfsigDom)
    subplot(3,2,3)
    semilogx(sfpref(:,:,k),XsigVsf(:,:,k),'.-k')
    set(gca,'XTick',round(sfpref(:,1)*100)/100)
    xlabel('peak sf preference'), ylabel('RF size (sigma)')
    ylim([0 max(XsigVsf(:))+.1])
    hold on
    
    subplot(3,2,4)
    loglog(sfpref(:,:,k),BWLinVsf(:,:,k),'.-k')
    set(gca,'XTick',round(sfpref(:,1)*100)/100)
    set(gca,'YTick',[.25 .5 1 2 4])
    xlabel('peak sf preference'), ylabel('SF delta bandwidth (sigma)')
    ylim([0 max(BWLinVsf(:))+.1])
    hold on
    
    subplot(3,2,6)
    loglog(sfpref(:,:,k),(BWVsf(:,:,k)),'.-k')
    set(gca,'XTick',round(sfpref(:,1)*100)/100)
    xlabel('peak sf preference'), ylabel('SF ratio bandwidth (sigma)')
    ylim([0 max(BWVsf(:))+.1])
    hold on
end

%%

%Each dot is a different pool from the T(sf,x-x') population. Each pool
%generates a complex cell RF with a corresponding spatial and sf tuning curve.  
%The y-axis is the measured width of the RF (2sigma). The x axis is the
%prediction from the sf tuning curve.  Each hue is a different window size
%of T, in the x dimension.  Each intensity is a different window size of T
%in the sf dimension.

%sig(x) = 1/(2*pi*sig(w))

figure

for i = 1:length(sfsigDom)
    
    cdom = [1 0 0; 0 1 0; 0 0 1];
    subplot(1,2,1)
    loglog(1./(sfpref(:,j,i))/2,2*XsigVsf(:,j,i),'o-','LineWidth',2,'MarkerSize',2,'Color',cdom(j,:)/i)
    set(gca,'XTick',[1/16 1/8 .25 .5 1 2],'YTick',[1/16 1/8 .25 .5 1 2])
    ylabel('actual RF width')
    xlabel('predicted RF width (2sig): .5/(sfpref)')
    xlim([.05 1.5]), ylim([.05 1.5])
    hold on
    
    subplot(1,2,2)
    loglog(2./(2*pi*BWLinVsf(:,j,i)),2*XsigVsf(:,j,i),'o-','LineWidth',2,'MarkerSize',2,'Color',cdom(j,:)/i)
    set(gca,'XTick',[1/16 1/8 .25 .5 1 2],'YTick',[1/16 1/8 .25 .5 1 2])
    ylabel('actual RF width')
    xlabel('predicted RF width (2sig): 2/(2*pi*BWsig)')
    xlim([.05 1.5]), ylim([.05 1.5])
    hold on
    
    clear leg
    for i = 1:length(xsigDom)
        leg{i} = num2str(xsigDom(i))
    end
    legend(leg)
end

subplot(1,2,1)
loglog([.05 2],[.05 2],'k'), axis square
subplot(1,2,2)
loglog([.05 2],[.05 2],'k'), axis square





%%
sigR = .9;
RFout = exp(-xdom.^2/(2*sigR^2));

F0out = abs(fft(RFout));
F0dom = linspace(0,length(xdom)/range(xdom),length(F0out));

xdomI = linspace(xdom(1),xdom(end),length(xdom)*30000);
RFoutI = interp1(xdom,RFout,xdomI,'spline');
[dum idma] = max(RFoutI);
[dum idsig] = min(abs(RFoutI(idma:end)-0.6065));
sigRF = idsig*(xdomI(end)-xdomI(end-1))

F0out = fftshift(F0out)/max(F0out);
TCdomI = linspace(F0dom(1),F0dom(end),length(F0dom)*30000);
TCoutI = interp1(F0dom,F0out,TCdomI,'spline');
TCoutI = TCoutI;
[dum idma] = max(TCoutI);
[dum idsigHi] = min(abs(TCoutI(idma:end)-.6065));
sigTC = idsigHi*(TCdomI(end)-TCdomI(end-1))

1/(2*pi*sigTC)
