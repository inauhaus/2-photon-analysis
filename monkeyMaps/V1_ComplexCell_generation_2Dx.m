function [XsigVsf BWLinVsf sfpref F1vsBand orisigVsf] = V1_ComplexCell_generation_2D(sfBands,xsig,sfsig,sf2sigModel)

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

xCN = 100; %number of samples in spatial shift domain
sfN = 101; %number of samples in spatial frequency shift domain
xN = 100; %number of samples in RF spatial domain. Gabors are made here.

xCenterdom = linspace(-.5,.5,xCN);

%sfprefdom = log2(logspace(log10(.1),log10(20),sfN)); 
sfprefdom = (linspace(-15,15,sfN));

xdom = linspace(-3,3,xN);
ydom = xdom;


%Fourier energy domain in Cartesian coordinates
[xkmesh ykmesh] = meshgrid(sfprefdom,sfprefdom);
 
%xCentermesh and yCentermesh are retinotopy; center of RF.
[xCentermesh yCentermesh] = meshgrid(xCenterdom,xCenterdom);
yCentermesh = flipud(yCentermesh);

%Set spatial domain of each RF
[xmesh ymesh] = meshgrid(xdom,ydom);
ymesh = flipud(ymesh);

sfprefmesh = sqrt(xkmesh.^2 + ykmesh.^2);

orimesh = atan2(ykmesh,xkmesh)*180/pi;
id = find(orimesh<0);
orimesh(id) = orimesh(id)+180;

sigmesh = getSizefromSFpref(sf2sigModel,sfprefmesh);


% ma = max(abs(SimplePop),[],3);
% %ma2 = max(abs(SimplePop2),[],3);
% for i = 1:size(SimplePop,3)
%
%     SimplePop(:,:,i) = SimplePop(:,:,i)./ma;
%     %SimplePop2(:,:,i) = SimplePop2(:,:,i)./ma;
%
% end

% figure
% for i = 1:20:xN
%
%     %dum = SimplePop(i,100,:);
%     %dum2 = SimplePop2(i,100,:);
%     %Env1 = squeeze(sqrt(dum.^2 + dum2.^2));
%
%     Env1 = squeeze(Env(i,150,:));
%     plot(xdom,Env1,'k')
%     hold on
% end



%%
%xsigDom = .2; %Gaussian pooling function in degrees of visual space
%sfsigDom = .01;  %Gaussian pooling octave

clear BWVsf BWLinVsf XsigVsf sfpref F1vsBand

figure(88)
clf

sfsigx = sfsig; %SF pooling bandwidth
sfsigy = sfsig/10;

xsigx = xsig; %Spatial pooling width
xsigy = xsig/50;

AR = 1.7;
sigmeshx = sigmesh; %Receptive field size
sigmeshy = sigmesh;



GxyWeight = exp(-xCentermesh.^2/(2*xsigx^2)) .* exp(-yCentermesh.^2/(2*xsigy^2));
idGxyW = find(GxyWeight>.1*max(GxyWeight(:)));

cdom = ['b' 'g' 'r' 'c' 'k'];

figN = round(rand(1)*100);



for i = 1:length(sfBands)
    
    RFout = 0;
    F0out = 0;
    F1out = 0;
    
    GsfWeight = exp(-(abs(xkmesh)-sfBands(i)).^2/(2*sfsigx^2)) .* exp(-(abs(ykmesh)).^2/(2*sfsigy^2));
    idGsfW = find(GsfWeight>.1*max(GsfWeight(:)) & xkmesh >= 0);
    
    %figure,imagesc(GsfWeight)
    
    for jdum = 1:length(idGsfW) %loop through each location of the Fourier plane
        j = idGsfW(jdum);
        %%
        for kdum = 1:length(idGxyW) %loop through each spatial (retinotopic) shift
            k = idGxyW(kdum);
            
            %rmesh = sqrt((xmesh-xCentermesh(k)).^2 + (ymesh-yCentermesh(k)).^2);
            %Env = exp((-rmesh.^2)./(2*sigmesh(j).^2)); %envelope of Gabor            
            
            xshift = xCentermesh(k);
            yshift = yCentermesh(k);
            
            Env = exp((-(xmesh-xshift).^2)./(2*sigmeshx(j).^2)); %envelope of Gabor
            Env = Env.*exp((-(ymesh-yshift).^2)./(2*sigmeshy(j).^2)); %no rotation of Gabor???
            
            ori = orimesh(j);
            %ori = 0;
            xrotmesh = xmesh.*cos(ori*pi/180) + ymesh.*sin(ori*pi/180);% counter-clockwise rotation
            xrotmesh = xrotmesh - xshift*cos(ori*pi/180) - yshift*sin(ori*pi/180);
            %xrotmesh = (xmesh-xCentermesh(k)).*cos(ori*pi/180) + (ymesh-yCentermesh(k)).*sin(ori*pi/180);
            
            
            %xrotmesh = xmesh.*cos(ori*pi/180) + ymesh.*sin(ori*pi/180);% counter-clockwise rotation
            %xrotmesh = xrotmesh - xshift*cos(ori*pi/180) - yshift*sin(ori*pi/180);
            %xrotmesh = (xmesh-xCentermesh(k)).*cos(ori*pi/180) + (ymesh-yCentermesh(k)).*sin(ori*pi/180);
            
            SimpleRF = Env.*cos(xrotmesh.*sfprefmesh(j)*2*pi);
            SimpleRF2 = Env.*sin(xrotmesh.*sfprefmesh(j)*2*pi);
            
%             figure(22),
%             clf,
%             pause(.01),
%             imagesc(SimpleRF2)
              
            %ma = max(abs([SimpleRF(:); SimpleRF2(:)]));

            
%             ma = norm(SimpleRF);     
%             ma2 = norm(SimpleRF2); 
%             SimpleRF = SimpleRF/ma;             
%             SimpleRF2 = SimpleRF2/ma2;
%             if ma == 0 | ma2 == 0
%                 SimpleRF = zeros(size(SimpleRF));
%                 SimpleRF2 = zeros(size(SimpleRF2));
%             end
          

            
            
            %octaves
            %Gweight = exp(-xCentermesh.^2/(2*xsig^2)) .* exp(-(log2(sfprefmesh)-log2(sfBands(i))).^2/(2*sfsig^2));

            
            
            %             if j == 40
            %             figure(22),
            %             subplot(16,16,k)
            %             imagesc(SimpleRF)
            %             drawnow
            %             end
            
            %imagesc(xCenterdom,sfprefdom,Gweight(:,:,50)),

            %dum = get(gca,'YTick');
            %set(gca,'YTickLabel',2.^dum);
            
            
            %In = (Gweight.*SimplePop).^2;
            %%%%
            
            %     dum = Gweight(:,:,50);
            %     [idy idx] = find(dum == max(dum(:)));
            %     RFout = squeeze(Env(idy(1),idx(1),:));
            %     figure, plot(RFout)
            %%%%
            
            Gweight =  GsfWeight(j) * GxyWeight(k);
            RFdum = (SimpleRF.^2 + SimpleRF2.^2)*Gweight^2; %They will get sq rooted after adding them up at the end of this for loop
            RFout = RFout + RFdum;
            
            %F0outdum = abs(fft(fft(SimpleRF2,[],1),[],2)); 
            F0outdum = abs(fft(fft(SimpleRF2,[],1),[],2)) + abs(fft(fft(SimpleRF,[],1),[],2));
            
           
            F0ma = max(F0outdum(:));
            if F0ma ~= 0
                F0outdum = F0outdum/F0ma*Gweight;
            end
            F0out = F0out + F0outdum;
            
            F1outdum = fft(fft(SimpleRF2,[],1),[],2);
            if F0ma ~= 0
                F1outdum = F1outdum/F0ma*Gweight;
            end
            F1out = F1out + F1outdum;

            
            %figure(19),imagesc(fftshift(fftshift(F0out,1),2))
            %pause(.001), drawnow
            
        end
        %%
    end

    figure,imagesc(SimpleRF2)
    
    RFout = sqrt(RFout);
    RFout = RFout/max(RFout(:));
    
    %Get Fourier matrix tuning curve for each cell (F0)
    F0out = F0out/pi; %F0 after half wave rectification of 0 DC sinewave
    
    %Get F1 modulation for each cell
    F1out = abs(F1out); %Mag of resultant vector at each SF. High sfs will interfere deconstructively... low sfs will be constructive.
    F1out = F1out/2; %F1 after half wave rectification of a sinewave
    F1F0 = F1out./F0out*2; %Mod ratio [0 2]
    
    
    %Normalize F0.  Need to do this after computing F1/F0
    F0ma = max(F0out(:));
    [idprefy idprefx] = find(F0out == F0ma); idprefy = idprefy(1); idprefx = idprefx(1);
    F0out = F0out/F0ma;
    
    %F1vsBand(i) = F1out(idpref)*sfBands(i); %Multiply by the sf to capture carrier amplitude. Otherwise, peak energy falls from fewer sample points in the envelope.
    F1vsBand(i) = F1F0(idprefy,idprefx); %Mod ratio [0 2]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%In a half-wave rectification of the sine wave, A*sin(x), the F1 harmonic has an
    %%amplitude of A/2, and the DC is A/pi. So, max "F1/F0" is pi.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Fourier dom of the Receptive field. Not the "Fourier shift space"
    F0dom = linspace(0,length(xdom)/range(xdom),length(F0out(1,:))); %0 to samplerate
    
    F0FourierDom = fftshift(fftshift(F0out,1),2);
    pixperdeg = 1/(xdom(2)-xdom(1));
    [sfpref(i) BWLinVsf(i) orisigVsf(i) varacc_sf varacc_ori oritc orifit sftc sfdom sffit sffitdom] = sforiFitfrom2D_2(F0FourierDom,pixperdeg,0,'LinGauss'); %BWLin is sigma
    
    [param oriffit varacc sigma] = oriFitfrom2D(F0FourierDom);
    orisigVsf(i) = sigma;
  
    %Get Sfreq bandwidth (sig)
%    F0dom = linspace(0,length(xdom)/range(xdom),length(F0out(1,:))); %0 to samplerate
%     idtrunc = 1:length(F0dom(1,:))/2;
%     F0outx = mean(F0out);
%     F0outx = F0outx/max(F0outx);
%     F0outx = F0outx(idtrunc);
%     F0dom = F0dom(idtrunc);

    F0dom = sfdom;
    F0outx = sftc;
    TCdomI = linspace(F0dom(1),F0dom(end),length(F0dom)*30000);
    TCoutI = interp1(F0dom,F0outx,TCdomI,'spline');
    
    [dum idma] = max(TCoutI);
    idma = idma(1);
    sfpref(i) = TCdomI(idma);
    
    [sflo sfhi] = gethcolco(TCdomI,TCoutI,0.61);
    
    BWLinVsf(i) = (sfhi-sflo)/2; %sigma
    BWLinVsf(i) = (sfhi-sfpref(i)); %sigma
    BWVsf(i) = (sfhi/sflo);
    
%      [param AllffitF AllvaraccountF] = Gaussfit2Dpolar(F0FourierDom);
%      orisigVsf(i) = param(2);
%     BWLinVsf(i) = param(4);
%      sfpref(i) = param(3);

    
    F1FourierDom = fftshift(fftshift(F1out,1),2);

    figure(figN)
    subplot(2,3,1), contour(xdom,ydom,RFout,[.5 .5],cdom(i)),
    title('RF envelope')
    hold on
    subplot(2,3,2), contour(F0FourierDom,[.5 .5],cdom(i)),
    title('Fourier Energy')
    hold on
    subplot(2,3,3), contour(F1FourierDom,[.5 .5],cdom(i)),
    title('F1 Amp')
    hold on
    
    subplot(2,3,4),
    %plot(oritc/max(oritc),['.-' cdom(i)])
    hold on
    %hold on, plot(ffit)
    xlabel('orientation')


    
    %F1vsBand(i) = max(F1out);
    
    figure(88)
    RFoutx = mean(RFout); RFoutx = RFoutx/max(RFoutx);
    subplot(3,2,1), plot(xdom,RFoutx), hold on
    plot([xdom(1) xdom(end)],[.61 .61],'--','Color',[.5 .5 .5])
    xlim([-2 2])
    xlabel('deg')
    
    F0outx = F0outx/max(F0outx);
    subplot(3,2,2), plot(F0dom,F0outx,'b'), hold on
    plot([F0dom(1,1) F0dom(1,end)],[.61 .61],'--','Color',[.5 .5 .5])
    xlim([0 F0dom(end)])
    xlabel('c/deg')
    
    %Get RF size (sig)
    xdomI = linspace(xdom(1),xdom(end),length(xdom)*30000);
    RFoutI = interp1(xdom,RFoutx,xdomI,'spline');
    [dum idma] = max(RFoutI);
    [dum idsig] = min(abs(RFoutI(idma:end)-0.6065));
    XsigVsf(i) = idsig*(xdomI(end)-xdomI(end-1));
    
    
    
    


    
end

figure(figN)
subplot(2,3,5)
plot(sfBands,orisigVsf,'.-k')
xlabel('sfBand'), ylabel('ori bandwidth')
ylim([0 45])

figure, plot(sfBands,F1vsBand,'-o')
xlabel('Center of SF Band')
ylabel('F1 amplitude')
%%

cdom = [1 0 0; 0 1 0; 0 0 1];
%BW (2sig) from poster is ~1.4 c/deg
figure(88)

subplot(3,2,3)
semilogx(sfpref,XsigVsf,'.-k')
set(gca,'XTick',round(sfpref(:,1)*100)/100)
xlabel('peak sf preference'), ylabel('RF size (sigma)')
ylim([0 max(XsigVsf(:))+.1])
hold on

subplot(3,2,4)
loglog(sfpref,BWLinVsf,'.-k')
set(gca,'XTick',round(sfpref(:,1)*100)/100)
set(gca,'YTick',[.25 .5 1 2 4])
xlabel('peak sf preference'), ylabel('SF delta bandwidth (sigma)')
ylim([0 max(BWLinVsf(:))+.1])
hold on


%%

%Each dot is a different pool from the T(sf,x-x') population. Each pool
%generates a complex cell RF with a corresponding spatial and sf tuning curve.  
%The y-axis is the measured width of the RF (2sigma). The x axis is the
%prediction from the sf tuning curve.  Each hue is a different window size
%of T, in the x dimension.  Each intensity is a different window size of T
%in the sf dimension.

%sig(x) = 1/(2*pi*sig(w))

figure

cdom = [1 0 0; 0 1 0; 0 0 1];

subplot(1,2,1)
loglog(1./(sfpref)/2,2*XsigVsf,'o-','LineWidth',2,'MarkerSize',2,'Color',cdom(1,:))
set(gca,'XTick',[1/16 1/8 .25 .5 1 2],'YTick',[1/16 1/8 .25 .5 1 2])
ylabel('actual RF width')
xlabel('predicted RF width (2sig): .5/(sfpref)')
xlim([.05 1.5]), ylim([.05 1.5])
hold on

subplot(1,2,2)
loglog(2./(2*pi*BWLinVsf),2*XsigVsf,'o-','LineWidth',2,'MarkerSize',2,'Color',cdom(1,:))
set(gca,'XTick',[1/16 1/8 .25 .5 1 2],'YTick',[1/16 1/8 .25 .5 1 2])
ylabel('actual RF width')
xlabel('predicted RF width (2sig): 2/(2*pi*BWsig)')
xlim([.05 1.5]), ylim([.05 1.5])
hold on



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
