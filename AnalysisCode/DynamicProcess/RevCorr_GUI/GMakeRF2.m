function RF = GMakeRF2(kern)

global Analyzer ACQinfo G_RChandles DM

%kern is generated from getTCfromRevCorr4

%%%%

%%%%%%%%%%%%%%%%%%%

oridom = DM.oridom;
sfdom = DM.sfdom;
phasedom = DM.phasedom;

Ncell = length(kern);

%plot sf curves

stimsize = getparam('x_size');

Wcm = stimsize/360*(2*pi*Analyzer.M.screenDist);
pixpercm = Analyzer.M.xpixels/Analyzer.M.screenXcm;
Npix = Wcm*pixpercm;

stimsize = 1/(sfdom(1))/2;
%stimsize = getparam('x_size')
%phi = linspace(-stimsize/2,stimsize/2,Npix);
phi = linspace(0,stimsize,Npix);
phi = phi-mean(phi);
%%
Dec = 1;
RFall = 0;
clear RFiradon
figure
for p = 1:Ncell
    
    Frespdum = kern{p};
    
    
    %Frespdum = Frespdum - min(Frespdum(:));
    %Frespdum = squeeze(Frespdum(:,:,1)+Frespdum(:,:,2)-Frespdum(:,:,3)-Frespdum(:,:,4));
    %Frespdum = squeeze(Frespdum(:,:,1)+Frespdum(:,:,2));
    %Frespdum = randn(size(Frespdum));
    
    
    Frespdum = Frespdum-prctile(Frespdum(:),50);
    id = find(Frespdum(:) < 0);
    Frespdum(id) = 0;
    
    
    
    for ori = 1:length(oridom)
        
        RF{p} = 0;
        for sf = 1:length(DM.sfdom)-1
            
            for phase = 1:length(DM.phasedom)
                
                xphi = sfdom(sf)*phi*2*pi/Dec; %radians
                
                
                RF{p} = RF{p} + Frespdum(ori,sf,phase)*sin(xphi-phasedom(phase)*pi/180);
                
            end
            
        end
        RFiradon(:,ori) = RF{p};  %projection for each exis
    end
    
    %tc = mean(Frespdum,3);
    %         tc = squeeze(mean(Frespdum(:,2:4),2));
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    %[dum id] = max(var(RFiradon));
    %plot(RFiradon(:,id))
    %ylim([-2 2])
    
    RF{p} = iradon(RFiradon,oridom);
    %hold on
    imagesc(RF{p})
    axis off
    
    %plot(mean(RFiradon,2),colors{c})
    
    %plot(squeeze(Frespdum(2,2,:)))
    
    RFall = RFall+RF{p};
    
    
    %imagesc(RF,[prctile(RF(:),1) prctile(RF(:),99)]), axis image
        
    %axis off
    drawnow
    
end

%figure,imagesc(RFall)