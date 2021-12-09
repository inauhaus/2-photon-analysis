function [sig_x sig_ori] = SFxyPoolingModel(fo,sigEInh,sigEPool,AR,EIratio)

global alpha

%fo = linspace(.5,20,50);

sigx = alpha*fo/(2*pi); %scale invariant model of SF bandwidth

% sigE = sqrt(sigx.^2 + sigEPool.^2)/AR;
% oriE = atan2(sigE,fo)*180/pi;
%figure,loglog(fo,oriE)
%hold on



x = linspace(-180,180,10000);
%sigE = 40;
sigIdomdum = (4:2:20);
sigIdomdum = linspace(.01,1,2);
sigIdomdum = sigEInh;
%sigIdomdum = 1
clear sig_x sig_y sig_ori

for i = 1:length(sigIdomdum) %loop each inhibitory width
    for j = 1:length(fo) %loop SF domain
        
        sigI_x = sqrt((sigx(j)).^2 + (sigEPool+sigIdomdum(i)).^2); %convolution
        sigE_x = sqrt((sigx(j)).^2 + (sigEPool).^2);
        
        Gex_x = exp(-(x).^2/(2*sigE_x^2));
        Gex_x = Gex_x/sum(Gex_x);
        Gin_x = exp(-(x).^2/(2*sigI_x^2));
        Gin_x = Gin_x/sum(Gin_x);
        
        V = Gex_x - Gin_x/EIratio; %voltage response as a function of y dim
        id = find(V<0);
        Srate_x = V;
        Srate_x(id) = 0;
        Srate_x = Srate_x/sum(Srate_x);
        
        sig_x(i,j) = sqrt(sum(Srate_x.*x.^2));  %use second moment for bandwidth
 
        
        %%%%%%%%%%%%%%%%
        
        %sigI_y = sqrt((sigx(j)/AR).^2 + (sigEPool+sigIdomdum(i)).^2); %convolution
        sigI_y = sqrt((sigx(j)/AR).^2 + (sigEPool+sigIdomdum(i)).^2); %convolution
        sigE_y = sqrt((sigx(j)/AR).^2 + (sigEPool).^2);
        
        Gex_y = exp(-(x).^2/(2*sigE_y^2));
        Gex_y = Gex_y/sum(Gex_y);
        Gin_y = exp(-(x).^2/(2*sigI_y^2));
        Gin_y = Gin_y/sum(Gin_y);
        
        V = Gex_y - Gin_y/EIratio; %voltage response as a function of y dim
        id = find(V<0);
        Srate_y = V;
        Srate_y(id) = 0;
        Srate_y = Srate_y/sum(Srate_y);
        
        sig_y(i,j) = sqrt(sum(Srate_y.*x.^2));  %use second moment for bandwidth
        
        %[dum idsig] = min(abs(Srate/max(Srate)-0.6065));
        %sig(i,j) = abs(x(idsig));
        sig_ori(i,j) = atan2(sig_y(i,j),fo(j))*180/pi;  
      
        
    end

    %subplot(5,5,i)
    %loglog(fo,oriE)
end
% figure
% subplot(1,2,1)
% loglog(fo,sig_ori')
% set(gca,'XTick',[.25 .5 1 2 4 8])
% set(gca,'YTick',[8 12 19 29 45])
% ylim([8 50]); xlim([.2 8.2])
% 
% subplot(1,2,2)
% loglog(2*fo/pi,sig_x')
% set(gca,'XTick',[.25 .5 1 2 4 8])
% set(gca,'YTick',[.25 .5 1 2 4 8])
% ylim([.25 8]); xlim([.25 8])
% hold on, plot([.25 8],[.25 8])
%%
% sig = .1;
% dx = linspace(0,.5,10);
% x = linspace(-10,10,1000);
% 
% figure
% for i = 1:length(dx) %loop each inhibitory width
%     sig = dx(i)/2;
%     Gex = exp(-(x+dx(i)/2).^2/(2*sig^2));
%     Gin = exp(-(x-dx(i)/2).^2/(2*sig^2));
%     G = Gex-Gin;
%     G = G/norm(G);
% 
%     
%     gw = abs(fft2(G));
%     
%     semilogx(gw/max(gw))
%     %plot(G)
%     hold on
% end
