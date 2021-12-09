function [percS rfit] = getOpsinInput(kernPop)

%Expecting each element of kernPop to be a ori x green x blue kernel for
%neuron


gdom = [-1 -.5 -.25 0 .25 .5 1];
bdom = [.1 .3 1];

[gdomM bdomM] = getConeContrast(gdom,bdom); %Convert LED space to opsin space

H = [abs(gdomM(:)) abs(bdomM(:))];


for i = 1:length(kernPop)

    oriTC = squeeze(mean(mean(kernPop{i},2),3));
    oriTC = oriTC/sum(oriTC);
    kernPopW = 0;
    for ori = 1:length(oriTC)       
        kernPopW = kernPopW + kernPop{i}(ori,:,:)*oriTC(ori);        
    end
    y = squeeze(mean(kernPop{i},1))';
    y = y(:);
    xhat = inv(H'*H)*H'*y;  %Linear regression
    xhat = phi(xhat);
    yhat = H*xhat;
    yhat = reshape(yhat,size(gdomM));
    y = reshape(y,size(gdomM));
    
    
    figure(9)
    subplot(ceil(length(kernPop)/8),8,i)   

    cdom = [0 0 0; 1 0 0; 0 0 1];
    for c = 1:length(bdom)
        plot(gdom,y(c,:),'Color',cdom(c,:))
        hold on
        plot(gdom,yhat(c,:),'--','Color',cdom(c,:))
    end
    percS(i) = round(xhat(2)/(xhat(1)+xhat(2))*100);
    title(['%S=' num2str(percS(i))])
    
    mi = min([y(:); yhat(:)]);
    ma = max([y(:); yhat(:)]);
    ylim([mi(1) ma(1)+eps])
    
    rfitdum = corrcoef(y(:),yhat(:));
    rfit(i) = rfitdum(1,2);
    
    figure(10)
    subplot(ceil(length(kernPop)/8),8,i)
    
    [Mdomdum Sdomdum] = meshgrid(-1:.01:1,-1:.01:1);
    fitIm = abs(Mdomdum)*xhat(1) + abs(Sdomdum)*xhat(2);  %Make 2D model
    %fitIm = ((Mdomdum.^2)*pM(E_vecloc) + (Sdomdum.^2)*pS(E_vecloc) + (Sdomdum.*Mdomdum)*pMS(E_vecloc));
    %fitIm = sqrt(fitIm);
    imagesc(Mdomdum(1,:),Sdomdum(:,1)',fitIm)
    xlabel('M'), ylabel('S')
    axis square
    axis off
    title(['%S=' num2str(percS(i))])
    
end

    
    