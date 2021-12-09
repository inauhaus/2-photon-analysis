function SF = getSFtuning

global maskS cellS

%%
figure

sfdom = getdomain('s_freq');

oridom = getdomain('ori');

Ncell = length(cellS.mukern);
for i = 1:Ncell
    
    SFx = cellS.mukern{i}(:,1);  %receptive field for ori = 0;
    SFy = cellS.mukern{i}(:,2);  %receptive field for ori = 90;
    mux = mean(SFx);
    muy = mean(SFy);
    SFraw = SFx;
    if muy>mux
        SFraw = SFy;
    end
        
    [param ffit varacc sigma] = Gaussfit(log2(sfdom),SFraw',0);
    
    param(1) = param(1)+log2(sfdom(1));
    
    SF.pos(i) = param(1);  %xposition of SF (degrees)
    
    SF.sig(i) = param(2); %xWidth of SF (1sigma degrees)
    
    SF.ffit(i,:) = ffit; %store the fit
    
    SF.varacc(i) = varacc;  %store the quality of the fit (variance accounted for)
    
    figure(98)
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),i) 
    plot(SFraw)
    hold on
    plot(SF.ffit(i,:))
    axis tight
    title(num2str(2^SF.pos(i)))
    axis off

    
end

