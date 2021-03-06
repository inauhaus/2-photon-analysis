function kern = photRevCorr4(CHs,bwCell1)

%3 uses the cell mask to create the time courses

%4 allows you to have a 'dtau' that is bigger than the acquisition period, so that
%it averages.  The code is simpler than 3, but slightly slower

global ACQinfo

%%%%

masklabel = bwlabel(bwCell1);
celldom = unique(masklabel);
celldom = celldom(1:end);

%%%%

dtau = 300;
taudom = -dtau:dtau:1000;  %This needs to have an element at 0

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %ms per acquired frame
ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

[oriseq sfseq] = stimseq(0:pepgetnoconditions-1);  %Get stimulus sequence
%oriseq = stimseq(0:4);  %Get stimulus sequence

%oridom = unique([sfseq{1} sfseq{2} sfseq{3}])
oridom = unique([oriseq{1} oriseq{2} oriseq{3}]);

% if round(oridom(end)) == 999
%     oridom = oridom(1:end-1);
% end


kernsmooth = zeros(length(oridom)-1,1);
kernsmooth(1:3) = [.5 1 .5];
kernsmooth = kernsmooth*ones(1,length(taudom));

Ncell = length(celldom);

figure
for p = 1:Ncell

    [idcelly idcellx] = find(masklabel == celldom(p));
    idcell = find(masklabel(:) == celldom(p));

    CoM = [mean(idcelly) mean(idcellx)];  %center of mass

    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    countmat{p} = zeros(length(oridom),length(taudom));
    kern{p} = zeros(length(oridom),length(taudom));

    for cond = 1:pepgetnoconditions
        pepsetcondition(cond-1)

        for rep = 1:length(pepgetnorepeats)
            pepsetrepeat(rep-1)
            
            clear tcourse
            for z = 1:length(CHs{cond,rep}{1}(1,1,:))
                CHsdum = CHs{cond,rep}{1}(:,:,z);
                %tcourse(z) = (mean(CHsdum(idcell)) - mean(CHsdum(:)))/std(CHsdum(:));
                tcourse(z) = mean(CHsdum(idcell));
            end

            %         fdom = linspace(0,1000/acqPeriod,length(tcourse)+1);
            %         fdom = fdom(1:end-1);
            %         figure,plot(fdom,abs(fft(tcourse-mean(tcourse))))

            %tcourse = LFPfilt(tcourse(:)',0,1000/acqPeriod,inf,.5);
            %tcourse = BandErase(tcourse(:)',1000/acqPeriod,.66,.43);  %Breathing artifact

            %         hold on,plot(fdom,abs(fft(tcourse)),'r')
            %         asdf

            %     hW = 21;
            %     hh = zeros(1,length(tcourse));
            %     hh(1:hW) = ones(1,hW);
            %     hh = -hh/sum(hh);
            %     hh(ceil(hW/2)) = hh(11) + 1;
            %     tcourse = ifft(fft(tcourse).*abs(fft(hh')));

            tcourse = zscore(tcourse);

            %Tf = 1000/pepParam('refresh');  %Frame period in ms (CRT)
            Tf = 1000/59.55;  %Frame period in ms  (LCD monitor)

            hper = pepgetparam('h_period');
            hper = hper(1);
            %hper = 1;
            Tupdate = Tf*hper;

            tdom = (0:length(tcourse)-1)*acqPeriod;

            tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000 - 50;   %time domain of the pixel relative to onset of first stimulus

            oriseqdum = oriseq{cond}(1:hper:end);

            for ori = 1:length(oridom)

                id = find(oriseqdum == oridom(ori));

                stimes = (id-1)*Tupdate; %Stimulus times

                for i = 1:length(stimes)
                    
                    for tauid = 1:length(taudom)
                    
                        idx = find(tdom_pix>stimes(i)+taudom(tauid)-dtau/2 & tdom_pix<stimes(i)+taudom(tauid)+dtau/2);

                        kern{p}(ori,tauid) = kern{p}(ori,tauid) + mean(tcourse(idx));
                        countmat{p}(ori,tauid) = countmat{p}(ori,tauid) + 1;
                    
                    end
                end
            end
        end
    end
    kern{p} = kern{p}./countmat{p};
    
    %normer = ones(length(kern{p}(:,1)),1)*mean(kern{p});
    %kern{p} = kern{p}-normer;
    
    kernplot = kern{p}(1:end-1,:);
    %kernblank = kern{p}(end,:);
    figure,imagesc(kernplot)
    asdf
    
    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));

    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    imagesc(taudom,oridom(1:end-1),kernplot)
    drawnow
end


%plot tuning curveof each pixel
tc1 = 0; tc2 = 0;
figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

    kernplot = kern{p}(1:end-1,:);
    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));

    [y x] = find(kernplot == max(kernplot(:)));
    
    %x = 7;

    tc = mean(kernplot(:,x+tc1:x+tc2),2);    
    tc = tc-min(tc);
    
    [OMag(p) OAng(p)] = orifind(tc,oridom(1:end-1));
    
    plot(oridom(1:end-1),tc), hold on, plot(oridom(1:end-1),tc,'o'), xlabel('ori')

end


%Now plot color image of tuning
OMag = OMag-min(OMag);
OMag = OMag/max(OMag);

mag = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
ang = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
for p = 2:Ncell

    kernplot = kern{p}(1:end-1,:);
    [y x] = find(kernplot == max(kernplot(:)));
    tc = mean(kernplot(:,x+tc1:x+tc2),2);

    
    idcell = find(masklabel(:) == celldom(p));
    mag(idcell) = OMag(p);
    ang(idcell) = OAng(p);
    
end

figure,
imagesc(ang,'AlphaData',mag), colormap hsv



%time courses
tc1 = 0; tc2 = 0;
figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

    kernplot = kern{p}(1:end-1,:);
    blankresp = kern{p}(end,:);
    
    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));

    [y x] = find(kernplot == max(kernplot(:)));
    orthy = y + round(length(kernplot(:,1))/2);
    if orthy > length(kernplot(:,1))
        orthy = orthy-length(kernplot(:,1));
    end
    
    %x = 7;

    tcourse_ma = mean(kernplot(y+tc1:y+tc2,:),1);   
    tcourse_orth = mean(kernplot(orthy+tc1:orthy+tc2,:),1);
    
%     tcourse_ma = mean(kernplot);
%     tcourse_orth = blankresp;
    
    plot(tcourse_ma,'-o'), hold on, plot(tcourse_orth,'-or')

end





function [OMag OAng] = orifind(G,oridomain)

R = sum(G.*exp(1i*oridomain*pi/90));
OAng = angle(R);                 %-pi to pi
OAng = OAng + pi*(1-sign(OAng+eps));  %0 to 2pi
OAng = OAng*90/pi;               %0 to 180
OMag = abs(R)/(sum(G));