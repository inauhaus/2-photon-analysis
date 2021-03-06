function kern = getrandorikernel_old(CHs,bwCell1)

%2allows you to have a 'dtau' that is bigger than the acquisition period, so that
%it averages.  The code is simpler than 3, but slightly slower

global ACQinfo Analyzer

%%%%

masklabel = bwlabel(bwCell1);
celldom = unique(masklabel);
celldom = celldom(1:end);

%%%%

dtau = 300;
taudom = -0:dtau:1000;  %This needs to have an element at 0

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %ms per acquired frame
ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];
load(['F:\neurostuff\log_files\' expt])
rseeds = eval(Analyzer.L.param{1}{2});

dim = [length(domains.oridom) length(domains.sfdom) length(domains.phasedom) length(domains.colordom)];

for i = 1:length(rseeds)
    eval(['seq = rseed' num2str(rseeds(i)) ';'])
    
    colorseq{i} = ceil(seq/(dim(1)*dim(2)*dim(3)));  %block N
    
    V3D = seq - dim(1)*dim(2)*dim(3)*(colorseq{i}-1);  %vectorized location in 3D block
    
    phaseseq{i} = ceil(V3D/(dim(1)*dim(2)));  %depth in block
    
    V2D = V3D - dim(1)*dim(2)*(phaseseq{i}-1); %vectorized location in matrix
    
    sfseq{i} = ceil(V2D/(dim(1)));
    
    oriseq{i} = V2D - dim(1)*(sfseq{i}-1);
    
end

for i = 1:length(rseeds)
    
    colorseq{i} = domains.colordom(colorseq{i});   
    phaseseq{i} =  domains.phasedom(phaseseq{i});  
    sfseq{i} =  domains.sfdom(sfseq{i});     
    oriseq{i} =  domains.oridom(oriseq{i});  
 
end
    

% if round(oridom(end)) == 999
%     oridom = oridom(1:end-1);
% end

%%%%%%%%%%%%%%%%%%%
oridom = domains.oridom;
sfdom = domains.sfdom;

kernsmooth = zeros(length(oridom),1);
kernsmooth(1:3) = [.5 1 .5];
kernsmooth = kernsmooth*ones(1,length(taudom));

Ncell = length(celldom);
NT = getnotrials;

figure
for p = 1:Ncell

    [idcelly idcellx] = find(masklabel == celldom(p));
    idcell = find(masklabel(:) == celldom(p));

    CoM = [mean(idcelly) mean(idcellx)];  %center of mass

    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    countmat{p} = zeros(length(oridom),length(sfdom),length(taudom));
    kern{p} = zeros(length(oridom),length(sfdom),length(taudom));

    for T = 1:NT

        clear tcourse
        for z = 1:length(CHs{T}{1}(1,1,:))
            CHsdum = CHs{T}{1}(:,:,z);
            %tcourse(z) = (mean(CHsdum(idcell)) - mean(CHsdum(:)))/std(CHsdum(:));
            tcourse(z) = mean(CHsdum(idcell));
            
        end


        %         fdom = linspace(0,1000/acqPeriod,length(tcourse)+1);
        %         fdom = fdom(1:end-1);
        %         figure,plot(fdom,abs(fft(tcourse-mean(tcourse))))

        tcourse = LFPfilt(tcourse(:)',0,1000/acqPeriod,1,.05);
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
        Tf = 1000/60;  %Frame period in ms  (LCD monitor)

        hper = getparam('h_per');
        hper = hper(1);
        %hper = 1;
        Tupdate = Tf*hper;

        tdom = (0:length(tcourse)-1)*acqPeriod;
        tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000 - 50;   %time domain of the pixel relative to onset of first stimulus
        
        [cond] = getcondrep(T);
%         oriseqdum = oriseq{cond}(1:hper:end);
%         sfseqdum = sfseq{cond}(1:hper:end);

        oriseqdum = oriseq{cond}(1:end);
        sfseqdum = sfseq{cond}(1:end);

        for ori = 1:length(oridom)
            for sf = 1:length(sfdom)

                id = find(oriseqdum == oridom(ori) & sfseqdum == sfdom(sf));
                stimes = (id-1)*Tupdate; %Stimulus times

                for i = 1:length(stimes)

                    for tauid = 1:length(taudom)

                        idx = find(tdom_pix>stimes(i)+taudom(tauid)-dtau/2 & tdom_pix<stimes(i)+taudom(tauid)+dtau/2);
                        
                        if ~isempty(idx)
                            kern{p}(ori,sf,tauid) = kern{p}(ori,sf,tauid) + mean(tcourse(idx));
                            countmat{p}(ori,sf,tauid) = countmat{p}(ori,sf,tauid) + 1;
                        end
                    end

                end
            end
        end
    end
    kern{p} = kern{p}./countmat{p};

    %normer = ones(length(kern{p}(:,1)),1)*mean(kern{p});
    %kern{p} = kern{p}-normer;
    
    kernplot = kern{p}(1:end,1:end,:);
    %kernblank = kern{p}(end,:);
    
    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));

    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    [ma idma] = max(kernplot,[],3);
    [bestoriid bestsfid] = find(ma == max(ma(:)));
    besttauid = idma(bestoriid,bestsfid);
    sforikern = kernplot(:,:,besttauid);    
    
    tcori = squeeze(mean(kernplot,2));
    tcsf = squeeze(mean(kernplot,1));
    %tcori = squeeze(kernplot(:,bestsfid,:));
    %tcsf = squeeze(kernplot(bestoriid,:,:));

    %imagesc(taudom,oridom(1:end-1),tcori)
    %imagesc(taudom,sfdom(1:end-1),tcsf)
    %imagesc(squeeze(kernplot(bestoriid,:,:)))
    imagesc(squeeze(mean(kernplot,2)))
    
    drawnow
end



tc1 = 0; tc2 = 0;
figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

    kernplot = kern{p}(1:end-1,1:end-1,:);
    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));
    
    [ma idma] = max(kernplot,[],3);
    [bestoriid bestsfid] = find(ma == max(ma(:)));
    tau = idma(bestoriid,bestsfid); 
    
    tcsf = mean(kernplot(:,:,tau+tc1:tau+tc2),3);
    tcsf = mean(tcsf(bestoriid+tc1:bestoriid+tc2,:),1);
    tcsf = squeeze(tcsf);    
    tcsf = tcsf-min(tcsf);
    
    plot(sfdom(1:end-1),tcsf), hold on, plot(sfdom(1:end-1),tcsf,'o'), xlabel('sfreq')

end

%plot tuning curveof each cell
tc1 = 0; tc2 = 0;
figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

    kernplot = kern{p}(1:end-1,1:end-1,:);
    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));
    
    [ma idma] = max(kernplot,[],3);
    [bestoriid bestsfid] = find(ma == max(ma(:)));
    tau = idma(bestoriid,bestsfid); 
    
    tcori = mean(kernplot(:,:,tau+tc1:tau+tc2),3);
    tcori = mean(tcori,2);
    tcori = squeeze(tcori);    
    tcori = tcori-min(tcori);
    
    [OMag(p) OAng(p)] = orifind(tcori,oridom(1:end-1));
    
    plot(oridom(1:end-1),tcori), hold on, plot(oridom(1:end-1),tcori,'o'), xlabel('ori')

end


%Now plot color image of tuning
OMag = OMag-min(OMag);
OMag = OMag/max(OMag);

mag = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
ang = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
for p = 2:Ncell

    kernplot = kern{p}(1:end-1,1:end-1,:);
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

    kernplot = kern{p}(1:end-1,1:end-1,:);
    blankresp = kern{p}(end,:);
    
    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));
    
    [ma idma] = max(kernplot,[],3);
    [bestoriid bestsfid] = find(ma == max(ma(:)));
    [mi idmi] = min(kernplot,[],3);
    [worstoriid worstsfid] = find(mi == min(mi(:)));

    %x = 7;

    tcourse_ma = squeeze(kernplot(bestoriid,bestsfid,:));   
    tcourse_orth = squeeze(kernplot(worstoriid,worstsfid,:)); 
    
%     tcourse_ma = mean(kernplot);
%     tcourse_orth = blankresp;
    
    plot(tcourse_ma,'-o'), hold on, plot(tcourse_orth,'-or')

end





function [OMag OAng] = orifind(G,oridomain)

R = sum(G'.*exp(1i*oridomain*pi/90));
OAng = angle(R);                 %-pi to pi
OAng = OAng + pi*(1-sign(OAng+eps));  %0 to 2pi
OAng = OAng*90/pi;               %0 to 180
OMag = abs(R)/(sum(G));