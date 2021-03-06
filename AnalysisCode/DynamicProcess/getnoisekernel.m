function [kern] = getnoisekernel(cellMat,synctimes,bwCell1,Ntau)

%3 takes the cell time courses as input, instead of all the images

global ACQinfo Analyzer

%%%%

masklabel = bwlabel(bwCell1);
celldom = unique(masklabel);
celldom = celldom(1:end);

%%%%

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %ms per acquired frame

ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

dtau = acqPeriod+50;
taudom = 0:dtau:dtau*Ntau;
%taudom2 = linspace(0,1600,9)

pixpercmX = Analyzer.M.xpixels/Analyzer.M.screenXcm;
pixpercmY = Analyzer.M.ypixels/Analyzer.M.screenYcm;

xcm = 2*pi*Analyzer.M.screenDist*getparam('x_size')/360;  %stimulus width in cm
xN = round(xcm*pixpercmX);  %stimulus width in pixels
ycm = 2*pi*Analyzer.M.screenDist*getparam('y_size')/360;   %stimulus height in cm
yN = round(ycm*pixpercmY);  %stimulus height in pixels

xN = round(xN/getparam('x_zoom'));  %Downsample for the zoom
yN = round(yN/getparam('y_zoom'));

%%%%%%%%%%%%%%%%%%%

Ncell = length(celldom);

NT = getnotrials;
tauN = 12;
figure
for p = 1:Ncell

    [idcelly idcellx] = find(masklabel == celldom(p));

    CoM = [mean(idcelly) mean(idcellx)];  %center of mass

    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);
    %kern{p} = zeros(1,length(cellMat{1}(p,:)));
    
    kern{p} = zeros(yN,xN,tauN);


    for T = 1:NT
        
        cond = getcondrep(T);
        Stim = getNoiseStim(cond);
        
        tcourse = cellMat{T}(p,:);
        %tcourse = mean(cellMat{T});
        
        %mu = mean(tcourse);
        tcourse = LFPfilt(tcourse(:)',0,1000/acqPeriod,2,.05);
        
        %tcourse = tcourse+mu;  %replace the mean
        
        tcourse = zscore(tcourse);
        
        tdom = (0:length(tcourse)-1)*acqPeriod;
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;
        
        stimes = synctimes{T}*1000; %stimulus times
        
        for k = 1:length(stimes)
            [dum id{k}] = min(abs(stimes(k)-tdom_pix));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for st = 1:length(stimes)
            
            tpiece = tcourse(id{st}(1):id{st}(1)+tauN-1)';
            
            for tau = 1:length(tpiece)
                
                if st~=1
                    kern{p}(:,:,tau) = kern{p}(:,:,tau) + tpiece(tau)*((Stim(:,:,st)-Stim(:,:,st-1)));
                else
                    kern{p}(:,:,tau) = kern{p}(:,:,tau) + tpiece(tau)*Stim(:,:,st);
                end
            end
        end
        
    end
    
    dim = size(kern{p});
    
    clear va
    for i = 1:dim(3)
        
        dum = kern{p}(:,:,i);
        va(i) =  std(dum(:));
        
    end

    [dum idx] = max(va);

    kernplot = mean(kern{p}(:,:,3:5),3);

    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    
    imagesc(kernplot)
    %plot(va)
    
    drawnow
end
