function kern = getrandoriRF(CHs,bwCell1)

%42allows you to have a 'dtau' that is bigger than the acquisition period, so that
%it averages.  The code is simpler than 3, but slightly slower

global ACQinfo

%%%%

masklabel = bwlabel(bwCell1);
celldom = unique(masklabel);
celldom = celldom(1:end);

%%%%

dtau = 300;
taudom = -0:dtau:1000;  %This needs to have an element at 0

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %ms per acquired frame
ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

[oriseq sfseq phaseseq] = stimseq(0:pepgetnoconditions-1);  %Get stimulus sequence
%oriseq = stimseq(0:4);  %Get stimulus sequence

oridom = unique([oriseq{1} oriseq{2} oriseq{3}]);
sfdom = unique([sfseq{1} sfseq{2} sfseq{3}]);
phasedom = unique([phaseseq{1} phaseseq{2} phaseseq{3}]);

% if round(oridom(end)) == 999
%     oridom = oridom(1:end-1);
% end

%%%%%%%%%%%%%%%%%%%


kernsmooth = zeros(length(oridom)-1,1);
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

    countmat{p} = zeros(length(oridom),length(sfdom),length(phasedom),length(taudom));
    kern{p} = zeros(length(oridom),length(sfdom),length(phasedom),length(taudom));

    for T = 1:NT
        pepsettimetag(T-1);
        cond = pepgetcondition+1;
        cond = cond(1);

        clear tcourse
        for z = 1:length(CHs{T}{1}(1,1,:))
            CHsdum = CHs{T}{1}(:,:,z);
            %tcourse(z) = (mean(CHsdum(idcell)) - mean(CHsdum(:)))/std(CHsdum(:));
            tcourse(z) = mean(CHsdum(idcell));
            
        end

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
        sfseqdum = sfseq{cond}(1:hper:end);
        phaseseqdum = phaseseq{cond}(1:hper:end);

        for ori = 1:length(oridom)
            for sf = 1:length(sfdom)
                for phase = 1:length(phasedom)
                    
                    id = find(oriseqdum == oridom(ori) & sfseqdum == sfdom(sf) & phaseseqdum == phasedom(phase));
                    stimes = (id-1)*Tupdate; %Stimulus times
                    
                    for i = 1:length(stimes)
                        
                        for tauid = 1:length(taudom)
                            
                            idx = find(tdom_pix>stimes(i)+taudom(tauid)-dtau/2 & tdom_pix<stimes(i)+taudom(tauid)+dtau/2);
                            
                            if ~isempty(idx)
                                kern{p}(ori,sf,phase,tauid) = kern{p}(ori,sf,phase,tauid) + mean(tcourse(idx));
                                countmat{p}(ori,sf,phase,tauid) = countmat{p}(ori,sf,phase,tauid) + 1;
                            end
                        end
                    end
                end
            end
        end
    end
    %kern{p} = kern{p}./countmat{p};
    
    kern{p} = kern{p}(1:end-1,1:end-1,1:end-1,:);  %get rid of blanks
    kern{p} = kern{p}(:,1:end-1,:,:);  %Get rid of last spatial frequency
    countmat{p} = countmat{p}(1:end-1,1:end-1,1:end-1,:);  %get rid of blanks
    countmat{p} = countmat{p}(:,1:end-1,:,:);
    
    sfdom2 = sfdom(1:end-2);
    oridom2 = oridom(1:end-1);
    phasedom2 = phasedom(1:end-1);
    
     
    kern{p} = kern{p}(1:2:end,:,:,:) + kern{p}(2:2:end,:,:,:);
    countmat{p} = countmat{p}(1:2:end,:,:,:) + countmat{p}(2:2:end,:,:,:);
    kern{p} = kern{p}(:,1:2:end,:,:) + kern{p}(:,2:2:end,:,:);
    countmat{p} = countmat{p}(:,1:2:end,:,:) + countmat{p}(:,2:2:end,:,:);
    
    sfdom2 = (sfdom2(1:2:end) + sfdom2(2:2:end))/2;
    oridom2 = (oridom2(1:2:end) + oridom2(2:2:end))/2;
    phasedom2 = phasedom2(1:end-1);

    kern{p} = kern{p}./countmat{p};
    
    RF{p} = MakeRF(squeeze(kern{p}(:,:,:,2)),oridom2,sfdom2,phasedom2);

    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    imagesc(RF{p})
    
    drawnow
    
%     [ma idma] = max(kernplot,[],3);
%     [bestoriid bestsfid] = find(ma == max(ma(:)));
%     besttauid = idma(bestoriid,bestsfid);
%     sforikern = kernplot(:,:,besttauid);    
%     
%     tcori = squeeze(mean(kernplot,2));
%     tcsf = squeeze(mean(kernplot,1));
    %tcori = squeeze(kernplot(:,bestsfid,:));
    %tcsf = squeeze(kernplot(bestoriid,:,:));

    %imagesc(taudom,oridom(1:end-1),tcori)
    %imagesc(taudom,sfdom(1:end-1),tcsf)
    %imagesc(squeeze(kernplot(bestoriid,:,:)))
    %imagesc(sforikern)

end

