function getCountMat(trialdom,cellMat)

%A modification of Ggetrandposkernel2 to just get the number of
%presentations


global ACQinfo Analyzer cellS G_RChandles maskS

%%%%

tauN = str2num(get(G_RChandles.kernelLength,'string'));

%%%Get rid of the glia:

nID = getNeuronMask;  %get the index values for the neurons
masklabel = bwlabel(maskS.neuronmask,4);
celldom = unique(masklabel);

%%%%

ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)


%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);

hper = getparam('h_per');
%hper = 1; %in my stimulus code, the sequences have only one value for each
%presentation, so you can't set hper to one

%expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt]
%load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt])
%load(['F:\neurostuff\log_files\' expt])
rseeds = eval(Analyzer.L.param{1}{2});


logfileroot = get(G_RChandles.logfilePath,'string');
expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];

load([logfileroot Analyzer.M.anim '\' expt])


Tf = 1000/frate;  %Frame period in ms (frate obtained from log file)
Tupdate = Tf*hper;


for s = 1:length(rseeds)

    if exist(['rseed' num2str(s)])
        eval(['bwS = rseed' num2str(s) '.bwseq;']);
        eval(['xS = rseed' num2str(s) '.xseq;']);
        eval(['yS = rseed' num2str(s) '.yseq;']);
        eval(['oriS = rseed' num2str(s) '.oriseq;']);
        eval(['colorS = rseed' num2str(s) '.colorseq;']);

        bwseq{s} = domains.bwdom(bwS);
        xseq{s} =  domains.xdom(xS);
        yseq{s} =  domains.ydom(yS);
        oriseq{s} =  domains.oridom(oriS);
        colorseq{s} =  domains.colordom(colorS);
    end

end


%%%%%%%%%%%%%%%%%%%
oridom = domains.oridom;
xdom = domains.xdom;
ydom = domains.ydom;
bwdom = domains.bwdom;
colordom = domains.colordom;

Ncell = length(nID);

for p = 1:1

    pID = nID(p);

    [idcelly idcellx] = find(masklabel == celldom(p));

    CoM = [mean(idcelly) mean(idcellx)];  %center of mass

    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    countmat{p} = zeros(length(oridom),length(xdom),length(ydom),length(bwdom),length(colordom),length(taudom));

    for trialid = 1:length(trialdom)
        
            T = trialdom(trialid);
            [cond rep] = getcondrep(T);
            
            if ~isempty(cellMat{cond})

                tcourse = squeeze(cellMat{cond}(pID,:,rep)); %only used to get the length here.

                tdom = (0:length(tcourse)-1)*acqPeriod;
                %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain of the pixel relative to onset of first stimulus
                tdom_pix = tdom + tau_xy;

                %idbase = find(tdom_pix<synctimes{cond,rep}(1)*1000 | tdom_pix>(synctimes{cond,rep}(end)*1000+500));

                seedno = Analyzer.loops.conds{cond}.val{1};

                oriseqdum = oriseq{seedno}(:);
                xseqdum = xseq{seedno}(:);
                yseqdum = yseq{seedno}(:);
                bwseqdum = bwseq{seedno}(:);
                colorseqdum = colorseq{seedno}(:);


                for ori = 1:length(oridom)
                    for x = 1:length(xdom)
                        for y = 1:length(ydom)
                            for bw = 1:length(bwdom)
                                for color = 1:length(colordom)

                                    id = find(oriseqdum == oridom(ori) & xseqdum == xdom(x) & yseqdum == ydom(y) & bwseqdum == bwdom(bw) & colorseqdum == colordom(color));

                                    %stimes = (id-1)*Tupdate; %Stimulus times
                                    %stimes = cellS.synctimes{cond,rep}(id)*1000;
                                    stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times (ms)

                                    for i = 1:length(stimes)

                                        [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1))));
                                        %idx1 = find(tdom_pix>=stimes(i)+taudom(1)-dtau/2 & tdom_pix<stimes(i)+taudom(1)+dtau/2);
                                        idx1 = idx1(1);
                                        tpiece = idx1:idx1+length(taudom)-1;

                                        if tpiece(1)>0 & tpiece(end)<length(tcourse)
                                            countmat{p}(ori,x,y,bw,color,:) = countmat{p}(ori,x,y,bw,color,:) + 1;
                                        end

                                    end
                                end
                            end
                        end
                    end
                end
            end

    end

end

kerncount{1} = reshape(countmat{1},[length(oridom) length(xdom) length(bwdom) length(colordom) length(taudom)]);  %get rid of 'y' dimension
for p = 1:Ncell
    kerncount{p} = kerncount{1};
end

cellS.kernCount = kerncount;
