function [F1 R Rfit] = CondF1_cellMask2(varargin)

global ACQinfo cellS Analyzer idExamp


%% Get trial delimiter: First and last acquisition frame during the stimulus
Flim = getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning
Flim(2) = Flim(1) + 1000*getparam('stim_time');

tdom = (0:size(cellS.cellMat{1},2)-1)*ACQinfo.msPerLine*ACQinfo.linesPerFrame;
[dum Flim(1)] = min(abs(tdom-Flim(1)));
[dum Flim(2)] = min(abs(tdom-Flim(2)));

%% Further truncate the end of the trial to the last full cycle
try
    STIMfrate = Analyzer.framerate; %ms/period
catch
    STIMfrate = 60;
    'Frame rate does not exist in Analyzer. Setting to 60Hz'
end

T = getParamVal('t_period')/STIMfrate*1000; %ms/cycle

maxCycle = max(unique(floor(tdom(1:Flim(2))/T)));
maxID = find(tdom/T>=maxCycle);
maxID = maxID(1);
Flim(2) = maxID;

%% Create the truncated time domain and F1 phase domain
tdom = tdom - tdom(Flim(1));  %Truncate time domain to stimulus window
tdom = tdom(Flim(1):Flim(2));

phi = 2*pi*tdom/T; %F1 phase domain

%% Fourier domain

fs = 1000/(ACQinfo.msPerLine*ACQinfo.linesPerFrame); %samples/sec
fdom = linspace(0,fs,length(tdom)+1); 
fdom = fdom(1:end-1);


%%
W = round(T/(tdom(2)-tdom(1))*1.5); %
shiftN = round(W/4);

expo = exp(1i*phi);

%ce = ce.*hann(length(ce))'; %Smooth in temporal frequency a bit
for c = 1:length(cellS.cellMat)
    
    dumMat = cellS.cellMat{c}(:,Flim(1):Flim(2),:);  %Grab chunk in stimulus window
    
    R{c} = mean(dumMat,3); %Mean over repeats
    
    muR{c} = mean(R{c},2);
    
    R{c} = R{c}-muR{c}*ones(1,length(tdom)); %Subtract mean from each time course
    
    for i = 1:size(R{c},1)
        polyorder = 1; LPsig = 3;
        R{c}(i,:) = processTcourse(R{c}(i,:),LPsig,polyorder);
    end
    %figure
    for i = 1:size(R{c},1) %loop each cell

        [Rfit{c}(i,:) F1{c}(i)] = getFourierSeries(R{c}(i,:),1,phi);

        
        Rdum = R{c}(i,:);
        Rdum = Rdum-mean(Rdum);
       
        Rdum = Rdum/norm(Rdum);
        expo = expo/(norm(expo)/sqrt(2));
        snr{c}(i) = abs(Rdum(:)'*expo(:)); %varies between 0 and 1;
        
%         subplot(ceil(sqrt(size(R{c},1))),ceil(sqrt(size(R{c},1))),i)
%         plot(R{c}(i,:)), hold on, plot(Rfit{c}(i,:))
%         %[r p] = corrcoef(R{c}(i,:),Rfit{c}(i,:));
%         %title(['r= ' num2str(r(1,2))])
%         title(['r= ' num2str(snr{c}(i))])
        
    end
    
    %figure
    for i = 1:size(R{c},1) %loop each cell

        %coh{c}(i) = getPhaseCoherence(R{c}(i,:),phi,W);
        coh{c}(i) = getPhaseCoherence(R{c}(i,:),phi,W,shiftN,1);
        
        %subplot(ceil(sqrt(size(R{c},1))),ceil(sqrt(size(R{c},1))),i)

        %plot(R{c}(i,:)),
        %title(['coh=' num2str(coh{c}(i))])
    end
   
     %Could compute using the peak response of the fit:
%     [dum id] = max(Rfit{c}');
%     Rdelay = rem(phi(id),2*pi);
%     F1{c} = exp(1i*Rdelay(:));

%     Rw = fft(R{c}');
%     Rw(1:2,:) = 0;
%     Rw((F1id-1)*8:end,:) = 0;
%     R{c} = real(ifft(Rw)');


    %F1{c} = Rfit{c}*ce(:);  %Don't use transpose b/c it takes complex conjugate
            
    %%How well do the harmonics fit the response:
    for i = 1:size(R{c},1)
        [cc pp] = corrcoef(R{c}(i,:),Rfit{c}(i,:));
        r(i) = cc(1,2);
        p(i) = pp(1,2);
    end
    %     id = find(r<thresh);
    %     %id = find(p>thresh);
    %     F1{c}(id) = NaN;
    
    %% Set threshold based on coherence

    
    if ~isempty(varargin)
        id = find(coh{c}<varargin{1});
        length(id)
        F1{c}(id) = NaN;
    end
    
%     if ~isempty(varargin)
%         id = find(snr{c}<varargin{1});
%         length(id)
%         F1{c}(id) = NaN;
%     end
    
%     if ~isempty(varargin)
%         id = find(r<varargin{1});
%         F1{c}(id) = NaN;
%     end
    
    

%     for i = 1:size(R{c},1)
%         [cc pp] = corrcoef(R{c}(i,:),Rfit{c}(i,:));
%         if coh{c}(i)>.9
%             figure(1)
%             clf
%             plot((R{c}(i,:)))
%             hold on
%             plot(Rfit{c}(i,:))
%             
%             drawnow
%             pause(.5)
%         end
%     end
    
    %%
    
    %%
    
    %Get SNR by comparing F1 energy to distribution local frequency energy
    %     harmE = (abs(Rw(F1id,:)));
    %     W = floor(F1id/2);
    %     localID = [F1id-W:F1id-1 F1id+1:F1id+W];
%     LocalE = abs(Rw(localID,:)); %Get local energy around the harmonic
%     LocalEmu = mean(LocalE);
%     LocalEsig = std(LocalE);
%     SNR = (harmE - LocalEmu)./LocalEsig;    
%     id = find(SNR<3)
%     F1{c}(id) = NaN;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end

%%
figure
for i = 1:length(idExamp)
    subplot(length(idExamp),1,i)
       
    ft = abs(fft(R{1}(idExamp(i),:)'));
    plot(fdom(1:20),ft(1:20),'-k','LineWidth',1.5)
    hold on
    
    ft = abs(fft(R{2}(idExamp(i),:)'));    
    plot(fdom(1:20),ft(1:20),'Color',[.5 .5 .5],'LineWidth',1.5)
    hold on
    
    ft = abs(fft(R{3}(idExamp(i),:)'));    
    plot(fdom(1:20),ft(1:20),'--','Color',[0 0 0],'LineWidth',1.5)
    hold on
    
    ft = abs(fft(R{4}(idExamp(i),:)'));    
    plot(fdom(1:20),ft(1:20),'--','Color',[.5 .5 .5],'LineWidth',1.5)
    hold on
        
    set(gca,'XTick',[1000/T 2000/T 3000/T])
    set(gca,'XTickLabel',{'F1','F2','F3'})
    axis tight
end
legend('Right','Left','Up','Down')


function y = processTcourse(y,LpSig,polyorder)

%Called from Ggetrandposkernel2 and Ggetrevcorrkernel2

id = find(isnan(y));
y(id) = nanmean(y);

%mu = mean(y);
%First subtract a low-order polynomial fit:
id = find(y>median(y)+2*std(y));
ytrunc = y;
%ytrunc(id) = median(y); %Remove outliers (large transient responses) before fitting
yfit = polyfitter(ytrunc,polyorder);
y = y-yfit';

%Linear smoothing
if ~isempty(LpSig)
    hh = abs(fft(fspecial('gaussian',size(y),LpSig)));
    y = ifft(fft(y).*hh);
end

%Nonlinear highpass filter to subtract estimated baseline

% h = fspecial('gaussian', [1 length(y)], 10); %Use stats from heavily smoothed version
% ydum = ifft(fft(y).*abs(fft(h)));
% y = y - ordfilt2(ydum, 20, ones(1,600));

% id = find(y<prctile(y,50) & y>prctile(y,0));
% y = y-median(y); 
% y = y/std(y(id));

%y = zscore(y);
%y = y.*sign(y);

function xfit = polyfitter(x,order)

dom = (0:length(x)-1)';

H = ones(length(dom),order+1);  %last column for DC
for i = 1:order
    H(:,i) = dom.^i;
end

p = inv(H'*H)*H'*x';
xfit = H*p;
