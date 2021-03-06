function hh = makeTemporalfilter(varargin)

global G_RChandles ACQinfo

togstateHP = get(G_RChandles.HPflag,'Value');
togstateLP = get(G_RChandles.LPflag,'Value');

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 

if ~isempty(varargin)
    N = getTrialLength(varargin{1});
else
    N = getTrialLength(1); %Just use the first trial
end

if togstateHP == 1
    Hwidth = str2double(get(G_RChandles.Hwidth,'string'));
    Hwidth = round(Hwidth/acqPeriod); %convert to samples
    ind = get(G_RChandles.HPWind,'value');

    switch ind
        case 1 %Gaussian
            dom = (1:N)-round(N/2);
            H = exp(-dom.^2/(2*Hwidth^2));
            H = -H/sum(H);
            H(round(N/2)) = 1+H(round(N/2));
        case 2 %Hann
            H = zeros(1,N);
            Hd = hann(Hwidth);
            Hd = -Hd./sum(Hd(:));
            Hd(round(Hwidth/2)) = 1+Hd(round(Hwidth/2));
            H(1:Hwidth) = Hd;
        case 3 %Disc
            H = zeros(1,N);
            Hd = -ones(1,Hwidth)/Hwidth;
            Hd(round(Hwidth/2)) = 1+Hd(round(Hwidth/2));
            H(1:Hwidth) = Hd;
    end
    if togstateLP == 0
        hh = abs(fft(H));   %Eliminate phase information
    end
end

if togstateLP == 1
    Lwidth = str2double(get(G_RChandles.Lwidth,'string'));
    ind = get(G_RChandles.LPWind,'value');

    switch ind
        case 1            
            Lwidth = Lwidth/acqPeriod; %convert to samples (can be decimal for Gaussian)
            dom = (1:N)-round(N/2);
            L = exp(-dom.^2/(2*Lwidth^2));
            L = L/sum(L);
        case 2            
            Lwidth = round(Lwidth/acqPeriod); %convert to samples
            L = zeros(1,N);
            Ld = hann(Lwidth);
            Ld = Ld./sum(Ld(:));
            L(1:Lwidth) = Ld;            
        case 3            
            Lwidth = round(Lwidth/acqPeriod); %convert to samples
            L = zeros(1,N);
            Ld = ones(1,Lwidth)/Lwidth;
            L(1:Lwidth) = Ld;            
    end
    if togstateHP == 0
        hh = abs(fft(L(:)))';   %Eliminate phase information
    else
        hh = abs(fft(L(:)).*fft(H(:)))';   %Take mag because phase gives a slight shift.
    end
end

if ~or(togstateLP,togstateHP)
    hh = [];
end



