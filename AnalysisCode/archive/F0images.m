function [Data1 Data2] = F0images(Tlim)

%Compute mean images for each condition, within time limit vector 'Tlim' (ms)

global ACQinfo

[dum dum2 sync] = GetTrialData(1,1);  %Just to get some sync info
imsize = length(sync(:,1));

low = min(sync(:))
high = max(sync(:));
mid = (low(1)+high(1))/2;

nc = pepgetnoconditions;
nr = pepgetnorepeats;

Data1 = cell(1,nc);
Data2 = cell(1,nc);


%Get sample period (ms/pixel)
sp = ACQinfo.msPerLine/ACQinfo.pixelsPerLine %(msPerLine is actually sec/line)

for c = 1:nc
    Data1{c} = 0;
    Data2{c} = 0;
    for r = 1:nr
        
        [Ddum1 Ddum2 sync] = GetTrialData(c,r);

        sync = sign(sync-mid);

        %idstart = find(sync(:)<0);  %Syncs are negative pulses?
        %idstart = idstart(1);
        idstart = 1;
        
        idstart = idstart+Tlim(1)/sp;
        idstop = idstart+Tlim(2)/sp;
        
        imstart = ceil(idstart/imsize);
        imstop = ceil(idstop/imsize);

        Ddum1 = mean(Ddum1(:,imstart:imstop),2);
        Ddum2 = mean(Ddum2(:,imstart:imstop),2);

        Data1{c} = Data1{c} + Ddum1;
        Data2{c} = Data2{c} + Ddum2;
        
    end
end

rows = ACQinfo.linesPerFrame;
cols = ACQinfo.pixelsPerLine;

for c = 1:nc
    Data1{c} = reshape(Data1{c},rows,cols);
    Data2{c} = reshape(Data2{c},rows,cols);
end


