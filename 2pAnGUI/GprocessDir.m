function [dirmap] = GprocessDir(f0,hh)

global Analyzer bsflag symbolInfo G_handles

idsym = symbolInfo.ID(1);
idsym2 = symbolInfo.ID(2);

dom2index = get(G_handles.secCollapse,'value') - 1;  %Use zero to indicate the mean

bflag = stimblank(getnoconditions); %if a blank exists in this experiment
if bflag
    f0blank = f0{end};
    f0(end) = [];
end

for i = 1:length(f0)
    dirdom(i) = Analyzer.loops.conds{i}.val{idsym};
    dom2cond(i) = Analyzer.loops.conds{i}.val{idsym2};
end

idslice = 1:length(f0);
if dom2index>0
    dom2 = eval(Analyzer.L.param{idsym2}{2});
    idslice = find(dom2(dom2index) == dom2cond);
end

dim = size(f0{1});
Tens = zeros(dim(1),dim(2),length(idslice),'single'); %preallocate
for k = 1:length(idslice)
    Tens(:,:,k) = f0{idslice(k)};
end

%I used to normalize by the blank if baseline subtraction was not checked.
%This is now done on the Tensor beforehand in CondTensor3
peak_dF = max(Tens,[],3);  %used for d

mi = min(Tens,[],3);
for k = 1:length(Tens(1,1,:))
    if bflag == 1  %if baseline subtraction button is not set
        Tens(:,:,k) = Tens(:,:,k)-f0blank;
    else
        %Tens(:,:,k) = Tens(:,:,k)-mi;
    end
end


su = sum(abs(Tens),3);
for k = 1:length(Tens(1,1,:))
    Tens(:,:,k) = Tens(:,:,k)./su;
end

id = find(Tens(:)<0);
Tens(id) = 0;

% id = find(Tens(:)<0);
% Tens(id) = 0;

dirmap = zeros(size(f0{1}));
for k = 1:length(idslice)
    %dirmap = dirmap + f0{k}*exp(1i*2*dir(k)*pi/180);    %Linear combination
    dirmap = dirmap + Tens(:,:,k)*exp(1i*dirdom(idslice(k))*pi/180);    %Linear combination
end

%if a filter exists, use it...
if ~isempty(hh)
    id = find(isnan(dirmap));
    dirmap(id) = 0;
    dirmap = ifft2(abs(fft2(hh)).*fft2(dirmap));    
end

