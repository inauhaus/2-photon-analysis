function [f0m1 f0m2 signals] = Gf0meanimage(Tlim,b,varargin);
%Tlim is a 2 element vector coresponding to the times within the repeats
%that the images are averaged.
%b is the baseline subtraction vector.  

global bcond f0
%% Compute mean image across all conditions and repeats

%Find blank condition:
nc = pepgetnoconditions;
bcond = [];
for(i=0:nc-1)
    pepsetcondition(i)
    if(pepblank)       %Identify blank
       bcond = i; 
       break
    end
end

f0 = cell(1,nc);
sig1 = cell(1,nc);

for(c=0:nc-1)
    if c == bcond
        f0{c+1} = [];
        sig1{c+1} = [];
    else
        [f0{c+1} sig1{c+1}] = Gf0image(c,Tlim,b,varargin);
    end
end


%% Now average all the repeats for the signals
if length(varargin)==2
    for(c=1:nc)
        if c == bcond+1
            signals{c} = [];
        else
            nr = length(f0{c});
            sig2 = addtrunc(sig1{c},nr); %sig2 is a matrix where each row is a pixel
            sig2 = sig2./nr;
            signals{c} = sig2;
        end
    end
else
    signals = 0;
end

nr = pepgetnorepeats;
repend = floor(nr/2)
repstart = floor(nr/2)+1

%process first half of repeats (images)
for(c=1:nc)
    if c == bcond+1
        f0m1{c} = [];
    else
        img = f0{c}{1};
        for(r=2:repend)
            img = img+f0{c}{r};
        end
        img = img/repend;
        f0m1{c} = img;
    end
end

%process second half of repeats (images)
for(c=1:nc)
    if c == bcond+1
        f0m2{c} = [];
    else
        img = f0{c}{repstart};
        for r = (repstart+1):nr
            img = img+f0{c}{r};
        end
        img = img/(nr-repend);
        f0m2{c} = img;
    end
end


function y = addtrunc(x,nr)
%Truncates all signals in x to the length of shortest one, and then adds them.
for i = 1:nr
    N(i) = length(x{i}(1,:));
end
shortest = min(N);  %Length of shortest repeat

y = zeros(length(x{1}(:,1)),shortest); %(No. of pixels) x (length of shortest repeat)
for i = 1:nr
    y = y + x{i}(:,1:shortest);
end
   
    
