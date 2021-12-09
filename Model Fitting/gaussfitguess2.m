function G = gaussfitguess2

%Double check Initial guesses

global dom RF;

f = RF;
[M,idx] = max(f);
mu = dom(idx);
sig = range(dom)/4;
%sig = length(f)/8;

mi = prctile(f,.1);
ma = max(f);

%G = [idx length(f)/6 max(RF)-min(RF) min(RF)];
G = [mu sig ma-mi mi];
%G = [length(f)/2 1 0 0];