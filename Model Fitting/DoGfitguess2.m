function G = DoGfitguess2

%Double check Initial guesses

global RF logdom;

f = RF;

[ma id] = max(f);
mi = prctile(f,.1);

sig1 = logdom(id);
Amp1 = ma;

sig2 = sig1*.8;
Amp2 = Amp1;

base = mi;

mu1 = 0;

%G = [idx length(f)/6 max(RF)-min(RF) min(RF)];
%G = [sig1 Amp1 sig2 Amp2 base];
G = [sig1 Amp1 sig2 Amp2 mu1];
%G = [length(f)/2 1 0 0];