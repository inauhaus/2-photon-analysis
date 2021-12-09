function [coh] = getPhaseCoherence(sig,phi,W,shiftN,varargin)

if isempty(varargin)
    [e,v] = dpss(W,1);
else
    e = varargin{1};
end



%Generate time frequency matrices

%lazy way of computing the number of shifts
N = length(sig);
for i = 1:N

    stdom = 1+(i-1)*shiftN:(i-1)*shiftN+W;

    if stdom(end)>N
        break
    end

    Nshift = i;

end

win = hann(W);
if isempty(varargin)
    win = ones(W,1);
end

ce = exp(1i*phi(:));
%%
z = 1;
for i = 1:Nshift

    stdom = 1+(i-1)*shiftN:(i-1)*shiftN+W;
    
    sigdum = sig(stdom);
    sigdum = sigdum(:)-mean(sigdum);
    %cedum = ce(stdom).*exp(-1i.*phi(shiftN*(i-1)+1));
    cedum = ce(stdom);

    for s = 1:length(e(1,:))
        
        HH(z) = mean(win.*sigdum.*e(:,s).*cedum);

%         figure,plot(win.*sigdum.*e(:,s))
%         hold on
%         plot(imag(cedum)*150)

        z = z+1;
    end

end

coh = abs(sum(HH))/sum(abs(HH));

