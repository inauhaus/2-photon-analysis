function [Proj p] = CircCorr_Coherence(theta1,theta2,pflag)

%assumes 0 to 2pi

R = sum(exp(1i*theta1).*exp(-1i*theta2))/length(theta1);

%Mag = abs(R);
Proj = real(R);


%This is equivalent
%Proj = (cos(theta1(:)')*cos(theta2(:)) + sin(theta1(:)')*sin(theta2(:))) / length(theta1);

if pflag
    
    Niter = 1000;
    clear MagBoot ProjBoot
    for i = 1:Niter
        [dum shuffid] = sort(rand(length(theta1),1));
        
        theta2dum = theta2(shuffid);
        
        R = sum(exp(1i*theta1).*exp(-1i*theta2dum))/length(theta1);
        
        MagBoot(i) =abs(R);
        ProjBoot(i) = real(R);
        
    end
    
    if Proj>0
        p  = length(find(ProjBoot > Proj))/length(ProjBoot);
    else
        p  = length(find(ProjBoot < Proj))/length(ProjBoot);
    end
    
else
    
    p = NaN;
    
end
    