function [oriCont_mu oriCont_sig sfCont_mu sfCont_sig] = depthContinuity(oriang,sfpref,dpmask,d)

Nsim = 100;

idR = find( dpmask{d(1)}.*dpmask{d(2)} == 1 );

ori1 = oriang{d(1)}(idR);
ori2 = oriang{d(2)}(idR);

dori = abs(oridiff(ori1*pi/180,ori2*pi/180)); %degrees

for i = 1:Nsim
    [dum shuffid] = sort(rand(1,length(idR)));
    oriangshuff1 = ori1(shuffid);
    [dum shuffid] = sort(rand(1,length(idR)));
    oriangshuff2 = ori2(shuffid);

    dorishuff = abs(oridiff(oriangshuff1*pi/180,oriangshuff2*pi/180)); %degrees
    mudorishuff(i) = mean(dorishuff);

end

oriCont_mu = mean(mudorishuff)/mean(dori);
oriCont_sig = std(mudorishuff/mean(dori));

pori = length(find(mudorishuff<mean(dori))) / Nsim;


%%

sf1 = sfpref{d(1)}(idR);
sf2 = sfpref{d(2)}(idR);

dsf = abs(log2(sf1./sf2));  %octaves

for i = 1:Nsim
    [dum shuffid] = sort(rand(1,length(idR)));
    sfprefshuff1 = sf1(shuffid);
    [dum shuffid] = sort(rand(1,length(idR)));
    sfprefshuff2 = sf2(shuffid);
    
    dsfshuff = abs(log2(sfprefshuff1./sfprefshuff2));  %octaves

    mudsfshuff(i) = mean(dsfshuff);

end

sfCont_mu = mean(mudsfshuff)/mean(dsf);
sfCont_sig = std(mudsfshuff/mean(dsf));

psf = length(find(mudsfshuff<mean(dsf))) / Nsim;




function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;