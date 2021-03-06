function [slope base] = RtoG_trx

global ACQinfo bw

nc = pepgetnoconditions;
pepsetrepeat(0);

covmat_12 = 0;
covmat_22 = 0;
revcor1 = 0;
revcor2 = 0;

z = 0;
for c = 0:nc-1

    pepsetcondition(c);
    nr = pepgetnorepeats;
    for r = 0:nr-1
        pepsetrepeat(r)
        CHs = GetTrialData([1 1 0]);
        for k = 1:length(CHs{1}(1,1,:))
            dum1 = CHs{1}(:,:,k);
            dum2 = CHs{2}(:,:,k);
            
            covmat_12 = covmat_12 + dum2(:);
            covmat_22 = covmat_22 + dum2(:).^2;  %recursively compute inverse of covariance matrix 
            
            revcor1 = revcor1 + dum1(:).*dum2(:); %recursively compute cross-correlation
            revcor2 = revcor2 + dum1(:);
            
            z = z + 1;
        end
    end
end

covmat_11 = z;
covmat_12 = -covmat_12;
covmat_21 = covmat_12;

covmatden = covmat_11.*covmat_22 - (-covmat_21).^2;

dim = size(bw);

slope = covmat_11.*revcor1 + covmat_12.*revcor2;
base = covmat_21.*revcor1 + covmat_22.*revcor2;

slope = slope./covmatden;
base = base./covmatden;

slope = reshape(slope,dim(1),dim(2));
base = reshape(base,dim(1),dim(2));


