function param = NLsatfitter(x0)

%options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
%param = fminsearch('expfitter2_handle',x0,options);
param = fminsearch('NLsatfitter_handle',x0);


%param = fminsearch('NLsat_softfitter_handle',x0);