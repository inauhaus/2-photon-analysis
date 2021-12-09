%% Ori Tuning vs Local neighborhood

%...as  a function of map CV
id = find(OMag~=0 & ~isnan(OMag) & ~isnan(oriCV_Pop));
figure, 
subplot(2,2,1), scatter(oriCV_Pop(id),1-OMag(id),'.k')
xlabel('orimap CV'),ylabel('ori cell CV')
[R p] = corrcoef(oriCV_Pop(id),1-OMag(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(orisig) & orisig<60 & ~isnan(oriCV_Pop));
subplot(2,2,3), scatter(oriCV_Pop(id),orisig(id),'.k')
xlabel('orimap CV'),ylabel('ori cell Sigma')
[R p] = corrcoef(oriCV_Pop(id),orisig(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

%...as  a function of map CV (peaks)
id = find(OMag~=0 & ~isnan(OMag) & ~isnan(oriCV_pks));
subplot(2,2,2), scatter(sqrt(oriCV_pks(id)),1-OMag(id),'.k')
xlabel('orimap CV (peaks)'),ylabel('ori cell CV')
[R p] = corrcoef(sqrt(oriCV_pks(id)),1-OMag(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(orisig) & orisig<60 & ~isnan(oriCV_pks));
subplot(2,2,4), scatter(sqrt(oriCV_pks(id)),orisig(id),'.k')
xlabel('orimap CV (peaks)'),ylabel('ori cell Sigma')
[R p] = corrcoef(sqrt(oriCV_pks(id)),orisig(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])


%% Sfreq Tuning vs Local neighborhood

sfBWhold = sfBW;
id = find(LPness>.7);
sfBW(id) = NaN;

sfBW_Pophold = sfBW_Pop;
id = find(LPness_Pop>.7);
sfBW_Pop(id) = NaN;


figure
%as a function of sfmap "variance"
id = find(~isnan(sfvar_Pop) & ~isnan(sfBW));
subplot(2,4,1), scatter(sqrt(sfvar_Pop(id)),sfBW(id),'.k')
xlabel('sf map stdev'),ylabel('sf cell bandwidth')
[R p] = corrcoef(sqrt(sfvar_Pop(id)),sfBW(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(Qfac) & ~isinf(Qfac) & ~isnan(sfvar_Pop));
subplot(2,4,5), scatter(sqrt(sfvar_Pop(id)),Qfac(id),'.k')
xlabel('sf map stdev'),ylabel('sf cell Q factor')
[R p] = corrcoef(sqrt(sfvar_Pop(id)),Qfac(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

%as a function of sfmap "variance" (peaks)
id = find(~isnan(sfvar_pks) & ~isnan(sfBW));
subplot(2,4,2), scatter(sqrt(sfvar_pks(id)),sfBW(id),'.k')
xlabel('sf map stddev (peaks)'),ylabel('sf cell bandwidth')
[R p] = corrcoef(sqrt(sfvar_pks(id)),sfBW(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(Qfac) & ~isinf(Qfac) & ~isnan(sfvar_pks));
subplot(2,4,6), scatter(sqrt(sfvar_pks(id)),Qfac(id),'.k')
xlabel('sf map stddev (peaks)'),ylabel('sf cell Q factor')
[R p] = corrcoef(sqrt(sfvar_pks(id)),Qfac(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

%as a function of sf map BW
id = find(~isnan(sfBW_Pop) & ~isnan(sfBW));
subplot(2,4,3), scatter(sfBW_Pop(id),sfBW(id),'.k')
xlabel('sf map bandwidth'),ylabel('sf cell bandwidth')
[R p] = corrcoef(sfBW_Pop(id),sfBW(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])
hold on, plot([0 3],[0 3])

id = find(~isnan(Qfac) & ~isinf(Qfac) & ~isnan(sfBW_Pop));
subplot(2,4,7), scatter(sfBW_Pop(id),Qfac(id),'.k')
xlabel('sf map bandwidth'),ylabel('sf cell Q factor')
[R p] = corrcoef(sfBW_Pop(id),Qfac(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

%as a function of sf map Qfac
id = find(~isnan(Qfac_Pop) & ~isinf(Qfac_Pop) & ~isnan(sfBW));
subplot(2,4,4), scatter(Qfac_Pop(id),sfBW(id),'.k')
xlabel('sf map Qfac'),ylabel('sf cell bandwidth')
[R p] = corrcoef(Qfac_Pop(id),sfBW(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(Qfac_Pop) & ~isinf(Qfac_Pop) & ~isnan(Qfac) & ~isinf(Qfac));
subplot(2,4,8), scatter(Qfac_Pop(id),Qfac(id),'.k')
xlabel('sf map Qfac'),ylabel('sf cell Q factor')
[R p] = corrcoef(Qfac_Pop(id),Qfac(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

sfBW = sfBWhold;
sfBW_Pop = sfBW_Pophold;
%% Cell-by-cell

sfBWhold = sfBW;
id = find(LPness>.7);
sfBW(id) = NaN;

sfBW_Pophold = sfBW_Pop;
id = find(LPness_Pop>.7);
sfBW_Pop(id) = NaN;

figure
id = find(~isnan(orisig) & orisig<60 & ~isnan(Qfac) & ~isinf(Qfac));
subplot(2,2,1), scatter(orisig(id),Qfac(id),'.k')
xlabel('ori fitted sigma'),ylabel('sf Q factor')
[R p] = corrcoef(orisig(id),Qfac(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(OMag) & ~isnan(Qfac) & ~isinf(Qfac));
subplot(2,2,2), scatter(1-OMag(id),Qfac(id),'.k')
xlabel('ori CV'),ylabel('sf Q factor')
[R p] = corrcoef(1-OMag(id),Qfac(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(orisig) & orisig<60 & ~isnan(sfBW));
subplot(2,2,3), scatter(orisig(id),sfBW(id),'.k')
xlabel('ori fitted sigma'),ylabel('sf bandwidth')
[R p] = corrcoef(orisig(id),sfBW(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(OMag) & ~isnan(sfBW));
subplot(2,2,4), scatter(1-OMag(id),sfBW(id),'.k')
xlabel('ori CV'),ylabel('sf bandwidth')
[R p] = corrcoef(1-OMag(id),sfBW(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

sfBW = sfBWhold;
sfBW_Pop = sfBW_Pophold;
%% Map comparisons

sfBW_Pophold = sfBW_Pop;
id = find(LPness_Pop>.8);
sfBW_Pop(id) = NaN;

figure
id = find(~isnan(sfvar_Pop) & ~isnan(oriCV_Pop));
subplot(2,4,1), scatter(sfvar_Pop(id),oriCV_Pop(id),'.k')
ylabel('ori map CV'),xlabel('sf map variance')
[R p] = corrcoef(oriCV_Pop(id),sfvar_Pop(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(sfvar_Pop) & ~isnan(oriCV_pks));
subplot(2,4,5), scatter(sfvar_Pop(id),sqrt(oriCV_pks(id)),'.k')
ylabel('ori map CV (peaks)'),xlabel('sf map variance')
[R p] = corrcoef(oriCV_pks(id),sfvar_Pop(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(sfvar_pks) & ~isnan(oriCV_Pop));
subplot(2,4,2), scatter(sqrt(sfvar_pks(id)),oriCV_Pop(id),'.k')
ylabel('ori map CV'),xlabel('sf map variance (peaks)')
[R p] = corrcoef(oriCV_Pop(id),sfvar_pks(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(sfvar_pks) & ~isnan(oriCV_pks));
subplot(2,4,6), scatter(sqrt(sfvar_pks(id)),sqrt(oriCV_pks(id)),'.k')
ylabel('ori map CV (peaks)'),xlabel('sf map variance (peaks)')
[R p] = corrcoef(oriCV_pks(id),sfvar_pks(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(Qfac_Pop) & ~isnan(oriCV_Pop) & ~isinf(Qfac_Pop));
subplot(2,4,3), scatter(Qfac_Pop(id),oriCV_Pop(id),'.k')
ylabel('ori map CV'), xlabel('sf map Qfac')
[R p] = corrcoef(Qfac_Pop(id),oriCV_Pop(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(Qfac_Pop) & ~isnan(oriCV_pks) & ~isinf(Qfac_Pop));
subplot(2,4,7), scatter(Qfac_Pop(id),sqrt(oriCV_pks(id)),'.k')
ylabel('ori map CV (peaks)'), xlabel('sf map Qfac')
[R p] = corrcoef(Qfac_Pop(id),oriCV_pks(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(sfBW_Pop) & ~isnan(oriCV_Pop));
subplot(2,4,4), scatter(sfBW_Pop(id),oriCV_Pop(id),'.k')
ylabel('ori map CV'), xlabel('sf map BW')
[R p] = corrcoef(sfBW_Pop(id),oriCV_Pop(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(~isnan(sfBW_Pop) & ~isnan(oriCV_pks));
subplot(2,4,8), scatter(sfBW_Pop(id),sqrt(oriCV_pks(id)),'.k')
ylabel('ori map CV (peaks)'), xlabel('sf map BW')
[R p] = corrcoef(sfBW_Pop(id),oriCV_pks(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

sfBW_Pop = sfBW_Pophold;

%%
id = find(~isnan(F1F0) & ~isnan(sfpref));
figure, scatter(F1F0(id),sqrt(sfpref(id)),'.k')
xlabel('modulation ratio (F1/F0)'), ylabel('sfreq preference')
[R p] = corrcoef(F1F0(id),sfpref(id));
title(['R = ' num2str(R(1,2)) '   p = ' num2str(p(1,2))])

id = find(F1F0<1 & ~isnan(sfpref));
sfprefcmplx_mu = mean(sfpref(id));
sfprefcmplx_sig = std(sfpref(id))/sqrt(length(id));
id = find(F1F0>1 & ~isnan(sfpref));
sfprefsimple_mu = mean(sfpref(id));
sfprefsimple_sig = std(sfpref(id))/sqrt(length(id));
figure,bar([sfprefcmplx_mu sfprefsimple_mu])
hold on
errorbar([sfprefcmplx_mu sfprefsimple_mu],[sfprefcmplx_sig sfprefsimple_sig],'.k')

