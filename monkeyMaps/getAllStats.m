function stats = getAllStats(Q)

%% Marginal stats

Q.meas.Dnorm = Q.meas.Dnorm(:);
Q.meas.width = Q.meas.size;
Q.SI_D.width_2Cmplx = Q.SI_D.size_2Cmplx;

pstrAll = {'sfpref','width','BWlin','BWlog','oriBW','F1F0','Dnorm'};

for i = 1:length(pstrAll)-1
    
    pstr = pstrAll{i};
    eval(['stats.' pstr '.mean = nanmean(Q.meas.' pstr ');']);
    eval(['stats.' pstr '.median = nanmedian(Q.meas.' pstr ');']);
    eval(['stats.' pstr '.SD = nanstd(Q.meas.' pstr ');']);
    eval(['stats.' pstr '.meanlog = nanmean(log2(Q.meas.' pstr '));']);
    eval(['stats.' pstr '.SDlog = nanstd(log2(Q.meas.' pstr '));']);
    
end

idD = find(~isnan(Q.meas.Dnorm)); %These use less cells because it fits to black and white
stats.Dnorm.mean = mean(Q.meas.Dnorm(idD));
stats.Dnorm.median = median(Q.meas.Dnorm(idD));
stats.Dnorm.SD = std(Q.meas.Dnorm(idD));
stats.Dnorm.meanlog = mean((log2(Q.meas.Dnorm(idD))));
stats.Dnorm.SDlog = std((log2(Q.meas.Dnorm(idD))));


stats.Dnorm.quartiles = [prctile(Q.meas.Dnorm(idD),0) prctile(Q.meas.Dnorm(idD),25) prctile(Q.meas.Dnorm(idD),50) prctile(Q.meas.Dnorm(idD),75) prctile(Q.meas.Dnorm(idD),100)];

%% make the table
cols = {'sfpref','width','BWlin','BWlog','oriBW','F1F0','Dnorm'};
rows = {'mean';  'median'; 'SD'; 'meanlog'; 'SDlog'};

for i = 1:length(cols)
   
    eval([cols{i} ' = round([stats.' cols{i} '.mean stats.' cols{i} '.median stats.' cols{i} '.SD stats.' cols{i} '.meanlog stats.' cols{i} '.SDlog]*100)/100;']);
    
end

BWlog(4:5) = NaN;

% sfpref = round([stats.sfpref.mean stats.sfpref.std stats.sfpref.geomean stats.sfpref.geostd]*100)/100;
% width = round([stats.size.mean stats.size.std stats.size.geomean stats.size.geostd]*100)/100;
% BWlin = round([stats.BWlin.mean stats.BWlin.std stats.BWlin.geomean stats.BWlin.geostd]*100)/100;
% BWlog = round([stats.BWlog.mean stats.BWlog.std NaN NaN]*100)/100;
% oriBW = round([stats.oriBW.mean stats.oriBW.std stats.oriBW.geomean stats.oriBW.geostd]*100)/100;
% F1F0 = round([stats.F1F0.mean stats.F1F0.std stats.F1F0.geomean stats.F1F0.geostd]*100)/100;

Tmarg = table(sfpref(:),width(:),BWlin(:),BWlog(:),oriBW(:),F1F0(:),Dnorm(:),'RowNames',rows','VariableNames',cols);

figure, uitable('Data',Tmarg{:,:},'ColumnName',Tmarg.Properties.VariableNames,...
    'RowName',Tmarg.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1],'BackgroundColor',[1 1 1; 1 1 1],'FontName','Times');
%%
% Get the table in string form.
figure
TString = evalc('disp(Tmarg)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);


%% Ratio with SF preference

Q.meas.invwidth = Q.meas.width.^-1;

pstrAll = {'invwidth','BWlin','oriBW','BWlog','F1F0','Dnorm'};

%Don't automate for Dnorm
for i = 1:3 %Use log for these
    
    eval(['stats.' pstrAll{i} '.meanLogRatio = nanmean(log2((Q.meas.' pstrAll{i}  ')./Q.meas.sfpref));']);
    eval(['stats.' pstrAll{i}  '.SDLogRatio = nanstd(log2((Q.meas.' pstrAll{i}  ')./Q.meas.sfpref));'])
    eval(['stats.' pstrAll{i}  '.meanRatio = nanmean((Q.meas.' pstrAll{i} ')./Q.meas.sfpref);'])
    eval(['stats.' pstrAll{i}  '.SDRatio = nanstd((Q.meas.' pstrAll{i}  ')./Q.meas.sfpref);'])

end

for i = 4:5  %only sf log for these
    eval(['stats.' pstrAll{i} '.meanLogRatio = mean(Q.meas.' pstrAll{i} '- log2(Q.meas.sfpref));'])
    eval(['stats.' pstrAll{i} '.SDLogRatio = std(Q.meas.' pstrAll{i} '- log2(Q.meas.sfpref));'])
    eval(['stats.' pstrAll{i} '.meanRatio = mean(Q.meas.' pstrAll{i} '- (Q.meas.sfpref));'])
    eval(['stats.' pstrAll{i} '.SDRatio = std(Q.meas.' pstrAll{i} '- log2(Q.meas.sfpref));'])


end

%Do this instead


stats.F1F0.OverSF_quartiles = [prctile(Q.meas.F1F0./Q.meas.sfpref,0) prctile(Q.meas.F1F0./Q.meas.sfpref,25) prctile(Q.meas.F1F0./Q.meas.sfpref,50) prctile(Q.meas.F1F0./Q.meas.sfpref,75) prctile(Q.meas.F1F0./Q.meas.sfpref,100)];

%Need to get the subset
stats.Dnorm.meanLogRatio = mean(Q.meas.Dnorm(idD) - log2(Q.meas.sfpref(idD)));
stats.Dnorm.SDLogRatio = std(Q.meas.Dnorm(idD) - log2(Q.meas.sfpref(idD)));
stats.Dnorm.meanRatio = mean(Q.meas.Dnorm(idD)./Q.meas.sfpref(idD));
stats.Dnorm.SDRatio = std(Q.meas.Dnorm(idD)./Q.meas.sfpref(idD));
stats.Dnorm.OverSF_quartiles = [prctile(Q.meas.Dnorm(idD)./Q.meas.sfpref(idD),0) prctile(Q.meas.Dnorm(idD)./Q.meas.sfpref(idD),25) prctile(Q.meas.Dnorm(idD)./Q.meas.sfpref(idD),50) prctile(Q.meas.Dnorm(idD)./Q.meas.sfpref(idD),75) prctile(Q.meas.Dnorm(idD)./Q.meas.sfpref(idD),100)];

%% correlation with SF preference

pstrAll = {'invwidth','BWlin','oriBW','BWlog','F1F0','Dnorm'};

idW = find(~isnan(Q.meas.sfpref.*Q.meas.width)); %for the width comparison

for i = 1:1
       
    eval(['[r p] = corrcoef(Q.meas.' pstrAll{i} '(idW),Q.meas.sfpref(idW));'])
    eval(['stats.' pstrAll{i} '.r_sf_linlin = r(1,2);'])
    eval(['stats.' pstrAll{i} '.p_sf_linlin = p(1,2);'])
    eval(['[r p] = corrcoef(log2(Q.meas.' pstrAll{i} '(idW)),log2(Q.meas.sfpref(idW)));'])
    eval(['stats.' pstrAll{i} '.r_sf_loglog = r(1,2);'])
    eval(['stats.' pstrAll{i} '.p_sf_loglog = p(1,2);'])

end

%Don't automate forDnorm
for i = 2:3
       
    eval(['[r p] = corrcoef(Q.meas.' pstrAll{i} ',Q.meas.sfpref);'])
    eval(['stats.' pstrAll{i} '.r_sf_linlin = r(1,2);'])
    eval(['stats.' pstrAll{i} '.p_sf_linlin = p(1,2);'])
    eval(['[r p] = corrcoef(log2(Q.meas.' pstrAll{i} '),log2(Q.meas.sfpref));'])
    eval(['stats.' pstrAll{i} '.r_sf_loglog = r(1,2);'])
    eval(['stats.' pstrAll{i} '.p_sf_loglog = p(1,2);'])

end

for i = 4:5 %only sf log for these
       
    eval(['[r p] = corrcoef(Q.meas.' pstrAll{i} ',Q.meas.sfpref);'])
    eval(['stats.' pstrAll{i} '.r_sf_linlin = r(1,2);'])
    eval(['stats.' pstrAll{i} '.p_sf_linlin = p(1,2);'])
    eval(['[r p] = corrcoef(log2(Q.meas.' pstrAll{i} '),(Q.meas.sfpref));'])
    eval(['stats.' pstrAll{i} '.r_sf_loglog = r(1,2);'])
    eval(['stats.' pstrAll{i} '.p_sf_loglog = p(1,2);'])

end


stats.BWlog.r_sf_linlin = NaN;
stats.BWlog.p_sf_linlin = NaN;


[r p] = corrcoef(Q.meas.Dnorm(idD),Q.meas.sfpref(idD));
stats.Dnorm.r_sf_linlin = r(1,2);
stats.Dnorm.p_sf_linlin = p(1,2);
[r p] = corrcoef((Q.meas.Dnorm(idD)),log2(Q.meas.sfpref(idD))); %don't want to take the log of Dnorm
stats.Dnorm.r_sf_loglog = r(1,2);
stats.Dnorm.p_sf_loglog = p(1,2);

%% linear fit with SF preference

Q.meas.invwidth = Q.meas.size.^-1;

pstrAll = {'invwidth','BWlin','oriBW','BWlog','F1F0','Dnorm'};

%Don't automate forDnorm
for i = 1:3
          
    eval(['params = regress(Q.meas.' pstrAll{i} ',[Q.meas.sfpref ones(length(Q.meas.sfpref),1)]);'])
    eval(['stats.' pstrAll{i} '.lineFit_sf_linlin = params;'])
    eval(['params = regress(log2(Q.meas.' pstrAll{i} '),[log2(Q.meas.sfpref) ones(length(Q.meas.sfpref),1)]);'])
    eval(['stats.' pstrAll{i} '.lineFit_sf_loglog = params;'])

end

for i = 4:5 %only sf log for these
          
    eval(['params = regress(Q.meas.' pstrAll{i} ',[Q.meas.sfpref ones(length(Q.meas.sfpref),1)]);'])
    eval(['stats.' pstrAll{i} '.lineFit_sf_linlin = params;'])
    eval(['params = regress((Q.meas.' pstrAll{i} '),[log2(Q.meas.sfpref) ones(length(Q.meas.sfpref),1)]);'])
    eval(['stats.' pstrAll{i} '.lineFit_sf_loglog = params;'])

end


stats.BWlog.lineFit_sf_linlin = NaN;

params = regress(Q.meas.Dnorm(idD),[Q.meas.sfpref(idD) ones(length(Q.meas.sfpref(idD)),1)]); %don't want to take the log of Dnorm
stats.Dnorm.lineFit_sf_linlin = params;
params = regress((Q.meas.Dnorm(idD)),[log2(Q.meas.sfpref(idD)) ones(length(Q.meas.sfpref(idD)),1)]);
stats.Dnorm.lineFit_sf_loglog = params;


%% correlation with pooling model prediction

Q.SI_D.invwidth_2Cmplx = Q.SI_D.width_2Cmplx.^-1;


pstrAll = {'invwidth','BWlin','oriBW','BWlog','F1F0','Dnorm'};

for i = 1:1
       
    eval(['[r p] = corrcoef(Q.meas.' pstrAll{i} '(idW),Q.SI_D.' pstrAll{i} '_2Cmplx(idW));'])
    eval(['stats.' pstrAll{i} '.rpool_sf_linlin = r(1,2);'])
    eval(['stats.' pstrAll{i} '.ppool_sf_linlin = p(1,2);'])
    eval(['[r p] = corrcoef(log2(Q.meas.' pstrAll{i} '(idW)),log2(Q.SI_D.' pstrAll{i} '_2Cmplx(idW)));'])
    eval(['stats.' pstrAll{i} '.rpool_sf_loglog = r(1,2);'])
    eval(['stats.' pstrAll{i} '.ppool_sf_loglog = p(1,2);'])

end

%Don't automate forDnorm
for i = 2:3
       
    eval(['[r p] = corrcoef(Q.meas.' pstrAll{i} ',Q.SI_D.' pstrAll{i} '_2Cmplx);'])
    eval(['stats.' pstrAll{i} '.rpool_sf_linlin = r(1,2);'])
    eval(['stats.' pstrAll{i} '.ppool_sf_linlin = p(1,2);'])
    eval(['[r p] = corrcoef(log2(Q.meas.' pstrAll{i} '),log2(Q.SI_D.' pstrAll{i} '_2Cmplx));'])
    eval(['stats.' pstrAll{i} '.rpool_sf_loglog = r(1,2);'])
    eval(['stats.' pstrAll{i} '.ppool_sf_loglog = p(1,2);'])

end

for i = 4:5 %only sf log for these
       
    eval(['[r p] = corrcoef(Q.meas.' pstrAll{i} ',Q.SI_D.' pstrAll{i} '_2Cmplx);'])
    eval(['stats.' pstrAll{i} '.rpool_sf_linlin = r(1,2);'])
    eval(['stats.' pstrAll{i} '.ppool_sf_linlin = p(1,2);'])
    eval(['[r p] = corrcoef((Q.meas.' pstrAll{i} '),log2(Q.SI_D.' pstrAll{i} '_2Cmplx));'])
    eval(['stats.' pstrAll{i} '.rpool_sf_loglog = r(1,2);'])
    eval(['stats.' pstrAll{i} '.ppool_sf_loglog = p(1,2);'])

end


stats.BWlog.rpool_sf_linlin = NaN;
stats.BWlog.ppool_sf_linlin = NaN;


[r p] = corrcoef(Q.meas.Dnorm(idD),Q.meas.sfpref(idD));
stats.Dnorm.rpool_sf_linlin = r(1,2);
stats.Dnorm.poiik_sf_linlin = p(1,2);
[r p] = corrcoef((Q.meas.Dnorm(idD)),log2(Q.meas.sfpref(idD))); %don't want to take the log of Dnorm
stats.Dnorm.rpool_sf_loglog = r(1,2);
stats.Dnorm.ppool_sf_loglog = p(1,2);

%% Make joint stats table

clear cols rows
%make the table
cols = {'invwidth'; 'BWlin'; 'BWlog'; 'oriBW'; 'F1F0'; 'Dnorm'};
rows = {'<difference>'; 'r' ; 'p'; 'slope'; 'intercept';'r_pool';'p_pool'};

pvalue_id = [3 7];
id = 1:length(rows);
id(pvalue_id) = [];

for i = 1:length(cols)
    eval([cols{i} ' = [stats.' cols{i} '.meanLogRatio; stats.' cols{i} '.r_sf_loglog; stats.' cols{i} '.p_sf_loglog; stats.' cols{i} '.lineFit_sf_loglog(1); stats.' cols{i} '.lineFit_sf_loglog(2); stats.' cols{i} '.rpool_sf_loglog; stats.' cols{i} '.ppool_sf_loglog; ];']);
    eval([cols{i} '(id) = round(' cols{i} '(id)*100)/100;'])

end

% pid = 3;
% id = 1:length(width);
% id(pid) = [];
% 
% width(id) = round(width(id)*100)/100;
% BWlin(id) = round(BWlin(id)*100)/100;
% BWlog(id) = round(BWlog(id)*100)/100;
% oriBW(id) = round(oriBW(id)*100)/100;
% F1F0(id) = round(F1F0(id)*100)/100;

Tjoint = table(eval(cols{1}),eval(cols{2}),eval(cols{3}),eval(cols{4}),eval(cols{5}),eval(cols{6}),'RowNames',rows,'VariableNames',cols);

figure, uitable('Data',Tjoint{:,:},'ColumnName',Tjoint.Properties.VariableNames,...
    'RowName',Tjoint.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

% Get the table in string form.
figure
TString = evalc('disp(Tjoint)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);



