function pval = getParamVal(psymbol)

global Analyzer

for i = 1:length(Analyzer.P.param)
    if strcmp(psymbol,Analyzer.P.param{i}{1})
    	idx = i;
        break;
    end
end

pval = Analyzer.P.param{idx}{3};
