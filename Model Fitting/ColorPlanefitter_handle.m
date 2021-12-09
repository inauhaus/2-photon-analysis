function err = ColorPlanefitter_handle(param)

global data domx domy;

alpha = param(1);
beta = param(2);

fit = abs(alpha*domx + beta*domy);

err = sum((fit(:)-data(:)).^2);