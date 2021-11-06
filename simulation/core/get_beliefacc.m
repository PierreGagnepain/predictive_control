function [ACC,BIC,AIC] = get_beliefacc(belief,y,opt_idx)

% Negative log-model accuracy of previous response (see tapas optim function)
negLogLl    = -nansum(y.*log(belief) + (1-y).*log(1-belief));
ACC         = -negLogLl;

% Calculate AIC and BIC (see optim.m in tapas)
d   = length(opt_idx);
ndp = sum(~isnan(y(:,1)));
AIC  = 2*negLogLl +2*d;
BIC  = 2*negLogLl +d*log(ndp);



