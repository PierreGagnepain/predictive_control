function [ACC,BIC,AIC] = get_beliefacc(belief,nu,y,opt_idx)

% beta function
betaprobdensity = @(alpha,beta,x) gamma(alpha+beta)./(gamma(alpha).*gamma(beta)).*x.^(alpha-1).*(1-x).^(beta-1);

alpha           = belief*nu;
beta            = nu - alpha;


% Negative log-model accuracy of previous response (see tapas optim
% function)
% y           = 0.95.*(y-0.5)+0.5;% shrink y
% negLogLl    = -nansum(log(betaprobdensity(alpha,beta,y)));
negLogLl    = -nansum(y.*log(belief) + (1-y).*log(1-belief));
ACC         = -negLogLl;

% Calculate AIC and BIC (see optim.m in tapas)
d   = length(opt_idx);
ndp = sum(~isnan(y(:,1)));
AIC  = 2*negLogLl +2*d;
BIC  = 2*negLogLl +d*log(ndp);



