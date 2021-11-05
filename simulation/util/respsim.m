function [resp] = respsim(belief,nu,sdintrusion,y,nitems)

% beta function
betaprobdensity = @(alpha,beta,x) gamma(alpha+beta)./(gamma(alpha).*gamma(beta)).*x.^(alpha-1).*(1-x).^(beta-1);

if length(y)<=1
    beta    = 0;
else
    beta    = 1;
end
 
    
if beta == 0
    
    pyes    = belief+randn(1)*sdintrusion;
    pno     = (1-belief)+randn(1)*sdintrusion;
    resp    = double(pyes>pno);
    resp    = resp(end);
    
elseif beta == 1
    bel             = belief;
    bel(end)        = bel(end)+(randn(1)*sdintrusion); % add noise to belief
    alpha           = bel*nu;
    beta            = nu - alpha;
    
    % Beta precision of current cue (for yes and no responses)
    pyes            = betaprobdensity(alpha(end),beta(end),0.95);
    pno             = betaprobdensity(alpha(end),beta(end),0.05);
    
    % Negative log-model accuracy of previous response (see tapas optim
    % function)
    y           = 0.95.*(y-0.5)+0.5;% shrink y
    
    if length(y) ~= length(belief) % HGF
        negLogLl    = -nansum(log(betaprobdensity(alpha(1:end-1),beta(1:end-1),y')));
    else
        negLogLl    = -nansum(log(betaprobdensity(alpha,beta,y')));
    end
    
    negLogLl_pyes   = -log(pyes);
    negLogLl_pno    = -log(pno);
    
    resp = negLogLl+negLogLl_pyes<negLogLl+negLogLl_pno;
end


