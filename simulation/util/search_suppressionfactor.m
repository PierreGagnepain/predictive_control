function MIN = search_suppressionfactor(input,options,meansup,prc_fun,parname)

% MIN = Error between real and simulated intrusion profile
% The grid search will select the suppression factor (meansup) that
% minimizes the error


% real data (4*TNT Sessions, average intrusion rating across whole
% population, Mary et al., Science, 2020)
realint = [0.3886    0.3109    0.2656    0.2474];

% generate simulated data
meansup = repmat(meansup,1,3);
output  = model_generation(input,options,parname,prc_fun,meansup);

% unpack output
fieldoutput = fieldnames(output);
for f = 1:length(fieldoutput)
    eval(sprintf('%s = output.%s;',fieldoutput{f},fieldoutput{f}))
end

% Compute Error between real and generated data, to select suppression
% factor that minimize the error

PRstate = [];PRitem = [];PRcomb = [];
for b = 1:size(Usim_comb,2)
    
    x       = realint';X       = [x,ones(4,1)];
    y       = sesavg_state(b,:)';
    B       = regress(y,X);
    yhat    = X*B;
    E       = sum((x-yhat).^2);

    PRstate(b)= E;
    
    y       = sesavg_item(b,:)';
    B       = regress(y,X);
    yhat    = X*B;
    E       = sum((x-yhat).^2);
    PRitem(b)= E;
    
    y       = sesavg_comb(b,:)';
    B       = regress(y,X);
    yhat    = X*B;
    E       = sum((x-yhat).^2);
    PRcomb(b)= E;
end

MIN         = [];
MIN_state   = nanmedian(PRstate);
MIN_item    = nanmedian(PRitem);
MIN_comb    = nanmedian(PRcomb);
MIN         = [MIN_state MIN_item MIN_comb];

