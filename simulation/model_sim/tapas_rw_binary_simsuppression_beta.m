function [belief_sim, Usim,precision,vargout] = tapas_rw_binary_simsuppression_beta(nitems,nreps,meansup,nsampling,parname,sdintrusion,nu,startprob,suppressiontype)

% Calculates the trajectories of the agent's representations under the RW
% ADAPTATION OF TAPAS CODE

c                   = tapas_rw_binary_config;
p                   = c.priormus;
p                   = tapas_rw_binary_transp([], p);% Transform paramaters back to their native space


% Unpack parameters
v_0 = startprob;%p(1);
al  = p(2);


% settings
n               = nitems*nreps;% Number of trials (including prior)


% sampling parameter space
if nitems > 1
    rowidx = 1;
elseif nitems == 1
    rowidx = 2;
end

% gaussian or uniform draw ?
draw = 1;

% with truncated normal distribution ?
trunc = 0;

for par = 1:length(parname.idx)
    
    if draw ==1
        priorsd     = sqrt(parname.prior_var(rowidx,par)); % see drawn from a Gaussian with prior sd (tapas_sampleModel.m)
        priormean   = parname.prior_mean(rowidx,par);
        ok = 0;
        
        if trunc == 1% truncated normal distribution with lower/upper limits to match subject data
            minmax  = [parname.prior_min(rowidx,par) parname.prior_max(rowidx,par)];
        else
            minmax  = [-inf inf];
        end
        
        samp_par = [];
        for i = 1:nsampling
            ok = 0;
            while ok~=1
                pval        = priormean+randn(1)*priorsd;
                if pval>minmax(1) & pval<minmax(2)
                    samp_par = [samp_par;pval];
                    ok = 1;
                    break;
                end
            end
        end
    elseif draw == 2
        a           = parname.prior_min(rowidx,par);
        b           = parname.prior_max(rowidx,par);
        pvasamp_parl        = unifrnd(a,b,[1 nsampling]);
    end
    pname       = sprintf('%s_sampling',parname.name{par});
    eval(sprintf('%s=samp_par;',pname))
    vargout{par} = eval(pname);
end


% suppression function
if strcmp(suppressiontype,'scal')
    supfun = @(a,b) a*b;
elseif  strcmp(suppressiontype,'sub')
    supfun = @(a,b) a-b;
end

belief_sim  = [];
Usim        = [];

suppressionfactor = meansup;

for i = 1:nsampling
    
    
    % update parameters
    al = tapas_sgm(al_sampling(i),1);
    
    
    % Initialize updated quantities
    v  = NaN(n,1);
    da = NaN(n,1);
    
    % input
    u = [];
    
    % Prior
    v(1) = v_0;
    
    % first response choose randomly (because in RW, belief is based on
    % previous trial, not current one)
    u(1)          = double(rand>0.5);
    
    % Pass through representation update loop
    for k = 2:1:n
        
                
        % Infer input based on belief
        suppressedbel   = supfun(v(k-1),suppressionfactor);
        suppressedbel   = max([0.1 suppressedbel]);
        belief          = [v(1:k-2);suppressedbel];
        u(k)            = respsim(belief,nu,sdintrusion,u,nitems);
        
        % Prediction error
        da(k)   = u(k)-v(k-1);
        
        % Value
        v(k)    = v(k-1)+al*da(k);
        
    end
    precision(:,i)      = ones(n,1);
    belief_sim(:,i)     = v(:,1);
    Usim(:,i)           = u;
end

