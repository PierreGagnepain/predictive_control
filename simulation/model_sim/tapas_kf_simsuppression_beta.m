function [belief_sim, Usim,precision,vargout] = tapas_kf_simsuppression_beta(nitems,nreps,meansup,nsampling,parname,sdintrusion,nu,startprob,suppressiontype)

% Calculates the trajectories of the agent's representations under the RW
% ADAPTATION OF TAPAS CODE

c                   = tapas_kf_config;
p                   = c.priormus;
p                   = tapas_kf_transp([], p);% Transform paramaters back to their native space


% Unpack parameters
g_0   = p(1);      % Initial gain
mu_0  = startprob;      % Initial hidden state mean
expom = exp(p(3)); % Process variance
pi_u  = p(4);      % Observation precision


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
        samp_par    = unifrnd(a,b,[1 nsampling]);
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
precision   = [];

suppressionfactor = meansup;

for i = 1:nsampling
    
    
    
    % Initialize updated quantities
    da = NaN(n,1); % Prediction error
    g  = NaN(n,1); % Kalman gain
    mu = NaN(n,1); % Hidden state mean
    u   = [];
    
    % Priors
    g(1)  = g_0;
    mu(1) = mu_0;
    
    % update parameters
    expom = exp(expom_sampling(i)); % Process variance
    pi_u  = pi_u_sampling(i);      % Observation precision
    
    u(1)  = double(rand>0.5);
    
    for k = 2:1:n
        %%%%%%%%%%%%%%%%%%%%%%
        % Effect of input u(k)
        %%%%%%%%%%%%%%%%%%%%%%
        
              
        % Infer input based on belief
        suppressedbel   = supfun(mu(k-1),suppressionfactor);
        suppressedbel   = max([0.1 suppressedbel]);
        belief          = [mu(1:k-2);suppressedbel];
        u(k)            = respsim(belief,nu,sdintrusion,u,nitems);
        
        % Prediction error
        da(k) = u(k)-mu(k-1);
        
        % Gain update
        g(k) = (g(k-1) +pi_u*expom)/(g(k-1) +pi_u*expom +1);
        
        % Hidden state mean update
        mu(k) = mu(k-1)+g(k)*da(k);
        
    end
    
    % Predicted value
    muhat = mu;
    
    % shrink belief (to avoid pb with observation model)
    muhat(muhat>0.95)=0.95;
    muhat(muhat<0.05)=0.05;
    
    precision(:,i)      = ones(n,1);
    belief_sim(:,i)     = muhat;
    Usim(:,i)           = u;
end

