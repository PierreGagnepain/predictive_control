function r = optim(r, prc_fun, obs_fun, opt_algo)

% Determine indices of parameters to optimize (i.e., those that are not fixed or NaN)
opt_idx = [r.c_prc.priorsas, r.c_obs.priorsas];
opt_idx(isnan(opt_idx)) = 0;
opt_idx = find(opt_idx);

% Number of perceptual and observation parameters
n_prcpars = length(r.c_prc.priormus);
n_obspars = length(r.c_obs.priormus);

% Construct the objective function to be MINIMIZED:
% The negative log-joint as a function of a single parameter vector
nlj = @(p) [negLogJoint(r, prc_fun, obs_fun, p(1:n_prcpars), p(n_prcpars+1:n_prcpars+n_obspars))];

% Use means of priors as starting values for optimization for optimized parameters (and as values
% for fixed parameters)
init = [r.c_prc.priormus, r.c_obs.priormus];

% Check whether priors are in a region where the objective function can be evaluated
[dummy1, dummy2, rval, err] = nlj(init);
if rval ~= 0
    rethrow(err);
end

% Do an optimization run
optres = optimrun(nlj, init, opt_idx, opt_algo, r.c_opt,r);

% Record optimization results
r.optim.init  = optres.init;
r.optim.final = optres.final;
r.optim.H     = optres.H;
r.optim.Sigma = optres.Sigma;
r.optim.Corr  = optres.Corr;
r.optim.negLl = optres.negLl;
r.optim.negLj = optres.negLj;
r.optim.LME   = optres.LME;
r.optim.accu  = optres.accu;
r.optim.comp  = optres.comp;

% Do further optimization runs with random initialization
if isfield(r.c_opt, 'nRandInit') && r.c_opt.nRandInit > 0
    for i = 1:r.c_opt.nRandInit
        % Use prior mean as starting value for random draw
        init = [r.c_prc.priormus, r.c_obs.priormus];

        % Get standard deviations of parameter priors
        priorsds = sqrt([r.c_prc.priorsas, r.c_obs.priorsas]);
        optsds = priorsds(opt_idx);

        % Add random values to prior means, drawn from Gaussian with prior sd
        rng('shuffle');
        init(opt_idx) = init(opt_idx) + randn(1,length(optsds)).*optsds;

        % Check whether initialization point is in a region where the objective
        % function can be evaluated
        [dummy1, dummy2, rval, err] = nlj(init);
        if rval ~= 0
            rethrow(err);
        end

        % Do an optimization run
        optres = optimrun(nlj, init, opt_idx, opt_algo, r.c_opt);

        % Record optimization if the LME is better than the previous record
        if optres.LME > r.optim.LME
            r.optim.init  = optres.init;
            r.optim.final = optres.final;
            r.optim.H     = optres.H;
            r.optim.Sigma = optres.Sigma;
            r.optim.Corr  = optres.Corr;
            r.optim.negLl = optres.negLl;
            r.optim.negLj = optres.negLj;
            r.optim.LME   = optres.LME;
            r.optim.accu  = optres.accu;
            r.optim.comp  = optres.comp;
        end
    end
end

% Calculate AIC and BIC
d = length(opt_idx);
if ~isempty(r.y)
    ndp = sum(~isnan(r.y(:,1)));
else
    ndp = sum(~isnan(r.u(:,1)));
end
r.optim.AIC  = 2*r.optim.negLl +2*d;
r.optim.BIC  = 2*r.optim.negLl +d*log(ndp);

end % function optim

% --------------------------------------------------------------------------------------------------
function [negLogJoint, negLogLl, rval, err] = negLogJoint(r, prc_fun, obs_fun, ptrans_prc, ptrans_obs)
% Returns the the negative log-joint density for perceptual and observation parameters

% Calculate perceptual trajectories. The columns of the matrix infStates contain the trajectories of
% the inferred states (according to the perceptual model) that the observation model bases its
% predictions on.
try
    [dummy, infStates] = prc_fun(r, ptrans_prc, 'trans');
catch err
    negLogJoint = realmax;
    negLogLl = realmax;
    % Signal that something has gone wrong
    rval = -1;
    return;
end

% Calculate the log-likelihood of observed responses given the perceptual trajectories,
% under the observation model
trialLogLls = obs_fun(r, infStates, ptrans_obs);
logLl = sum(trialLogLls, 'omitnan');
negLogLl = -logLl;

% Calculate the log-prior of the perceptual parameters.
% Only parameters that are neither NaN nor fixed (non-zero prior variance) are relevant.
prc_idx = r.c_prc.priorsas;
prc_idx(isnan(prc_idx)) = 0;
prc_idx = find(prc_idx);

logPrcPriors = -1/2.*log(8*atan(1).*r.c_prc.priorsas(prc_idx)) - 1/2.*(ptrans_prc(prc_idx) - r.c_prc.priormus(prc_idx)).^2./r.c_prc.priorsas(prc_idx);
logPrcPrior  = sum(logPrcPriors);

% Calculate the log-prior of the observation parameters.
% Only parameters that are neither NaN nor fixed (non-zero prior variance) are relevant.
obs_idx = r.c_obs.priorsas;
obs_idx(isnan(obs_idx)) = 0;
obs_idx = find(obs_idx);

logObsPriors = -1/2.*log(8*atan(1).*r.c_obs.priorsas(obs_idx)) - 1/2.*(ptrans_obs(obs_idx) - r.c_obs.priormus(obs_idx)).^2./r.c_obs.priorsas(obs_idx);
logObsPrior  = sum(logObsPriors);

negLogJoint = -(logLl + logPrcPrior + logObsPrior);

% Signal that all has gone right
err = [];
rval = 0;

end % function negLogJoint

% --------------------------------------------------------------------------------------------------
function optres = optimrun(nlj, init, opt_idx, opt_algo, c_opt,r)
% Does one run of the optimization algorithm and returns results

% The objective function is now the negative log joint restricted
% with respect to the parameters that are not optimized
obj_fun = @(p_opt) restrictfun(nlj, init, opt_idx, p_opt);

% Optimize
if ~isfield(r,'optim')
    disp(' ')
    disp('Optimizing...')
    optres = opt_algo(obj_fun, init(opt_idx)', c_opt);
    
    % Record initialization point
    optres.init = init;
    
    % Replace optimized values in init with arg min values
    final = init;
    final(opt_idx) = optres.argMin';
    optres.final = final;
    
    % Get the negative log-joint and negative log-likelihood
    [negLj, negLl] = nlj(final);
    
else
    optres.init     = r.optim.init;
    optres.final    = r.optim.final;
    negLj           = r.optim.negLj;
    negLl           = r.optim.negLl;
    optres.argMin   = r.optim.argMin;
    optres.valMin   = r.optim.valMin;
    optres.Sigma    = r.optim.Sigma;
    optres.T        = r.optim.T;
end

% Calculate the covariance matrix Sigma and the log-model evidence (as approximated
% by the negative variational free energy under the Laplace assumption).
disp(' ')
disp('Calculating the log-model evidence (LME)...')
d     = length(opt_idx);

% Numerical computation of the Hessian of the negative log-joint at the MAP estimate
options.init_h    = 1;
options.min_steps = 10;
H = tapas_riddershessian(obj_fun, optres.argMin, options);


% Use the Hessian from the optimization, if available,
% if the numerical Hessian is not positive definite
if any(isinf(H(:))) || any(isnan(H(:))) || any(eig(H)<=0)
    if isfield(optres, 'T')
        % Hessian of the negative log-joint at the MAP estimate
        % (avoid asymmetry caused by rounding errors)
        H = inv(optres.T);
        H = (H' + H)./2;
        % Parameter covariance
        Sigma = optres.T;
        
        % Parameter correlation
        Corr = tapas_Cov2Corr(Sigma);
        % Log-model evidence ~ negative variational free energy
        LME = -optres.valMin + 1/2*log(1/det(H)) + d/2*log(2*pi);
    else
        disp('Warning: Cannot calculate Sigma and LME because the Hessian is not positive definite.')
    end
else
    % Calculate parameter covariance
    Sigma = inv(H);
    Sigma = (Sigma' + Sigma)./2;
    
    % Parameter correlation
    Corr = tapas_Cov2Corr(Sigma);
    % Log-model evidence ~ negative variational free energy
    LME = -optres.valMin + 1/2*log(1/det(H)) + d/2*log(2*pi);
end

% Record results
optres.H = H;
optres.Sigma = Sigma;
optres.Corr = Corr;
optres.negLl = negLl;
optres.negLj = negLj;
optres.LME = LME;

% Calculate accuracy and complexity (LME = accu - comp)
optres.accu = -negLl;
optres.comp = optres.accu -LME;

end % function optimrun

% --------------------------------------------------------------------------------------------------
function val = restrictfun(f, arg, free_idx, free_arg)
% This is a helper function for the construction of file handles to
% restricted functions.
%
% It returns the value of a function restricted to subset of the
% arguments of the input function handle. The input handle takes
% *one* vector as its argument.
% 
% INPUT:
%   f            The input function handle
%   arg          The argument vector for the input function containing the
%                fixed values of the restricted arguments (plus dummy values
%                for the free arguments)
%   free_idx     The index numbers of the arguments that are not restricted
%   free_arg     The values of the free arguments

% Replace the dummy arguments in arg
arg(free_idx) = free_arg;

% Evaluate
val = f(arg);

end % function val