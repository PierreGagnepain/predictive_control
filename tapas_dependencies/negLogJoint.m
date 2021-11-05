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
logLl = nansum(trialLogLls);
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

return;