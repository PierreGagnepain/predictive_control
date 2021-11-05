function [rstate,rallitem,rcomb,ACC,BIC,AIC,store_item] = tapas_binary_combined_onestep(r,perceptualmodel,obsmodel,np,usecalibration)

% This is the main function to estimate model trajectories and parameters

% INPUT
% ====================================================
% "r" is the input structure as defined in tapas HGF toolbox and also
% contain item identity information (see start_computational_dcm.m for an example on how to configure "r")

% "perceptualmodel": tapas perceptual model
% "obsmodel": tapas observation model

% There is an option to use pre-calibrated prior different from tapas
% config that should be defined in the variable "np" if "usecalibration" is set to 1

% OUTPUT
% ====================================================
% rstate = "state" model 
% ritem = "item" model 
% rcomb = "combined" model (state and item models)

% each structure contains belief trajectories, fitted perceptual and
% observation parameters, and LME

varargin = 'trans';
u = r.u;
r.u = u(:,1); % to prevent error with softmax obs function

% Default perceptual model
% ~~~~~~~~~~~~~~~~~~~~~~~~
r.c_prc = eval(perceptualmodel);

% Default observation model
% ~~~~~~~~~~~~~~~~~~~~~~~~~
r.c_obs = eval(obsmodel);

% Default optimization algorithm
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r.c_opt = tapas_quasinewton_optim_config;


% work out parameters
opt_idx = [r.c_prc.priorsas, r.c_obs.priorsas];
opt_idx(isnan(opt_idx)) = 0;
opt_idx = find(opt_idx);
d       = length(opt_idx);


traj = [];

% define items
item    = unique(u(:,2));
nitem   = length(item);

% STATE belief
% ------------

% Estimate mode of posterior parameter distribution (MAP estimate)
if usecalibration
    r.c_prc.priormus = np.state.c_prc.priormus;
    r.c_prc.priorsas = np.state.c_prc.priorsas;
end
rstate      = optim(r,r.c_prc.prc_fun , r.c_obs.obs_fun, r.c_opt.opt_algo); % tapas optimization function (see tapas_fitModel)


% get trajectory from optimal parameters (i.e. at MAP estimate)
optimstate      = rstate.optim.final;
hessianstate    = rstate.optim.H;
[trajstate, infStates_state] = r.c_prc.prc_fun(r, optimstate, varargin);

% Separate perceptual and observation parameters
n_prcpars = length(r.c_prc.priormus);
ptrans_prc = rstate.optim.final(1:n_prcpars);
ptrans_obs = rstate.optim.final(n_prcpars+1:end);

rstate.p_prc.ptrans = ptrans_prc;
rstate.p_obs.ptrans = ptrans_obs;

muhat_state         = infStates_state(:,1,1); % first level belief

if size(infStates_state,3)>1 % for HGF
    uncertainty_state   = infStates_state(:,1,2); % uncertainty of first level belief
else % for RW (i.e. no uncertainty)
    uncertainty_state   = ones(size(infStates_state,1),1);
end

% store
rstate.belief               = muhat_state;
rstate.uncertainty_state    = uncertainty_state;


% Compute ACCURACY for state
% --------------------------
[ACC_state,BIC_state,AIC_state]     = get_beliefacc(muhat_state,r.u,opt_idx);


% Item belief
% ------------
muhat_item = [];uncertainty_item = [];
store_item = {}; lme_item = [];par_item = [];prc_item = [];obs_item = [];
for i = 1:nitem
    
    ritem   = [];
    iditem  = find(u(:,2) == item(i));
    y       = r.y(iditem,1);
    ritem   = r;
    ritem.y = y;
    ritem.u = u(iditem,1);
    
    % Estimate mode of posterior parameter distribution (MAP estimate)
    if usecalibration
        ritem.c_prc.priormus = np.item.c_prc.priormus;
        ritem.c_prc.priorsas = np.item.c_prc.priorsas;
    end
    ritem      = optim(ritem, r.c_prc.prc_fun, r.c_obs.obs_fun, r.c_opt.opt_algo);
    optimitem  = ritem.optim.final;
    
    % Separate perceptual and observation parameters
    n_prcpars = length(r.c_prc.priormus);
    ptrans_prc = ritem.optim.final(1:n_prcpars);
    ptrans_obs = ritem.optim.final(n_prcpars+1:end);
    
    ritem.p_prc.ptrans = ptrans_prc;
    ritem.p_obs.ptrans = ptrans_obs;
    
    % get trajectory from optimal parameters (i.e. at MAP estimate)
    [trajitem, infStates_item] = r.c_prc.prc_fun (ritem, optimitem, varargin);
    
    muhat_item (iditem,1)         = infStates_item(:,1,1);% first level belief
    
    if size(infStates_item,3)>1 % for HGF
        uncertainty_item(iditem,1)   = infStates_item(:,1,2); % uncertainty of first level belief
    else % for RW (i.e. no uncertainty)
        uncertainty_item(iditem,1)   = ones(size(infStates_item,1),1);
    end
    
    store_item{i}   = ritem; % for bayesian paramater averaging
    lme_item(i)     = ritem.optim.LME;
    par_item(i,:)   = ritem.optim.final(opt_idx);
    prc_item(i,:)   = ritem.optim.final(1:n_prcpars);
    obs_item(i,:)   = ritem.optim.final(n_prcpars+1:end);
end


% Unified item-model (i.e. combined item to get ACC)
% -------------------------------------------------

% incorporate combined belief
rallitem           = [];
rallitem           = r;
rallitem.belief    = muhat_item;

ptrans_prc = mean(prc_item);
ptrans_obs = mean(obs_item);

rallitem.p_prc.ptrans = ptrans_prc;
rallitem.p_obs.ptrans = ptrans_obs;

% Compute ACCURACY for item
% --------------------------
[ACC_item,BIC_item,AIC_item]     = get_beliefacc(muhat_item,r.u,opt_idx);


% Combined belief using precision (EQ 8&9 Weilnhammer et al.2018)
% ---------------------------------------------------------------

% take first 2 levels to combine beliefs
% here we only combine item and state belief after the first repetition
% (i.e. item beliefs can only influence state after we see it once)

% estimate weight of State vs Item ?
W    = [1 1]; % use similar weight (but interesting to consider different weight for state/item ?

% sahat first level == uncertainty (i.e. inverse of precision)
precison_item = 1./uncertainty_item;
precison_state = 1./uncertainty_state;

numerator    = (W(1)*muhat_state.*precison_state) + (W(2)*muhat_item.*precison_item); % EQ 8&9 Weilnhammer et al.2018
denominator  = precison_state + precison_item;
infStates    = [];
infStates    = numerator./denominator;

% replace first presentation with state
infStates(1:nitem) = muhat_state(1:nitem);


% Compute LME from combined belief
% ---------------------------------------------------------------
% incorporate combined belief
rcomb           = [];
rcomb           = r;
rcomb.belief    = infStates;

ptrans_prc      = nanmean([rstate.p_prc.ptrans;rallitem.p_prc.ptrans]); % not used
ptrans_obs      = nanmean([rstate.p_obs.ptrans;rallitem.p_obs.ptrans]);

rcomb.p_prc.ptrans = ptrans_prc;
rcomb.p_obs.ptrans = ptrans_obs;

% Compute ACCURACY for combined
[ACC_comb,BIC_comb,AIC_comb]     = get_beliefacc(rcomb.belief,r.u,opt_idx);


% All ACCURACY
% --------------------------------
ACC = [ACC_state,ACC_item,ACC_comb];
BIC = [BIC_state,BIC_item,BIC_comb];
AIC = [AIC_state,AIC_item,AIC_comb];


