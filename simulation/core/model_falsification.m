clear

% ==========================================
% DEFINE inputs
% ==========================================
input             = [];
input.nsim        = 200; % number of simulation (i.e. "virtual" participants)
input.nreps       = 8; % repetition of TNT suppression cue
input.nsess       = 4; % number of TNT sessions
input.nitems      = 18; % number of item
input.nsampling   = 100; % number of sampling per simulations
input.itempos     = repmat(1:input.nitems,1,input.nreps); % virtual trials

% ==========================================
% DEFINE PERCEPTUAL MODELS
% ==========================================
% models can be found in the "model_sim" directory
evalmodel   = {@tapas_hgf_binary_simsuppression_beta,...
    @tapas_rw_binary_simsuppression_beta,...
    @tapas_kf_simsuppression_beta};

% ==========================================
% DEFINE SAMPLING DISTRIBUTION
% ==========================================
[parname,modname] = get_parname();

% ==========================================
% OPTIONS
% ==========================================
options                 = [];
options.sdintrusion     = 0.1; % noise added to belief
options.nu              = exp(0); % nu parameter of beta function
options.startprob       = 0.5; % starting probability (assuming intrusion/nointrusion are equally probable at first trial)
options.suppressiontype = 'scal'; % scal = * and sub = -

if strcmp(options.suppressiontype,'scal')
    options.supfun = @(a,b) a*b;
elseif  strcmp(options.suppressiontype,'sub')
    options.supfun = @(a,b) a-b;
end

% ==========================================
% GRID SEARCH OF SUPPRESSION FACTOR
% ==========================================
% perform grid-search to fine-tune the suppression factor

xsup    = [];

run_gridsearch = 0;

if run_gridsearch
    supfact     = repmat([.8:0.01:1],10,1); % adjust depending of scaling or subtracting !!
    xsup_res    = [];
    store_MIN = {};store_MAX = {};
    for m = 1:length(evalmodel);
        MIN = [];
        prc_fun = evalmodel{m}; % perceptual model
        for s = 1:size(supfact,1)
            for ss = 1:size(supfact,2)
                meansup     = supfact(s,ss);
                mi = search_suppressionfactor(input,options,meansup,prc_fun,parname{m});
                MIN(s,:,ss) = mi;
            end
            disp(sprintf('griding factor#%d of model#%d',s,m))
        end
        
        % min
        [jk,idm] = min(squeeze(nanmean(MIN,1))');
        xsup_res(m,:) = supfact(1,idm);
        
        % store
        store_MIN{m} = MIN;
        
    end
    xsup = xsup_res;
    
else
    % pre-computed suppression factor
    % xsup = N perceptual models (hgf,rw,kf) * N observation models (state, items, combined)
    xsup = [0.9600    0.9400    0.9700
        0.9000    0.9300    0.9500
        0.8000    0.9400    0.9700];
end

% ==========================================
% Generate intrusion
% ==========================================


simulation      = {};meansuppression = [];
for r = 1:input.nsim
    
    % display
    disp(sprintf('running simulation #%d',r))
    
    
    % MODELS
    % ======================
    for m = 1:length(evalmodel)
        
        % SUPPRESSION FACTOR
        meansup  = xsup(m,:);%+randn(1)*0.01;
        meansuppression(r,:,m) = meansup;
        
        % Perceptual model
        prc_fun             = evalmodel{m};
        
        % Generate response
        output = model_generation(input,options,parname{m},prc_fun,meansup);
             
        % STORE
        % ======================
        
        % not saved (file too large for github)
        eval(sprintf('simulation{r}.%s.beliefstate_sim = output.beliefstate_sim;',modname{m}))
        eval(sprintf('simulation{r}.%s.beliefitem_sim = output.beliefitem_sim;',modname{m}))
        eval(sprintf('simulation{r}.%s.beliefcomb_sim = output.beliefcomb_sim;',modname{m}))
        eval(sprintf('simulation{r}.%s.Usim_state = output.Usim_state;',modname{m}))
        eval(sprintf('simulation{r}.%s.Usim_item = output.Usim_item;',modname{m}))
        eval(sprintf('simulation{r}.%s.Usim_comb = output.Usim_comb;',modname{m}))
        
        eval(sprintf('simulation{r}.%s.sampled_parameter = output.sampled_parameter;',modname{m}))
        eval(sprintf('simulation{r}.%s.sesavg_state = output.sesavg_state;',modname{m}))
        eval(sprintf('simulation{r}.%s.sesavg_item = output.sesavg_item;',modname{m}))
        eval(sprintf('simulation{r}.%s.sesavg_comb = output.sesavg_comb;',modname{m}))
        
        
    end
    
end

% simulation are save as simulation_belief2intrusion_beta_vfinale_v1.mat
% fs = 'simulation_belief2intrusion_beta_vfinale_v1_withlimits.mat';
% save(fs,'simulation','input')