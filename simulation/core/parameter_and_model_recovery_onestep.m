function parameter_and_model_recovery_onestep(idx)

% idx  = which of the 200 simulation to recover ?

% This function is intended to run on our computing grid server using
% parallel job

% if you run all simulations at once (i.e. idx = 1:200), it will take
% forever

% Output of these 200 simulations is already save in "result_simulation" directory

% OUTPUT:
% for each virtual participant, we store a "fittedsim" structure with fields:
% - estimated_par (estimated perceptual for each of the 100 random sampling)
% - LME_simulation (log-model evidence) for the 3 simulated and fitted "State",
% "Item", and "Combined" models. LME_simulation is a block diagonal concatenation of the corresponding 3*3
% matrix across the 3 main models (HGF, RW, and KF)

% please see /simulation/pipeline_fig2/figure2.m to plot the outcomes of these simulations 

maxNumCompThreads(4);
rng(idx)

% set path
mainpath        = .../predictive_control/; % define you own !
rootpath        = fullfile(mainpath,'simulation'); 
resultpath      = fullfile(rootpath,'result_simulation');
cd(rootpath)
addpath(genpath(mainpath))

% ==========================================
% SIMULATED DATA
% ==========================================
analysisname = 'simulation_belief2intrusion_beta_vfinale_v1';

% get model simulated input (already defined by model_falsification.m, just load them)
% -----------------------------------------------------------------------------
fn = fullfile(rootpath,'store', sprintf('%s.mat',analysisname));
load(fn)

simudata    = simulation{idx}; % simulation
itempos     = input.itempos;% virtual trials (see model_falsification.m)
nitems      = input.nitems;% virtual number of distinct items (see model_falsification.m)
nsampling   = input.nsampling;
[parname,modname] = get_parname();

% ==========================================
% DEFINE MODELS
% ==========================================

% models
toeval              = {'state','item','comb'};
perceptualmodel     = {};
perceptualmodel{1}  = 'tapas_hgf_binary_2lev_config';
perceptualmodel{2}  = 'tapas_rw_binary_config';
perceptualmodel{3}  = 'tapas_kf_config';
obsmodel            = 'tapas_beta_obs_config';

% name for storage
fnstore = fullfile(resultpath,sprintf('%s_onestep+traj_%d.mat',analysisname,idx));

% run only if file does not exist

if exist(fnstore) == 0
    
    
    % initialize
    fittedsim       = [];
    estimated_par   = {};
    ACC_simulation  = [];
    BIC_simulation  = [];
    belief_fit      = [];
    for ns = 1:nsampling
        
        
        ACC_all = nan(3,3);BIC_all = nan(3,3); modelinc = [0 3 6];
        for m = 1:length(perceptualmodel)
            
            acc = [];bic = [];
            for e = 1:length(toeval)
                
                % generated input
                U = eval(sprintf('simudata.%s.Usim_%s',modname{m},toeval{e}));
                U = U(:,ns);
                
                % FIT GENERATED DATA WITH ALL OBSERVATION MODELS
                % ----------------------------------------------
                
                                
                % prepare data
                r       = [];
                r.u     = [U,itempos'];
                r.y     = U;
                r.ign   = [];
                r.irr   = [];
                
                % FIT MODEL TRAJECTORIES, PARAMETER and GET LME
                % ----------------------------------------------
                % here use evalc to avoid annoying printed output
                evalc('[rstate,ritem,rcomb,ACC,BIC,AIC,store_item] = tapas_binary_combined_onestep(r,perceptualmodel{m},obsmodel,[],0)');
                acc  = [acc;ACC];
                bic  = [bic;BIC];
                
                % store belief
                eval(sprintf('belief_fit.(modname{m}).%s_state(ns,:)    = rstate.belief;',toeval{e}))
                eval(sprintf('belief_fit.(modname{m}).%s_item(ns,:)     = ritem.belief;',toeval{e}))
                eval(sprintf('belief_fit.(modname{m}).%s_comb(ns,:)     = rcomb.belief;',toeval{e}))
                
                % store parameters
                pidx = parname{m}.idx;
                if m ==1 % HGF 2 levels
                    
                    estimated_par{m}.om(ns,:,e)     = [rstate.p_prc.ptrans(pidx(1)) ritem.p_prc.ptrans(pidx(1)) rcomb.p_prc.ptrans(pidx(1))];
                    
                    for ii = 1:nitems
                        estimated_par{m}.om_item(ns,ii,e)     = store_item{ii}.p_prc.ptrans(pidx(1));
                    end
                elseif m ==2 % RW models
                    
                    estimated_par{m}.al(ns,:,e)     = [rstate.p_prc.ptrans(pidx(1)) ritem.p_prc.ptrans(pidx(1)) rcomb.p_prc.ptrans(pidx(1))];
                    for ii = 1:nitems
                        estimated_par{m}.al_item(ns,ii,e)     = store_item{ii}.p_prc.ptrans(pidx(1));
                    end
                elseif m ==3 % KF
                    
                    estimated_par{m}.expom(ns,:,e)     = [rstate.p_prc.ptrans(pidx(1)) ritem.p_prc.ptrans(pidx(1)) rcomb.p_prc.ptrans(pidx(1))];
                    estimated_par{m}.pi_u(ns,:,e)     = [rstate.p_prc.ptrans(pidx(2)) ritem.p_prc.ptrans(pidx(2)) rcomb.p_prc.ptrans(pidx(2))];
                    
                    for ii = 1:nitems
                        estimated_par{m}.expom_item(ns,ii,e)     = store_item{ii}.p_prc.ptrans(pidx(1));
                        estimated_par{m}.pi_u_item(ns,ii,e)     = store_item{ii}.p_prc.ptrans(pidx(2));
                    end
                    
                end    
                
            end
            
            % store
            ACC_all([1:3]+modelinc(m),[1:3]+modelinc(m))    = acc;
            BIC_all([1:3]+modelinc(m),[1:3]+modelinc(m))    = bic;
            
        end
        ACC_simulation(:,:,ns)  = ACC_all;
        BIC_simulation(:,:,ns)  = BIC_all;
        
        % display
        disp(sprintf('SIMULATION #%d sampling:%d done ...',idx,ns))
    end
    
    % store & save
    fittedsim.ACC_simulation    = ACC_simulation;
    fittedsim.BIC_simulation    = BIC_simulation;
    fittedsim.estimated_par     = estimated_par;
    fittedsim.belief_fit        = belief_fit;
    save(fnstore,'fittedsim')
    
end






