
mainpath    = '/Volumes/Transcend/ComputationLeone/projet_code/predictive_control/simulation';
simname     = 'simulation_belief2intrusion_beta_vfinale_v1';

% load sim data
load('/Volumes/dataemi6/REMEMBER/remember_all/giovanni_leone/computational_model/store/simulation_belief2intrusion_beta_vfinale_v1.mat')


% parameter
item_identity       = repmat(1:18,1,8);
[parname,modname]   = get_parname();
perceptualmodel     = {};
perceptualmodel{1}  = 'tapas_hgf_binary_2lev_config';
perceptualmodel{2}  = 'tapas_rw_binary_config';
perceptualmodel{3}  = 'tapas_kf_config';
priorpercp          = {};

priorpercp{1}       = [NaN 0 0 NaN -2.3026 0 NaN 0 0 0 -Inf NaN -3 0]; % HGF
priorpercp{2}       = [0 0]; % RW
priorpercp{3}       = [-2.9957    0.5000   -7.0000    0.5000]; % KF

Ctraj = [];MItraj = [];Cmaxtraj = [];
for s = 1:200
    smfn = fullfile(mainpath,'result_simulation',sprintf('%s_%d.mat',simname,s));
    load(smfn)
    
    for ss = 1:100
        
        belief_sim = [];
        belief_fit = [];
        for m = 1:length(modname)
            
            % SIMULATED
            % =============
            belief_simstate     = simulation{s}.(modname{m}).beliefstate_sim(:,ss)';
            U_simstate          = simulation{s}.(modname{m}).Usim_state(:,ss)';
            belief_simitem      = simulation{s}.(modname{m}).beliefitem_sim(:,ss)';
            U_simitem           = simulation{s}.(modname{m}).Usim_item(:,ss)';
            belief_simcomb      = simulation{s}.(modname{m}).beliefcomb_sim(:,ss)';
            U_simcomb           = simulation{s}.(modname{m}).Usim_comb(:,ss)';
            
            belief_sim          = [belief_sim,belief_simstate',belief_simitem',belief_simcomb'];
            
            % FITTED
            % =============
            r.ign   = [];
            r.irr   = [];
            r.u     = U_simstate';
            r.c_prc = eval(perceptualmodel{m});

            % STATE
            prc_pvec                    = priorpercp{m};
            for np = 1:length(parname{m}.name)
                pval                        = fittedsim.estimated_par{m}.(parname{m}.name{np})(ss,1,1);
                prc_pvec(parname{m}.idx(np)) 	= pval;
            end
            
            [trajstate, infStates_state]    = r.c_prc.prc_fun(r, prc_pvec, 'trans');
            muhat_state         = infStates_state(:,1,1); % first level belief
            
            if size(infStates_state,3)>1 % for HGF
                uncertainty_state   = infStates_state(:,1,2); % uncertainty of first level belief
            else % for RW (i.e. no uncertainty)
                uncertainty_state   = ones(size(infStates_state,1),1);
            end

            
            % ITEM
            muhat_item=[];uncertainty_item=[];
            for i = 1:18
                
                prc_pvec                    = priorpercp{m};
                for np = 1:length(parname{m}.name)
                    pval                        = fittedsim.estimated_par{m}.([parname{m}.name{np},'_item'])(ss,i,2);
                    prc_pvec(parname{m}.idx(np)) 	= pval;
                end
            
                                
                iditem  = find(item_identity==i);
                r.u     = U_simitem(iditem)';
                [trajstate, infStates_item]    = r.c_prc.prc_fun(r, prc_pvec, 'trans');
                
                
                muhat_item (iditem,1)         = infStates_item(:,1,1);% first level belief
                
                if size(infStates_item,3)>1 % for HGF
                    uncertainty_item(iditem,1)   = infStates_item(:,1,2); % uncertainty of first level belief
                else % for RW (i.e. no uncertainty)
                    uncertainty_item(iditem,1)   = ones(size(infStates_item,1),1);
                end
            end
            
            
            % estimate weight of State vs Item ?
            W    = [1 1]; % use similar weight (but interesting to consider different weight for state/item ?
            
            % sahat first level == uncertainty (i.e. inverse of precision)
            precison_item = 1./uncertainty_item;
            precison_state = 1./uncertainty_state;
            
            numerator    = (W(1)*muhat_state.*precison_state) + (W(2)*muhat_item.*precison_item); % EQ 8&9 Weilnhammer et al.2018
            denominator  = precison_state + precison_item;
            muhat_comb    = [];
            muhat_comb    = numerator./denominator;
            
            belief_fit          = [belief_fit,muhat_state,muhat_item,muhat_comb];
        end
        
        % corr similarity
        C       = corr(belief_sim,belief_fit);
        Ctraj   = cat(3,Ctraj,C);
        [jk,idm]=max(C,[],2);
        Cmax = zeros(9,9);
        for c = 1:size(C,1)
            Cmax(c,idm(c))=1;
        end
        Cmaxtraj   = cat(3,Cmaxtraj,Cmax);
        % MI
        MI = [];
        for p = 1:9
            for pp = 1:9
            a = belief_sim(:,p);
            b = belief_fit(:,pp);
            mi=condentropy(a,b);
            MI(p,pp)=mi;
            end
        end
        MItraj   = cat(3,MItraj,MI);
    end
    disp(s)
end
