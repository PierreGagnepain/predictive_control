function [fitted, simulated, sim_intrusions, CONFacc,CONFcorr] = extract_simulations(nsim,nsampling,nitems,...
    simfile,respath,fittedsim_name,PerceptMod,ObsMod,parameters)

load(simfile)
fitted          = struct();
simulated       = struct();
sim_intrusions  = {};
lme             = [];
CONFacc         = [];
CONFcorr        = [];
model_indx      = [1:3;4:6;7:9];
i = 0;
for s = 1:nsim
    
    clear fn
    fn = fullfile(respath, sprintf([fittedsim_name '_%d.mat'],s));
    
    fprintf('Extracting simulation %d\n',s);
    
    if exist(fn)
        load(fn)
        i = i + 1;
        
        % estimated parameter
        estimated_par       = fittedsim.estimated_par;
        
        
        % perceptual model loop
        for m = 1:length(PerceptMod)
            
            % simulated belief
            simstate.(PerceptMod{m})  = simulation{i}.(PerceptMod{m}).beliefstate_sim;
            simitem.(PerceptMod{m})   = simulation{i}.(PerceptMod{m}).beliefitem_sim;
            simcomb.(PerceptMod{m})   = simulation{i}.(PerceptMod{m}).beliefcomb_sim;
            
            Cns = [];
            for ns = 1:nsampling
                
                % Extract simulated and fitted hidden parameters
                %-----------------------------------------------
                for par = 1:length(parameters{m})
                    
                    estpar = estimated_par{m}.(parameters{m}{par});%(:,:,3);
                    
                    if m < 3
                        
                        fitted.(PerceptMod{m}).state.(parameters{m}{par})(i,ns)     = estpar(ns,1,1);
                        fitted.(PerceptMod{m}).item.(parameters{m}{par})(i,ns)      = estpar(ns,2,2);
                        
                        state                       = simulation{s}.(PerceptMod{m}).sampled_parameter.state{1}{1}(ns);
                        items                       = reshape(cell2mat(vertcat(simulation{s}.(PerceptMod{m}).sampled_parameter.item{:})),[],nitems);
                        items                       = items(1:nsampling,:);
                        
                        
                        itemsavg                    = mean(items(ns,:));
                        simulated.(PerceptMod{m}).state.(parameters{m}{par})(i,ns)   = state;
                        simulated.(PerceptMod{m}).item.(parameters{m}{par})(i,ns)   = itemsavg;
                        
                    elseif m == 3
                        % "state" parameters
                        fitted.(PerceptMod{m}).state.(parameters{m}{par})(i,ns)      = estpar(ns,1);
                        fitted.(PerceptMod{m}).item.(parameters{m}{par})(i,ns)      = estpar(ns,2);
                        
                        state                       = simulation{s}.(PerceptMod{m}).sampled_parameter.state{1}{1}(ns);
                        items = [];
                        for it = 1:nitems
                            items(:,it)             = simulation{i}.(PerceptMod{m}).sampled_parameter.item{it}{1};
                        end
                        itemsavg                    = mean(items(ns,:));
                        
                        itemsavg                    = mean(items(ns,:));
                        simulated.(PerceptMod{m}).state.(parameters{m}{par})(i,ns)  = state;
                        simulated.(PerceptMod{m}).item.(parameters{m}{par})(i,ns)   = itemsavg;
                        
                        
                        
                    end
                end
                
                % Extract Simulated intrusions
                %-----------------------------------------------
                for ob = 1:length(ObsMod)
                    
                    sim_intrusions{s}.(PerceptMod{m}).(ObsMod{ob}) = simulation{s}.(PerceptMod{m}).(['sesavg_' ObsMod{ob}]);
                end
                
                % Correlation fitted/simulated belief trajectory
                %-----------------------------------------------
                
                % fit belief
                state   = [fittedsim.belief_fit.(PerceptMod{m}).state_state(ns,:)',fittedsim.belief_fit.(PerceptMod{m}).state_item(ns,:)',fittedsim.belief_fit.(PerceptMod{m}).state_comb(ns,:)'];
                item    = [fittedsim.belief_fit.(PerceptMod{m}).item_state(ns,:)',fittedsim.belief_fit.(PerceptMod{m}).item_item(ns,:)',fittedsim.belief_fit.(PerceptMod{m}).item_comb(ns,:)'];
                comb    = [fittedsim.belief_fit.(PerceptMod{m}).comb_state(ns,:)',fittedsim.belief_fit.(PerceptMod{m}).comb_item(ns,:)',fittedsim.belief_fit.(PerceptMod{m}).comb_comb(ns,:)'];
                
                % corr
                C      = [];
                C(1,:) = corr(simstate.(PerceptMod{m})(:,ns),state,'type','pearson');
                C(2,:) = corr(simitem.(PerceptMod{m})(:,ns),item,'type','pearson');
                C(3,:) = corr(simcomb.(PerceptMod{m})(:,ns),comb,'type','pearson');
                
                
                Cns(:,:,ns) = C;
                
                
            end
            Cb          = [];
            Cb          = mean(Cns,3);
            
            % Model accuracy
            %-----------------------------------------------
            acc              = fittedsim.ACC_simulation(model_indx(m,:),model_indx(m,:),:);
            lme.(PerceptMod{m})(:,:,i) = sum(acc,3);
            
            
            % Create Confusion matrices
            %-----------------------------------------------
            conf        = zeros(3,3);
            conf_corr   = zeros(3,3);
            for j=1:3
                
                % acc
                l = lme.(PerceptMod{m})(j,:,i);
                [jk,mx]        = max(l,[],2);
                conf(j,mx)     = 1;
                
                % belief
                b = Cb(j,:);
                [jk,mx]        = max(b,[],2);
                conf_corr(j,mx)     = 1;
            end
            CONFacc.(PerceptMod{m})(:,:,i)       = conf;
            CONFcorr.(PerceptMod{m})(:,:,i)      = conf_corr;
        end
    end
end

