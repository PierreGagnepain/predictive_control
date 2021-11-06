function B = estimate_belief(condition,options)

% prepare data
nt_idx  = 1:length(condition.y);%find(~isnan(condition.y));% index of NT items
r       = [];
r.u     = [condition.y(nt_idx),condition.item_identity];
r.y     = condition.y(nt_idx);
r.ign   = [];
r.irr   = [];
np      = []; % in case want use different prior than those configured in tapas_config

% run tapas
[rstate,ritem,rcomb] = tapas_binary_combined_onestep(r,options.perceptualmodel,options.obsmodel,np,0);

B = rcomb.belief; % State/Item combined belief