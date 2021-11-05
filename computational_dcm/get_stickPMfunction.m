function SF = get_stickPMfunction(options,condition)    
    
SF      = [];
ncond   = length(options.input_combo);

for c = 1:ncond
    
    coidx   = find(ismember(condition.condidx,options.input_combo{c}));
    CO      = condition.onset(coidx); % onset
   
    
    %-Create stimulus functions (32 bin offset) - from spm_get_ons
    %======================================================================
    k         = options.k;
    T         = options.T;
    dt        = options.dt;
    TR        = 1; % unit is in secs.
    ons       = CO;
    
    if ~strcmp(options.parametric_modulator{c},'no')
        u         = eval(sprintf('condition.%s(coidx)',options.parametric_modulator{c}))';
    else
        u         = ons.^0;
    end
    
    dur       = repmat(options.event_dur,length(ons),1);
    
    %-And scale so sum(u*dt) = number of events, if event-related
    %----------------------------------------------------------------------
    if ~any(dur)
        u     = u/options.dt;
    end
    
    
    ton       = round(ons*TR/dt) + 33;               % onsets
    tof       = round(dur*TR/dt) + ton + 1;          % offset
    sf        = sparse((k*T + 128),size(u,2));
    ton       = max(ton,1);
    tof       = max(tof,1);
    for j = 1:length(ton)
        if size(sf,1) > ton(j)
            sf(ton(j),:) = sf(ton(j),:) + u(j,:);
        end
        if size(sf,1) > tof(j)
            sf(tof(j),:) = sf(tof(j),:) - u(j,:);
        end
    end
    sf        = cumsum(sf);                         % integrate
    sf        = sf(1:(k*T + 32),:);
    
    SF(:,c)   = sf;
end