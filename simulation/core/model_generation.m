function output = model_generation(input,options,parname,prc_fun,meansup)

% prc_fun   = perceptual model
% meansup   = suppression factor for states (1) and items (2)

% ==========================================
% UNPACK INPUT & OPTIONS
% ==========================================
fieldinput = fieldnames(input);
for f = 1:length(fieldinput)
    eval(sprintf('%s = input.%s;',fieldinput{f},fieldinput{f}))
end
fieldoptions = fieldnames(options);
for f = 1:length(fieldoptions)
    eval(sprintf('%s = options.%s;',fieldoptions{f},fieldoptions{f}))
end


% ==========================================
% SIMULATE TRAJECTORIES
% ==========================================
sampled_parameter   = [];

% Simulate State trajectory
% -------------------------
% nu = 1+randn(1,1)*sqrt(.08);
[beliefstate_sim, Usim_state,precision_state,parsampling] = prc_fun(nitems,nreps,meansup(1),nsampling,parname,sdintrusion,nu,startprob,suppressiontype);
sampled_parameter.state{1} = parsampling;

% Simulate Item trajectory
% -------------------------
beliefitem_sim = [];precision_item = [];Usim_item = [];
for i = 1:nitems
%     nu = 1+randn(1,1)*sqrt(.08);
    [beliefitem, Usim,precision,parsampling] = prc_fun(1,nreps,meansup(2),nsampling,parname,sdintrusion,nu,startprob,suppressiontype);
    
    beliefitem_sim(find(itempos==i),:)  = beliefitem;
    precision_item(find(itempos==i),:)  = precision;
    Usim_item(find(itempos==i),:)       = Usim;
    
    sampled_parameter.item{i} = parsampling;
end

% Combined belief
% -------------------------
numerator       = (beliefstate_sim.*precision_state) + (beliefitem_sim.*precision_item); % EQ 8&9 Weilnhammer et al.2018
denominator     = precision_state + precision_item;
beliefcomb_sim  = [];
beliefcomb_sim  = numerator./denominator;

% replace first presentation with state
beliefcomb_sim(1:nitems,:) = beliefstate_sim(1:nitems,:);

% get beliefcomb_sim response
Usim_comb = [];suppressionfactor=meansup(3);%no suppression of combined ?
for b = 1:size(beliefcomb_sim,2)
    y = [];
    for bb = 1:size(beliefcomb_sim,1)
        bel             = beliefcomb_sim(1:bb,b);
        suppressedbel   = supfun(bel(end),suppressionfactor);
        suppressedbel   = max([0.1 suppressedbel]); % limit suppressed belief
        bel(end)        = suppressedbel;
        [y(bb)]         = respsim(bel,nu,sdintrusion,y,nitems);
        Usim_comb(b,bb) = y(bb);
    end
end
Usim_comb = Usim_comb';

% Get intrusion patterns across session
sesavg_state    = [];sesavg_item = [];sesavg_comb = [];
for b = 1:size(Usim_comb,2)
    int_pattern = Usim_state(:,b);
    sesavg_state(b,:) = mean(reshape(int_pattern,[],nsess),1);
    
    int_pattern = Usim_item(:,b);
    sesavg_item(b,:) = mean(reshape(int_pattern,[],nsess),1);
    
    int_pattern = Usim_comb(:,b);
    sesavg_comb(b,:) = mean(reshape(int_pattern,[],nsess),1);
end

% OUTPUT
% ======================
output = [];

% belief trajectory
output.beliefstate_sim = beliefstate_sim;
output.beliefitem_sim = beliefitem_sim;
output.beliefcomb_sim = beliefcomb_sim;

% Response (intrusion/no-intrusion)
output.Usim_state = Usim_state;
output.Usim_item = Usim_item;
output.Usim_comb = Usim_comb;

% 4*TNT sessions intrusion proportions
output.sesavg_state = sesavg_state;
output.sesavg_item = sesavg_item;
output.sesavg_comb = sesavg_comb;

% sampled parameters
output.sampled_parameter = sampled_parameter;

