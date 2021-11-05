function [parname,modname] = get_parname()


% HGF
% ------------------------------------------
parname                         = {};
parname{1}.name                 = {'om'}; % omega
parname{1}.idx                  = 13; % paramater indexes in tapas priors config

% parameters for simulations are sampled from a gaussian distribution with prior_mean and prior_std (based on data)
% first rows is for state, second for items
parname{1}.prior_mean(1,1)      = -7.05; % state mean distribution
parname{1}.prior_var(1,1)       = 1.8;  % state std distribution
parname{1}.prior_mean(2,1)      = -1.5; % item mean distribution
parname{1}.prior_var(2,1)       = 0.74; % item std distribution

% if sampling from uniform (not used)
parname{1}.prior_min(1,1)      = -8.7; % state min
parname{1}.prior_max(1,1)      = -3; % state max
parname{1}.prior_min(2,1)      = -2.9; % item min
parname{1}.prior_max(2,1)      = -0.19; % item max

% RW
% ------------------------------------------
parname{2}.name = {'al'}; % alpha
parname{2}.idx  = 2;
parname{2}.prior_mean(1,1)   = -3.05; 
parname{2}.prior_var(1,1)    = 4.4; 
parname{2}.prior_mean(2,1)   = -1.2;
parname{2}.prior_var(2,1)    = 0.015;

% if sampling from uniform (not used)
parname{2}.prior_min(1,1)      = -4.26;
parname{2}.prior_max(1,1)      = 0;
parname{2}.prior_min(2,1)      = -1.16;
parname{2}.prior_max(2,1)      = -0.79;


% Kalman filter
% ------------------------------------------
parname{3}.name                 = {'expom','pi_u'};
parname{3}.idx                  = [3 4];
parname{3}.prior_mean(1,:)      = [-3.0335 0.4938]; 
parname{3}.prior_var(1,:)       = [0.0045  0.0001]; 
parname{3}.prior_mean(2,:)      = [-3.0566 0.47];
parname{3}.prior_var(2,:)       = [0.0020 0.0007];

% if sampling from uniform (not used)
parname{3}.prior_min(1,:)       = [-3.1148 0.49]; 
parname{3}.prior_max(1,:)       = [-2.9957 0.50]; 
parname{3}.prior_min(2,:)       = [-3.1383 0.427];
parname{3}.prior_max(2,:)       = [-2.9976 .5014];

% model name (for storage)
modname = {'HGF','RW','KF'};