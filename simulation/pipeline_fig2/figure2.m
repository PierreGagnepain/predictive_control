%% ========================================================================
%  DEFINE PATHS
%  ========================================================================
clear;clc;close all;
rootpath = '.../simulation/';addpath(genpath(rootpath))


%% ========================================================================
%  DEFINE INPUTS
%  ========================================================================
nsim        = 200; % Number of simulations;
nsampling   = 100; % Number of sampled parameters for each simulation;
nitems      = 18; % Number of word-object pairs (i.e., items);
nsess       = 4; % Number of sessions
PerceptMod  = {'HGF','RW','KF'}; % Perceptual models
ObsMod      = {'state','item','comb'}; % Observation models
parameters  = {{'om'}, {'al'}, {'expom' 'pi_u'}};


% Define paths and names
simfile              = fullfile(rootpath,'store',...
    sprintf('simulation_belief2intrusion_beta_vfinale_v1.mat')); % Simulations file
respath         = fullfile(rootpath,'result_simulation'); % Path of the models' fitting of simulated intrusions.
fittedsim_name  = 'simulation_belief2intrusion_beta_vfinale_v1_onestep+traj'; % Name of the models' fitting of simulated intrusions.



%% ========================================================================
%  EXTRACT SIMULATIONS
%  ========================================================================

[fitted, simulated, sim_intrusions, CONFacc,CONFcorr] = extract_simulations(nsim,nsampling,nitems,...
    simfile,respath,fittedsim_name,PerceptMod,ObsMod,parameters);

%% ========================================================================
%  PLOT EXAMPLE TRAJECTORIES - Fig. 2a
%  ========================================================================
% Example data (partecipant's subcode: REMEMBEREX019TM)

load('example_data')
plot_example_beliefs(example_data)

%% ========================================================================
%  MODEL FALSIFICATION PLOT - Figure 2b
%  ========================================================================

[Falsification, SimIntrusionPattern,real_sim_diff, real_sim_corr] = model_falsification_plot(...
    sim_intrusions,PerceptMod,ObsMod,nsim,nsampling,nsess);

%% ========================================================================
%  BELIEF/MODEL RECOVEY PLOT - Figure 2c
%  ========================================================================
modelrecovery_plot(CONFcorr,CONFacc,PerceptMod)

%% ========================================================================
%  PARAMETER RECOVERY - Figure 2d
%  ========================================================================
ObsMod       = {};
ObsMod       = {'state','item'};

[simfit_corr] = parameter_recovery_plot(simulated,fitted,PerceptMod,ObsMod,parameters,nsim,real_sim_diff, real_sim_corr);

