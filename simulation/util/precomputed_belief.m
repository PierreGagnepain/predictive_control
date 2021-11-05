function [traj, infStates] = precomputed_belief(r, p, varargin)

% pseudo-function to force using pre-computed belief as output of the
% perceptual model during optimization
traj        = [];
infStates   = r.belief;