function [traj, infStates] = tapas_kf(r, p, varargin)
% The scalar Kalman filter
%
% This function can be called in two ways:
% 
% (1) tapas_kf(r, p)
%   
%     where r is the structure generated by tapas_fitModel and p is the parameter vector in native space;
%
% (2) tapas_kf(r, ptrans, 'trans')
% 
%     where r is the structure generated by tapas_fitModel, ptrans is the parameter vector in
%     transformed space, and 'trans' is a flag indicating this.
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2016 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Transform paramaters back to their native space if needed
if ~isempty(varargin) && strcmp(varargin{1},'trans');
    p = tapas_kf_transp(r, p);
end

% Unpack parameters
g_0   = p(1);      % Initial gain
mu_0  = p(2);      % Initial hidden state mean
expom = exp(p(3)); % Process variance
pi_u  = p(4);      % Observation precision

% Add dummy "zeroth" trial
u = [0; r.u(:,1)];
n = length(u);

% Initialize updated quantities
da = NaN(n,1); % Prediction error
g  = NaN(n,1); % Kalman gain
mu = NaN(n,1); % Hidden state mean

% Priors
g(1)  = g_0;
mu(1) = mu_0;

% Pass through update loop
for k = 2:1:n
    if not(ismember(k-1, r.ign))
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Effect of input u(k)
        %%%%%%%%%%%%%%%%%%%%%%
        
        % Prediction error
        da(k) = u(k)-mu(k-1);
        
        % Gain update
        g(k) = (g(k-1) +pi_u*expom)/(g(k-1) +pi_u*expom +1);
        
        % Hidden state mean update
        mu(k) = mu(k-1)+g(k)*da(k);
    else
        da(k) = 0;
        g(k)  = g(k-1);
        mu(k) = mu(k-1);
    end
end

% Predicted value
muhat = mu;
muhat(end) = [];

% shrink belief (to avoid pb with observation model)
muhat(muhat>0.95)=0.95;
muhat(muhat<0.05)=0.05;

% Remove priors
da(1) = [];
g(1)  = [];
mu(1) = [];

% Create result data structure
traj = struct;

traj.g     = g;
traj.muhat = muhat;
traj.mu    = mu;
traj.da    = da;

% Create matrix (in this case: vector) needed by observation model
infStates = [traj.muhat, traj.mu];

return;
