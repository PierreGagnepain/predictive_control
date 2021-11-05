function computational_dcm(options,xY,condition,P)

% =========================================================================
% DCM
% =========================================================================
clear DCM
clear DCM_initial

% define DCM
DCM_initial.n = length(xY);      % number of regions
DCM_initial.v = size(xY(1).u,1); % number of time points
DCM_initial.Y.dt  = options.RT;
DCM_initial.Y.X0  = xY(1).X0; 

% reoder voi region according to DCM matrices (see define_dcmmatrix.m)
VOIidx              = [4 5 1 2 3]; % aMFG pMFG rHIP cHIP % PC

for i = 1:DCM_initial.n
    
    DCM_initial.Y.y(:,i)  = xY(VOIidx(i)).u;
    DCM_initial.Y.name{i} = xY(VOIidx(i)).name;
end

DCM_initial.Y.Q    = spm_Ce(ones(1,DCM_initial.n)*DCM_initial.v);

DCM_initial.delays = repmat(options.RT/2,DCM_initial.n,1);
DCM_initial.TE     = 0.03;

% Estimate dcm (loop through models)
% ----------------------------------------------------------------
nmodel      = length(options.dcm2run);
modindex    = options.dcm2run;

for dcmI = 1:nmodel
    
    clear DCM
    DCM = DCM_initial;
    
    % model index
    mI  = modindex(dcmI);
    
    % input
    DCM.U.dt   = options.dt;
    DCM.U.name = {};
    DCM.U.u    = [];
    
    % get stick function with parametric modulation
    if dcmI<29 % Computation & Bottom-up FAMILLIES
        options.parametric_modulator= {'no','no','PE','Belief'};% which DCM condition is parametrically modulated by compution ?
    elseif dcmI>28 % NO-Computation
        options.parametric_modulator= {'no','no','no','no'};
    end
    SF                      = get_stickPMfunction(options,condition);    
    
    % DCM condition names
    for nI = 1:length(options.input_combo)
        
        DCM.U.name{end + 1} = options.dcm_condlabel{nI};
    end
    
    % DCM condition stick functions
    for nI = 1:length(options.input_combo)
        
        DCM.U.u(:,nI)       = SF(33:end,nI); % to account for 32 bin offset (see get_stickPMfunction)
        
    end
    
    % DCM matrix
    DCM.a = [];
    DCM.a = options.A(:,:,mI);
    
    DCM.b = [];
    DCM.b = options.B{mI};
    
    DCM.c = [];
    DCM.c = options.C(:,:,mI);
    
  
    % dcm option
    DCM.options.nonlinear     = 0;
    DCM.options.two_state     = 0;
    DCM.options.stochastic    = 0;
    DCM.options.centre        = 0;
    DCM.options.nograph       = 1;
    DCM.options.endogenous    = 0;
    
    % estimate
    DCM = spm_dcm_estimate(DCM);
    
    dcm_name = sprintf('%s_DCM_%s',P.subcode,options.dcm_name);
    save(fullfile(P.dcm_dir,sprintf([dcm_name,'_mod%0d'],mI)),'DCM');
    
end