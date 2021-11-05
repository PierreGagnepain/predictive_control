function start_computational_dcm(subcode)

% Requires spm12 on your matlab path !!!

%% Define path
% =========================================================================

% here subcode = 'REMEMBEREX001'

P               = [];
P.rootpath      = '/..../predictive_control/';
P.funpath       = fullfile(P.rootpath,'computational_dcm');
cd(P.rootpath); addpath(genpath(P.funpath))
P.subpath       = fullfile(P.funpath,subcode); % subject data directory 
P.fmripath      = fullfile(P.subpath,'tntfile');% preprocessed fmri files
P.glmpath       = fullfile(P.subpath,'spm_glm'); % path containing first-level analysis (SPM.mat) using preprocessed native space fmri data defined above
P.invfile       = get_files(P.subpath,'y_r*_anatT1.nii'); % forward deformation field from SPM normalization
P.dcm_dir       = fullfile(P.funpath,'dcm_output');% directory DCM storage
P.voi_dir       = fullfile(P.funpath,'voi_output');% directory VOI storage
P.cond_dir      = fullfile(P.funpath,'cond_output');% directory condition information storage
P.subcode       = subcode;

% create storage directories if necessary
try mkdir(P.dcm_dir);end
try mkdir(P.voi_dir);end
try mkdir(P.cond_dir);end

%% Define ROI 
% =========================================================================
% define your own, here is an exemple
load(fullfile(P.funpath,'roicoordinate')); % roicoordinate contain XYZ (in mm) coordinates in MNI space (as well as ROI labels)
nvoxel = 30; % N (native space) voxels use to create VOI

%% Contrast used to select VOI: name of SPM contrast (as defined in SPM.mat), define one per ROI
% =========================================================================
cname    = {};
cname{1} = 'TH-NT'; % 'rHip'
cname{2} = 'TH-NT'; % 'cHip'
cname{3} = 'TH-NT'; % 'PC'
cname{4} = 'NT-TH'; % 'aMFG'
cname{5} = 'NT-TH'; % 'pMFG'

%% CREATE VOI
% =========================================================================
analysiname     = 'DCM_native'; % name for saving
runcreatevoi    = 1;
if runcreatevoi
    xY = native_dcmvoi(P,roicoordinate,cname,nvoxel,analysiname);
else
    fn = fullfile(P.voi_dir,sprintf('%s_VOI_%s.mat',P.subcode,analysiname));
    load(fn,'xY')
end

%% RUN computational-DCM
% =========================================================================
options                     = [];

% SPM/Condition
load(fullfile(P.glmpath,'SPM.mat'))
options.RT                  = 2.05; % fMRI Repetion time
options.dt                  = SPM.Sess(1).U(1).dt; % micro-time resolution
options.sessons             = cumsum(SPM.nscan)-SPM.nscan(1); % to correct onset of trials given TNT sessions
options.sessons             = options.sessons.*options.RT; % session onset (in sec)
options.condlabel           = {'IN','NI'}; % intrusion and non-intrusion (cond #1 and #2 below)
options.event_dur           = 3; % duration in sec of reminder cue
options.k                   = sum(SPM.nscan); % total number of scans
options.T                   = SPM.xBF.T; % parameter for HFR function

% tapas
options.perceptualmodel     = 'tapas_hgf_binary_2lev_config'; % for belief estimation in TAPAS toolbox
options.obsmodel            = 'tapas_beta_obs_config'; % for belief estimation in TAPAS toolbox

% DCM
[A,B,C]                     = define_dcmmatrix(); % function to get DCM.A, DCM.B, DCM.C matrices (define your own !)
options.A                   = A; % DCM.A intrinsic connections
options.B                   = B; % DCM.B modularity input
options.C                   = C; % DCM.C driving input
options.input_combo         = {}; % Define DCM conditions
options.input_combo{1}      = [1];      % Intrusions (driving input #1)
options.input_combo{2}      = [1 2];    % No-Think (driving input #2)
options.input_combo{3}      = [1];      % Prediction-error (only on intrusion trial)
options.input_combo{4}      = [1 2];    % Belief (for both intrusive and non-intrusive cues)
options.dcm_condlabel       = {'IN','NT','PE','BELIEF'}; % condition labels
options.dcm_name            = 'computational_v1';
options.dcm2run             = 1:44;%1:size(A,3); % what DCM model to run ? (depends on define_dcmmatrix function)

% Get ONSET (in sec) for each CONDITION, and Intrusion Binary vector form SPM file
% --------------------------------------------------------------------------------
I = getIntrusion(SPM,options); % 1st col of "I" is onset and 2nd col is intrusion rating


% Define Item identity (use your own function)
% -----------------------------------------
tnt_file                    = fullfile(P.subpath,'tnt.txt');
[y,item_identity,SO]        = get_intrusion(tnt_file);% SO (session & onset)

% 1st col. of "item_identity" is a vector of N no-think trials which encodes item
% identity (i.e. I1, I2, I3, ... I18). !! USE YOUR OWN (example provided)

% 2nd col. of "item_identity" is the index of the TNT sessions

% store stimuli info
condition                   = [];
condition.onset             = I(:,1); % stim onset
condition.item_identity     = item_identity;
condition.y                 = I(:,2); % intrusion response

% Estimate Belief (tapas): here we use the "Combined" belief
% ----------------------------------------------------------
Belief              = estimate_belief(condition,options);
pe                  = condition.y-Belief; % prediction-error
pe(pe<0)            = 0;% PE+

% store belief/PE trajectory along with stimuli info
nt_idx              = 1:length(condition.y);% index of NT items
condition.Belief    = nan(1,length(condition.onset));
condition.Belief(nt_idx) =  Belief;
condition.PE        = nan(1,length(condition.onset));
condition.PE(nt_idx)=  pe;
condition.condidx   = condition.y;
condition.condidx(condition.condidx == 1)   = 1;% INT
condition.condidx(condition.condidx == 0)   = 2;% NI

% save .mat
fn = fullfile(P.cond_dir,sprintf('%s_condition.mat',P.subcode));
save(fn,'condition')

% run DCM
% -----------------------------------------
computational_dcm(options,xY,condition,P)

