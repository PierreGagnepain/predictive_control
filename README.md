This is the Matlab code used for our paper on predictive memory control in PTSD, published in Nature Communications, by Leone, Postel, Mary, Fraisse, Vallee, Viader, de La Sayette, Peschanski, Dayan, Eustache, Gagnepain*

*corresponding author: pierre.gagnepain@inserm.fr

This code is composed of 4 folders:

===================
SIMULATION folder
-------------------

- MODEL FALSIFICATION

The goal of these simulations was to establish whether our computational models (please see our paper for model description) were able to generate the behavioral reduction in intrusion proportion that we normally observe across the four blocks of the TNT task (see Fig. 2).
We designed a virtual experimental setting with 144 suppression cues distributed across 4 TNT sessions, as in our real experiment. We started with a belief of .5 for the first trial, and at each new simulated trial, we generated a new belief based on the perceptual model considered and randomly drawn corresponding perceptual parameters. A suppression parameter was introduced to simulate memory suppression and to avoid the tilting of belief trajectories toward 1.

1) To launch the model falsification, please see /simulation/core/model_falsification.m.
2) The outcomes of model falsification can be downloaded from: https://drive.google.com/file/d/10t0vcP4WEWliBXBhM_fmBFpj69c7Fk24/view?usp=sharing
   
This downloaded file must then be stored in /simulation/store/simulation_belief2intrusion_beta_vfinale_v1.mat. This file contains:

	- the 200 virtual simulated participants (each of this virtual participant is encoded into a cell structure in simulation_belief2intrusion_beta_vfinale_v1) and repeated the virtual experiment 100 times using perceptual parameter randomly drawn from a Gaussian priors distribution tailored to match our own data (to sample plausible parameters)

  	- Binary rating generated for each of these 200 simulations were averaged across repeated sampling and summarized as intrusion proportion across the 4 artificial TNT sessions.

Please note that if you save your rerun "model_falsification" and save your own simulation_belief2intrusion_beta_vfinale_v1.mat, then you will also have to relaunch all parameter and model recovery analyses (see below).

3) We then compared these simulated intrusion proportion with real data. Comparaison plot used in figure 2 of our paper can be reproduced using /simulation/pipeline_fig2/figure2.m

- PARAMETER & MODEL RECOVERY

This section tests the ability of our models to recover generated beliefs (trajectory recovery), and their perceptual parameters (parameter recovery), and verify the reliability of the model selection criterion for identifying the true generative model within a set of competitive models, ensuring that this selection is not biased in favor of one particular model (model recovery).
Analyses of recovery tests the generative performance of a model, by verifying whether the fitting procedure produces meaningful trajectories and parameters, namely the true trajectories and parameters used to generate the data.

1) To launch parameter and model recovery, please see /simulation/core/parameter_and_model_recovery.m.


2) The outcomes of these analyses is stored in /simulation/result_simulation (already stored).
	- Each of the 200 virtual participants are stored as "simulation_belief2intrusion_beta_vfinale_v1_onestep+traj_XX" with XX being the virtual participant index
	- This file contains a structure name "fittedsim" with estimated perceptual parameter and Log-model evidence (Accuracy) for simulated and fitted models.
	- Please see /simulation/pipeline_fig2/figure2.m to plot the outcome of this analysis
	- Please note that each of the 200 simulation takes time to run and should not be launch on the same matlab. We use this code with our computation grid and queue submission system to launch the 200 jobs in parallel.

========================
COMPUTATIONAL DCM folder
------------------------

This folder contains the code to run computational DCM (requires SPM12 on your matlab path !!!)

1) The main starting code is start_computational_dcm.m.

2) This function requires your own data, but we provide one exemple in "REMEMBEREX001" folder.

3) This function will:
	- Define the DCM voi in native space, using the MNI coordinates defined in "roicoordinate.mat"
	- Estimate the belief trajectories from intrusion rating using the specified perceptual model (here tapas hgf 2 levels)
	- Run computational-DCM, using these trajectories as parametric modulator of the DCM modulatory stick function

This folder also contains a "plot_figures" folder with minimum dataset and code to reproduce figures 4-6

==========================
MODEL_FIT folder
--------------------------

Contains the main function to estimate "State", "Item", and "Combined" belief trajectories using intrusion rating given a perceptual and an observation model (tapas_binary_combined_onestep.m)


==========================
TAPAS_DEPENDENCIES
--------------------------
tapas functions used during model configuration and estimation
