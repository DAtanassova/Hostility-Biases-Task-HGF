# HostilityBiasTask_HGF
Repo for HGF model fitting for a hostility bias task with perceptually ambiguous outcomes.
The perceptual model "tapas_ehgf_ar1_binary_pu" is adapted from the original ehgf_ar1_binary model to handle perceptual ambiguity. The response model "tapas_logrt_linear_binary_val" is adapted from the original tapas_logrt_linear_binary and models trialwise RTs according to the following parameters:
- be1: captures the valence (effect of facial expression), increasing/decreasing RTs as a function on whether the displayed face was hostile (angry) or happy (non-hostile). Positive be1=slower RTs on hostile trials, negative be1=faster RTs on hostile trials.
- be2: captures the effect of surprise: increases/decreases RTs as a function of how unexpected the stimulus was. Positive be2=slower RTs on surprising trials.
- be3: captures the interaction between type of facial expression and surprise (tests whether surprise slows RTs  more  specifically when the condition is hostile, compared to when it isn't). Positive be3=surprise has a bigger RT-slowing effect in the hostile condition than elsewhere.
- zeta: decision noise. Smaller zeta=more deterministic RTs 

The original models and all supporting scripts can be accessed from the tapas repository (https://github.com/translationalneuromodeling/tapas/tree/master/HGF).

**Analysis Pipeline**
Running the script "HB_modelFit.m" will execute the model fitting, simulations and recovery, and will extract the fitted parameters and trialwise computational quantities (e.g. predictions muhat1, updates mu2, mu3, uncertainty in predictions sahat1, precision sa2, sa3, pwPEs, valence (hostile/non-hostile), trial number, RT). 

You need to specify the paths to the different models, and maintain a local version of the HGF supporting scripts (e.g. for functions "fitModel", "simModel", etc). The HGF models adapted for this task can be found in the subfolder "HB_HGF-models".
