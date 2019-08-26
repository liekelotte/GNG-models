# GNG-models
All files needed to perform computational modelling analysis of the valenced go/no-go task described first in Guitart-Masip et al., (2012), NeuroImage: "Go and no-go learning in reward and punishment: Interactions between affect"

run_all_models_ite.m is the script to run for the modelling analysis. The analysis requires the input of a datafile with datafields A (actions), S (states) and R (rewards). Each cell in these datafiles should contain the actions, states and rewards in order for each participant. The models that are fitted are specified in the run_all_models_ite.m file. The modelling file titles will tell you which models will be run:

a = alpha (learning rate)
b = beta (temp parameter)
2b = 2 betas (separate reward and punishment temperature parameters)
x = noise (lapse rate)
b = go bias (up the probability of performing a go from trial 1)
kwins = kappa for wins (boost learning rate on rewarded go trials only)
k = kappa for rewarded gos and punished no-gos (faster learning for rewarded go, slower learning for punished nogo)
ep = pavlovian parameter that fluctuates with the pavlovian value of each stimulus
epc = pavlovian parameter that is constant from the first time a reward or punishment is encountered in a state. will modulate 1 and -1 for rewards and punishments respectively 
