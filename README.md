# AudioModellingStimuli
Tools for creating, processing and analyzing audio stimuli for auditory modeling

--> run_model.m pulls the functions together and creates the neurograms  
--> generate_stimulus_script.m is there to play around and you can listen to the stimuli  
--> generate_stimulus.m is called in run_model.m  

TODO (30/08/2024): 
- comment run_model.m and run_modelv2.m (Sascha will do soon)
- implement IOU for mean rate plots (add legend in plot & append to csv data file)
- add sum of spikes to csv file
- compute NSI/SSIM accross layers (AN/CN)
- check NSI/SSIM reliability/validity/accuracy
- check oscilatory behaviour in FD task
- Problem: noise (Pa is too high dB/Hz)
- vcv-task/noise (ICRA noise dBA)
- no off-ramp so far, implementation with amt_fade
- test run all tasks from battery
- change line colour of SSI/NSIM plot in plot_neuorgram_comparison.m
- human derived decision layer (logistic regression to IMAP)
