# AudioModellingStimuli
Tools for creating, processing and analyzing audio stimuli for auditory modeling

--> run_modelv2.m pulls the functions together and creates the neurograms  
--> generate_stimulus_script.m is there to play around and you can listen to the stimuli  
--> generate_stimulus.m, simulate_auditory_response.m, plot_neurograms.m, plot_neurogram_comparison.m and manage_data.m is called in run_modelv2.m  

TODO (01/09/2024): 
- implement IOU for mean rate plots (add legend in plot & append to csv data file)
- add sum of spikes to csv file
- compute NSI/SSIM across layers (AN/CN)
- check NSI/SSIM reliability/validity/accuracy
- check oscilatory behaviour in FD task
- Problem: noise (Pa is too high dB/Hz)
- vcv-task/noise (ICRA noise dBA)
- no off-ramp so far, implementation with amt_fade
- test run all tasks from battery
- change line colour of SSI/NSIM plot in plot_neuorgram_comparison.m
- human-derived decision layer (logistic regression to IMAP)
