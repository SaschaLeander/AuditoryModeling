%% Function simulates the acitivity of neural output from the auditory nerve, cochlear nucleus and inferior colliculus

% Author Sascha Muhlinghaus (05/08/2024)
% contact: saschamuhlinghaus@gmail.com

function [r_mean, ihc, ic_sout_BE, cn_sout] = simulate_auditory_response(stim, Fs, cf, bm, fiberType, numH, numM, numL, fiber_num, cohc, cihc, BMF)

    % SIMULATE_AUDITORY_RESPONSE return cochlear nucleus and inferior
    % colliculus output
    %
    %   Usage: [ic_sout_BE, psth, cn_sout_contra] = simulate_auditory_response(stim, Fs, cf, bm, fiberType, numH, numM, numL, fiber_num, cohc, cihc, BMF)
    %
    %   Input parameters:
    %       stim:       audio sample 
    %       Fs:         sampling frequency
    %       cf:         cochlear frequency range
    %       bm:         basilar membrane tuning, 1 = normal, > 1 = broader 
    %       fibertype:  fiber types 
    %       numH:       high spontaneous rate fibres
    %       numM:       medium spontaneous rate fibres
    %       numL:       low spontaneous rate fibres
    %       fiber_num:  number of auditory nerve fibres
    %       cohc:       OHC function scaling factor: 1 denotes normal function
    %       cihc:       IHC function scaling factor: 1 denotes normal function  
    %       BMF:        CN/IN best modulation frequency
    %
    %   Output parameters:
    %       r_mean:     mean firing rates
    %       ihc:        inner hair cell firing 
    %       ic_sout_BE: inferior colliculus output
    %       cn_sout:    cochlear nucleus
    
    % inner hair cells and synapse out   
    [r_mean, ~, ihc] = zilany2014(stim, Fs, cf, 'bm', bm, 'fiberType', fiberType, ...
        'numH', numH, 'numM', numM, 'numL', numL, 'nrep', fiber_num, ...
        'cohc', cohc, 'cihc', cihc);
    
    % cochlear nucleus and inferior colliculus out  
    [ic_sout_BE, ~, cn_sout] = carney2015(r_mean, BMF, Fs);
    
end

