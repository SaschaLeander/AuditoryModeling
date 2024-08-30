%% Script runs stimulus, generation, zilany2014, carney2015 and stores neurogram similarity ==================================================================

% code written by Sam Jones and Sascha Muhlinghaus

%% set path, start modified amt ===========================================

addpath 'amtoolbox-full-1.5.0/amtoolbox-full-1.5.0'

amt_start;

amt_mex;

%% sanity check - ensure modified models work as expected

% demo_zilany2014; % run and compare to link below 
% https://amtoolbox.org/amt-1.5.0/doc/demos/demo_zilany2014.php 

% demo_carney2015; % run and compare to link below
% https://amtoolbox.org/amt-1.5.0/doc/demos/demo_carney2015.php 

%% system configuration ===================================================
function [nsimArray, ssimArray, varArray] = run_analysis(numCF_input, condition_input, dependent_var_input, limit_input, stepper_input, status_input, varargin)

    Fs = 100e3; % "model sampling frequency" [samples/sec]
    
    %% stimulus params =====================================================
    condition = condition_input;  % condition_flag: 'forward masking', 'backward masking', 'simultaneous', 'notch', 'frequency discrimination'
    signal_length_ms = 20; % signal_length_ms: length of the stimulus in (ms)
    switch condition
        case 'frequency discrimination'
            dB_signal = 70; % dB_signal: desired decibel level of the tone
            signal_freq = dependent_var_input; % tone_freq: frequency of the tone in Hz
    
        otherwise
            dB_signal = dependent_var_input; % dB_signal: desired decibel level of the tone
            signal_freq = 1000; % tone_freq: frequency of the tone in Hz
    end
    
    db_noise = 30; % varargin1: inputs depending on the chosen condition
    gap = 50; % varargin2: inputs depending on the chosen condition
    vcv_string = 'a b a';
    debug_flag = false;
    
    %% generate a neurogram for a specific hearing loss =======================
    
    ag_dbloss = [0 0 0 0 0 0 0]; % modify loss here (+ db, linked to ag_fs)
    ag_fs = [125 250 500 1e3 2e3 4e3 8e3]; % audiometric frequencies
    
    %% define auditory pathway params (AN, CN, IC) ========================
    
    species = 2; % 1 = cat; 2 = human
    numCF = numCF_input; % number of inner hair cells
    CF_range = [125 14e3]; % frequency range
    cf = audspace(CF_range(1),CF_range(2),numCF); % generate cfs
      
    fiberType = 1; % AN fiber type, 1 = L, 2 = M, 3 = H spontaneous rate, 4 = ratio defined below
    fiber_num = 20; % number of fibres
    numH = 10; % high spontaneous rate fibres
    numM = 5; % medium spontaneous rate fibres
    numL = 5; % low spontaneous rate fibres
    BMF = 100; % default CN/IN best modulation frequency
    
    %% run models and save data =============================================================
    limit = limit_input;
    stepper = stepper_input;
    
    % Which metric is computed
    flag = 'both';
    
    % Number of iterations for plotting neurogram
    plot_every = 20;
    
    % Dependent variable
    dependent_var = dependent_var_input; % Either dB_signal or signal_freq
    
    % Number of generated samples
    samples = abs(abs(dependent_var - limit) / stepper);
    
    fprintf('Number of samples: %d\n', samples);
    
    % Initialize arrays for similarity function
    nsimArray = [];
    ssimArray = [];
    
    varArray = [];
    
    % Choose degraded or reference condition based on input
    if status_input == 2
        fprintf('Status: %s\n', 'degraded');
        bm_deg = varargin{1}; % scale basilar membrane tuning, 1 = normal, > 1 = broader
        cohc_deg = varargin{2}; % OHC function scaling factor: 1 denotes normal function
        cihc_deg = varargin{3}; % IHC function scaling factor: 1 denotes normal function
        bm_current = bm_deg;
        cohc_current = cohc_deg;
        cihc_current = cihc_deg;
    else
        fprintf('Status: %s\n', 'control');
        bm_current = 1;
        cohc_current = 1;
        cihc_current = 1;
    end
    
    % Generate stimulus 1 (reference)
    switch condition
        case 'vcv'
            [stimulus1, ~, ~] = generate_stimulus(condition, ...
                signal_freq, signal_length_ms, dependent_var, debug_flag, vcv_string);
        case 'frequency discrimination'
            [stimulus1, ~, ~] = generate_stimulus(condition, ...
                dependent_var, signal_length_ms, dB_signal, debug_flag, db_noise, gap);
        otherwise
            [stimulus1, ~, ~] = generate_stimulus(condition, ...
                signal_freq, signal_length_ms, dependent_var, debug_flag, db_noise, gap);
    end
    
    % Model auditory response
    [r_mean1, ihc1, ic_sout_BE1, cn_sout_contra1] = simulate_auditory_response(stimulus1, ...
        Fs, cf, bm_current, fiberType, numH, numM, numL, fiber_num, cohc_current, cihc_current, BMF);
    
    % Plot reference neurogram
    plot_neurograms(stimulus1, r_mean1, ihc1, ic_sout_BE1, cn_sout_contra1, Fs, cf);
    
    for i = 1:samples
        
        % Update dependent variable 
        dependent_var = dependent_var + stepper;
        
        fprintf('dependent variable (dB/Hz): %d\n', dependent_var);
        varArray = [varArray; dependent_var];
           
        % Generate stimulus 2 (altered)
        switch condition
            case 'vcv'
                [stimulus2, ~, ~] = generate_stimulus(condition, ...
                    signal_freq, signal_length_ms, dependent_var, debug_flag, vcv_string);
            case 'frequency discrimination'
                [stimulus2, ~, ~] = generate_stimulus(condition, ...
                    dependent_var, signal_length_ms, dB_signal, debug_flag, db_noise, gap);
            otherwise
                [stimulus2, ~, ~] = generate_stimulus(condition, ...
                    signal_freq, signal_length_ms, dependent_var, debug_flag, db_noise, gap);
        end
        
        % Model auditory response for altered stimulus
        [r_mean2, ihc2, ic_sout_BE2, cn_sout_contra2] = simulate_auditory_response(stimulus2, ...
            Fs, cf, bm_current, fiberType, numH, numM, numL, fiber_num, cohc_current, cihc_current, BMF);
        
        if mod(i, plot_every) == 0
            % Plot neurogram 2 (altered)
            %plot_neurograms(stimulus2, r_mean2, ihc2, ic_sout_BE2, cn_sout_contra2, Fs, cf);
        end
        
        % Compute similarity scores
        [nsim_score, ssim_score] = compute_similarity(ic_sout_BE1, ic_sout_BE2, flag);
        nsimArray = [nsimArray; nsim_score];
        ssimArray = [ssimArray; ssim_score];
        
        % Save parameters and data to csv file
        switch condition
            case 'frequency discrimination'
                manage_data(i, dB_signal, dependent_var, cihc_current, cohc_current, ...
                    bm_current, condition, dependent_var_input, nsim_score, ssim_score);
            otherwise
                manage_data(i, dependent_var, signal_freq, cihc_current, cohc_current, ...
                    bm_current, condition, dependent_var_input, nsim_score, ssim_score);
        end
    end  
    % Plot comparison neurogram plot 
    plot_neurogram_comparison(stimulus1, r_mean1, ihc1, ic_sout_BE1, cn_sout_contra1, ...
        Fs, cf, stimulus2, r_mean2, ihc2, ic_sout_BE2, cn_sout_contra2, ...
        varArray, nsimArray, ssimArray);
end


%% Call function ==========================================================
% bm_params = [3,1,1,3];
% cihc_params = [1,0.3,1,0.3];
% cohc_params = [1,1,0.3,0.3];
% 
% 
% titles = {
%     ['Neurogram Similarity, numCF=', num_cf ,', BM = 3'], ...
%     ['Neurogram Similarity, numCF=', num_cf ,', cihc = .3'], ...
%     ['Neurogram Similarity, numCF=', num_cf ,', cohc = .3'], ...
%     ['Neurogram Similarity, numCF=', num_cf ,', full impairment bm = 3, cihc = .3, cohc = .3']};
% 
% condition_params = {'frequency discrimination', 'vcv'};
% 
% limit_dependent_var = [1500, 30];
% change_rate_dep_var = [10, -5];
% dependent_variable = [1000, 80];

function wrapper(bm, cihc, cohc, num_cf, condition, dependent_variable, limit_dependent_var, change_rate_dep_var)
    
    status = 1;
    
    [nsimArray_ref, ssimArray_ref, varArray_ref] = run_analysis(num_cf, ...
        condition, dependent_variable, limit_dependent_var, change_rate_dep_var, status);
    
    status = 2;

    [nsimArray_degraded, ssimArray_degraded, varArray_degraded] = run_analysis( num_cf, ...
        condition, dependent_variable, limit_dependent_var, change_rate_dep_var, status, bm, cihc, cohc);
end 

bm = 3;
cohc = 1;
cihc = 1;
num_cf = 50;
condition = 'frequency discrimination';
dependent_variable = 1000;
limit_dependent_var = 1500;
change_rate_dep_var = 10;
title_input = 'Neurogram Similarity, numCF=5, full impairment bm = 3, cihc = .3, cohc = .3';

wrapper(bm, cohc, cihc, num_cf, condition, dependent_variable, limit_dependent_var, change_rate_dep_var);

