%% about ==================================================================

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

function run_analysis(bm, cohc, cihc, title_text)

    Fs = 100e3; % "model sampling frequency" [samples/sec]
    
    %% stimulus parameter =====================================================
    condition = 'frequency discrimination';  % condition_flag: 'before', 'after', 'simultaneous', 'notch', 'frequency discrimination'
    signal_freq = 1000; % tone_freq: frequency of the tone in Hz
    signal_length_ms = 20; % signal_length_ms: length of the stimulus in (ms)
    dB_signal = 70; % dB_signal: desired decibel level of the tone
    
    db_noise = 30; % varargin1: inputs depending on the chosen condition 
    gap = 50; % varargin2: inputs depending on the chosen condition
    vcv_string = 'a b a';
    debug_flag = false;
    
    %% generate a neurogram for a specific hearing loss =======================
    
    ag_dbloss = [0 0 0 0 0 0 0]; % modify loss here (+ db, linked to ag_fs)
    ag_fs = [125 250 500 1e3 2e3 4e3 8e3]; % audiometric frequencies
    
    %% define auditory pathway parameters (AN, CN, IC) ========================
    
    species = 2; % 1 = cat; 2 = human
    numCF = 1000; % number of inner hair cells 
    CF_range = [125 14e3]; % frequency range
    cf = audspace(CF_range(1),CF_range(2),numCF); % generate cfs
    
    % scale basilar membrane tuning, 1 = normal, > 1 = broader
    bm_ref = 1; 
    bm_degraded = bm;
    % OHC function scaling factor: 1 denotes normal function
    cohc_ref = 1; 
    cohc_degraded = cohc; 
    % IHC function scaling factor: 1 denotes normal function
    cihc_ref = 1;
    cihc_degraded = cihc;
    
    
    fiberType = 1; % AN fiber type, 1 = L, 2 = M, 3 = H spontaneous rate, 4 = ratio defined below
    fiber_num = 20;% number of fibres
    numH = 10; % high spontaneous rate fibres
    numM = 5; % medium spontaneous rate fibres
    numL = 5; % low spontaneous rate fibres
    BMF = 100; % default CN/IN best modulation frequency   
    
    %% run models and save data =============================================================
    limit = 1500;
    stepper = 10; 
    
    % Which metric is computed
    flag = 'both';
    
    % Dependent variable 
    dependent_var = signal_freq; % Either dB_signal or signal_freq
    
    % Number of samples that will be generated
    samples = abs(dependent_var-limit) / stepper;
    
    % Initilize arrays for similarity function
    nsimArray_degraded = [];
    ssimArray_degraded = [];
    
    nsimArray_ref = [];
    ssimArray_ref = [];
    
    varArray = [];
    
    % Status refers to degraded or normal function
    status = 'REFERENCE';
    
    competitor_stimulus = signal_freq;
    
    for i = 1:samples
    
        if i == 1
            % Generate stimulus 1
            [stimulus1,Fs_stimulus1, length_stim1] = generate_stimulus(condition, ...
                signal_freq, signal_length_ms, dB_signal, debug_flag, db_noise, gap); 
    
            % Model auditory response 1a
            [r_mean1a, ihc1a, ic_sout_BE1a, cn_sout_contra1a] = simulate_auditory_response(stimulus1, ...
                Fs, cf, bm_ref, fiberType, numH, numM, numL, fiber_num, cohc_ref, cihc_ref, BMF);
            
            % Model auditory response 1b
            [r_mean1b, ihc1b, ic_sout_BE1b, cn_sout_contra1b] = simulate_auditory_response(stimulus1, ...
                Fs, cf, bm_ref, fiberType, numH, numM, numL, fiber_num, cohc_ref, cihc_ref, BMF);
    
            % Plot neurograms stim 1
            plot_neurograms(stimulus1, r_mean1a, ihc1a, ic_sout_BE1a, cn_sout_contra1a, Fs, cf);
            plot_neurograms(stimulus1, r_mean1b, ihc1b, ic_sout_BE1b, cn_sout_contra1b, Fs, cf);
        end
        
        % Independent variable 
        dependent_var = dependent_var + stepper;
        fprintf('Tone frequency: %d\n', dependent_var);
    
        % Generate stimulus 2
        [stimulus2,Fs_stimulus2, length_stim2] = generate_stimulus(condition, ...
            dependent_var, signal_length_ms, dB_signal, debug_flag, db_noise, gap); 
        
        % Model auditory response 2
        [r_mean2, ihc2, ic_sout_BE2, cn_sout_contra2] = simulate_auditory_response(stimulus2, ...
            Fs, cf, bm_ref, fiberType, numH, numM, numL, fiber_num, cohc_degraded, cihc_degraded, BMF);
        
        if mod(i,10) == 0
            % Plot neurogram 2
            plot_neurograms(stimulus2, r_mean2, ihc2, ic_sout_BE2, cn_sout_contra2, Fs, cf);
        end
    
       % [nsim_score_samdif, ssim_score_samdif] = compute_similarity(ic_sout_BE1a, ic_sout_BE2, flag);
       % [nsim_score_samsam, ssim_score_samsam] = compute_similarity(ic_sout_BE1a, ic_sout_BE1b, flag);
       % [nsim_score_difsam, ssim_score_difsam] = compute_similarity(ic_sout_BE2, ic_sout_BE1b, flag);
    
       fprintf('Status:  %s\n', status);
       % fprintf('NSI for stimulus 1&2: %d\n', nsim_score_samdif);
       % fprintf('NSI for stimulus 1&1: %d\n', nsim_score_samsam);
       % fprintf('NSI for stimulus 2&1: %d\n', nsim_score_difsam);

       [nsim_score_ref, ssim_score_ref] = compute_similarity(ic_sout_BE1a, ic_sout_BE2,flag);
       nsimArray_ref = [nsimArray_ref; nsim_score_ref]; 
       ssimArray_ref = [ssimArray_ref; ssim_score_ref]; 
    
       manage_data(i, dB_signal, dependent_var, cihc_ref, cohc_ref, bm_ref, condition, ...
           competitor_stimulus, nsim_score_ref, ssim_score_ref);
    end
    
    % Dependent variable 
    dependent_var = signal_freq; % Either dB_signal or signal_freq
    status = 'DEGRADED';
    
    for j = 1:samples
    
        if j == 1
            % Generate stimulus 1
            [stimulus1,Fs_stimulus1, length_stim1] = generate_stimulus(condition, ...
                signal_freq, signal_length_ms, dB_signal, debug_flag, db_noise, gap); 
        
            % Model auditory response 1a
            [r_mean1a, ihc1a, ic_sout_BE1a, cn_sout_contra1a] = simulate_auditory_response(stimulus1, ...
                Fs, cf, bm_degraded, fiberType, numH, numM, numL, fiber_num, cohc_degraded, cihc_degraded, BMF);
            
            % Model auditory response 1b
            [r_mean1b, ihc1b, ic_sout_BE1b, cn_sout_contra1b] = simulate_auditory_response(stimulus1, ...
                Fs, cf, bm_degraded, fiberType, numH, numM, numL, fiber_num, cohc_degraded, cihc_degraded, BMF);
    
           % Plot neurograms stim 1
            plot_neurograms(stimulus1, r_mean1a, ihc1a, ic_sout_BE1a, cn_sout_contra1a, Fs, cf);
            
            plot_neurograms(stimulus1, r_mean1b, ihc1b, ic_sout_BE1b, cn_sout_contra1b, Fs, cf);
        end
        
        % Independent variable 
        dependent_var = dependent_var + stepper;
    
        varArray = [varArray; dependent_var];
    
        % Generate stimulus 2
        [stimulus2,Fs_stimulus2, length_stim2] = generate_stimulus(condition, ...
            dependent_var, signal_length_ms, dB_signal, debug_flag, db_noise, gap); 
        
        % Model auditory response 2
        [r_mean2, ihc2, ic_sout_BE2, cn_sout_contra2] = simulate_auditory_response(stimulus2, ...
            Fs, cf, bm_degraded, fiberType, numH, numM, numL, fiber_num, cohc_ref, cihc_ref, BMF);
        
        if mod(i,20) == 0
            % Plot neurogram 2
            plot_neurograms(stimulus2, r_mean2, ihc2, ic_sout_BE2, cn_sout_contra2, Fs, cf);
        end
    
        % [nsim_score_samdif, ssim_score_samdif] = compute_similarity(ic_sout_BE1a, ic_sout_BE2, flag);
        % [nsim_score_samsam, ssim_score_samsam] = compute_similarity(ic_sout_BE1a, ic_sout_BE1b, flag);
        % [nsim_score_difsam, ssim_score_difsam] = compute_similarity(ic_sout_BE2, ic_sout_BE1b, flag);
        % 
        fprintf('Status:  %s\n', status);
        % fprintf('NSI for stimulus 1&2: %d\n', nsim_score_samdif);
        % fprintf('NSI for stimulus 1&1: %d\n', nsim_score_samsam);
        % fprintf('NSI for stimulus 2&1: %d\n', nsim_score_difsam);
    
        [nsim_score_deg, ssim_score_deg] = compute_similarity(ic_sout_BE1a, ic_sout_BE2,flag);
        nsimArray_degraded = [nsimArray_degraded; nsim_score_deg]; % Store the NSIM score in the array
        ssimArray_degraded = [ssimArray_degraded; ssim_score_deg]; % Store the SSIM score in the array
    
        manage_data(i+j, dB_signal, dependent_var, cihc_degraded, cohc_degraded, ...
            bm_degraded, condition, competitor_stimulus, nsim_score_deg, ssim_score_deg);
    
    end
    
    figure;
    plot(varArray, nsimArray_ref, 'r-', 'LineWidth', 2); % Plot the first line in red
    hold on; % Hold on to add the second line to the same plot
    plot(varArray, nsimArray_degraded, 'b-', 'LineWidth', 2); % Plot the second line in blue
    hold on; % Hold on to add the second line to the same plot
    plot(varArray, ssimArray_ref, 'g-', 'LineWidth', 2); % Plot the first line in green
    hold on; % Hold on to add the second line to the same plot
    plot(varArray, ssimArray_degraded, '-', 'LineWidth', 2); % Plot the second line in yellow
    hold off;
    
    % Add labels and title
    xlabel('Signal Frequency [Hz]');
    ylabel('Similarity Metric [0 - 1]');
    title(title_text);
    
    % Add legend
    legend('Reference NSIM', 'Degraded NSIM', 'Reference SSIM', 'Degraded SSIM');
    
    % Show grid
    grid on;

end

bm = 3;
cihc = 1;
cohc = 1;
title = 'Neurogram Similarity, numCF=1000, BM = 3';
run_analysis(bm, cihc, cohc, title)

bm = 1;
cihc = .3;
cohc = 1;
title = 'Neurogram Similarity, numCF=1000, cihc = .3';
run_analysis(bm, cihc, cohc, title)

bm = 1;
cihc = 1;
cohc = .3;
title = 'Neurogram Similarity, numCF=1000, cohc = .3';
run_analysis(bm, cihc, cohc, title)

bm = 3;
cihc = .3;
cohc = .3;
title = 'Neurogram Similarity, numCF=1000, full impairment bm = 3, cihc = .3, cohc = .3';
run_analysis(bm, cihc, cohc, title)