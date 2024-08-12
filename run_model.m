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

Fs = 100e3; % "model sampling frequency" [samples/sec]

%% stimulus parameter =====================================================

condition = 'vcv';  % condition_flag: 'before', 'after', 'simultaneous', 'notch', 'frequency discrimination'
signal_freq = 1000; % tone_freq: frequency of the tone in Hz
signal_length_ms = 20; % signal_length_ms: length of the stimulus in (ms)
dB_signal = 70; % dB_signal: desired decibel level of the tone

db_noise = 30; % varargin1: inputs depending on the chosen condition 
gap = 50; % varargin2: inputs depending on the chosen condition
freq_limit = 1500; % varargin3: inputs depending on the chosen condition
increment = 10; % varargin4: inputs depending on the chosen condition
vcv_string = 'a b a';
debug_flag = false;

%% generate a neurogram for a specific hearing loss =======================

ag_dbloss = [0 0 0 0 0 0 0]; % modify loss here (+ db, linked to ag_fs)
ag_fs = [125 250 500 1e3 2e3 4e3 8e3]; % audiometric frequencies

%% define auditory pathway parameters (AN, CN, IC) ========================

species = 2; % 1 = cat; 2 = human
numCF = 50; % number of inner hair cells 
CF_range = [125 14e3]; % frequency range
cf = audspace(CF_range(1),CF_range(2),numCF); % generate cfs

% scale basilar membrane tuning, 1 = normal, > 1 = broader
bm_ref = 1; 
bm_degraded = 3;
% OHC function scaling factor: 1 denotes normal function
cohc_ref = 1; 
cohc_degraded = 1; 
% IHC function scaling factor: 1 denotes normal function
cihc_ref = 1;
cihc_degraded = 1;


fiberType = 1; % AN fiber type, 1 = L, 2 = M, 3 = H spontaneous rate, 4 = ratio defined below
fiber_num = 20;% number of fibres
numH = 10; % high spontaneous rate fibres
numM = 5; % medium spontaneous rate fibres
numL = 5; % low spontaneous rate fibres
BMF = 100; % default CN/IN best modulation frequency   

%% run models and save data =============================================================

[stimulus,Fs_stimulus, n_samples] = generate_stimulus(condition, signal_freq, signal_length_ms, dB_signal, debug_flag, vcv_string); % Generate stimulus

nsimArray = [];

for i = 1:n_samples
    
    filename = 'test_run.mat';
    
    [r_mean_ref, ihc_ref, ic_sout_BE_ref, cn_sout_contra_ref] = simulate_auditory_response(stimulus, Fs, cf, bm_ref, fiberType, numH, numM, numL, fiber_num, cohc_ref, cihc_ref, BMF, filename);
    
    [r_mean, ihc, ic_sout_BE, cn_sout_contra] = simulate_auditory_response(stimulus, Fs, cf, bm_degraded, fiberType, numH, numM, numL, fiber_num, cohc_degraded, cihc_degraded, BMF, filename);
    % signal_freq = signal_freq + 10;
    % 
    % if i ~= n_samples
    %     signal_frequencies = [signal_frequencies; signal_freq];
    % end
    % fprintf('Signal Frequency: %s\n', mat2str(signal_frequencies)); % Print the entire array

    % plot output ============================================================

    % Plot the reference neurograms
    plot_neurograms(stimulus, r_mean_ref, ihc_ref, ic_sout_BE_ref, cn_sout_contra_ref, Fs, cf);
    
    % Plot the degraded neurograms
    plot_neurograms(stimulus, r_mean, ihc, ic_sout_BE, cn_sout_contra, Fs, cf);
    
    nsim_score = compute_nsim(ic_sout_BE_ref, ic_sout_BE);
    %fprintf('NSIM score: %f\n', nsim_score);
    nsimArray = [nsimArray; nsim_score]; % Store the NSIM score in the array
    fprintf('NSIM Array: %s\n', mat2str(nsimArray)); % Print the entire array
end

% % Create the plot
% figure;
% plot(signal_frequencies, nsimArray, 'r-', 'LineWidth', 2); % Plot the first line in green
% hold on; % Hold on to add the second line to the same plot
% plot(signal_frequencies, nsimArray2, 'b-', 'LineWidth', 2); % Plot the second line in blue
% hold off;
% 
% % Add labels and title
% xlabel('Signal Frequency');
% ylabel('NSIM score');
% title('NSIM Scores, numCF=50');
% 
% % Add legend
% legend('BM=1', 'BM=3');
% 
% % Show grid
% grid on;




