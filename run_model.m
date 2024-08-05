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

%% stimuli ================================================================

condition = 'backward masking';  % condition_flag: 'before', 'after', 'simultaneous', 'notch', 'frequency discrimination'
signal_freq = 1000; % tone_freq: frequency of the tone in Hz
signal_length_ms = 20; % signal_length_ms: length of the stimulus in (ms)
dB_signal = 75; % dB_signal: desired decibel level of the tone

db_noise = 30; % varargin1: inputs depending on the chosen condition 
gap = 50; % varargin2: inputs depending on the chosen condition
freq_limit = 1500; % varargin3: inputs depending on the chosen condition
increment = 10; % varargin4: inputs depending on the chosen condition

[wav1,Fs_wav1] = generate_stimulus(condition, signal_freq, signal_length_ms, dB_signal, db_noise, gap); % Generate stimulus

% [wav1,Fs_wav1] = audioread('sample.wav');

spl = 75;

Pref = 20*10^(-6);

if ~isnan(spl)
    wav1 = wav1*(Pref*10.^(spl/20)/rms(wav1));
end

stimulus = resample(wav1,Fs,Fs_wav1);

%% generate a neurogram for a specific hearing loss =======================

ag_dbloss = [0 0 0 0 0 0 0]; % modify loss here (+ db, linked to ag_fs)
ag_fs = [125 250 500 1e3 2e3 4e3 8e3]; % audiometric frequencies

%% define auditory pathway parameters (AN, CN, IC) ========================

species = 2; % 1 = cat; 2 = human
numCF = 50; % number of inner hair cells 
CF_range = [125 14e3]; % frequency range
cf = audspace(CF_range(1),CF_range(2),numCF); % generate cfs
bm = 1; % scale basilar membrane tuning, 1 = normal, > 1 = broader

fiberType = 1; % AN fiber type, 1 = L, 2 = M, 3 = H spontaneous rate, 4 = ratio defined below
fiber_num = 20;% number of fibres
numH = 10; % high spontaneous rate fibres
numM = 5; % medium spontaneous rate fibres
numL = 5; % low spontaneous rate fibres

cohc = 1; % OHC function scaling factor: 1 denotes normal function
cihc = 1; % IHC function scaling factor: 1 denotes normal function

BMF = 100; % default CN/IN best modulation frequency   

%% run models and save data =============================================================

filename = 'test_run.mat';

[r_mean, ihc, ic_sout_BE, cn_sout_contra] = simulate_auditory_response(stimulus, Fs, cf, bm, fiberType, numH, numM, numL, fiber_num, cohc, cihc, BMF, filename);

%% plot output ============================================================

plot_neurograms(stimulus, r_mean, ihc, ic_sout_BE, cn_sout_contra, Fs, cf);



