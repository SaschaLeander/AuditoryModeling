%% Function generates stimuli for IMAP test battery =======================

% written by Sascha Muhlinghaus 
% contact: saschamuhlinghaus@gmail.com

function [stim, fs_stim, length_stim, varargout] = generate_stimulus(condition_flag, tone_freq, signal_length_ms, dB_signal, debug_flag, varargin)
    % Input Parameters:
    % condition_flag: 'forward masking', 'backward masking', 'no notch',
    % 'notch', 'frequency discrimination', 'vcv'
    % tone_freq: frequency of the tone in Hz
    % signal_length_ms: length of the stimulus in (ms)
    % dB_signal: desired decibel level of the tone
    % debug_flag: if true, additional debugging information will be returned
    % varargin: inputs depending on the chosen condition 

    % Check the number of input arguments
    nArgs = length(varargin);
    
    % Sampling rate
    fs_stim = 100e3; 
    
    % Convert ms to seconds
    stim_length = signal_length_ms / 1000; 
    
    % Generate time vector for the signal length
    t_signal = (0:round(stim_length * fs_stim) - 1) / fs_stim;
    
    % Generate tone signal (sine wave)
    signal = scaletodbspl(sin(2*pi*tone_freq*t_signal),dB_signal);

    % Define the ramp duration in seconds
    rampDuration = 0.01; % 10 ms ramp duration
    
    % Convert the ramp duration to samples
    rampSamples = round(rampDuration * fs_stim);
    
    % Generate cosine ramp
    t_ramp = (0:rampSamples-1)' / rampSamples;
    ramp = 0.5 * (1 - cos(pi * t_ramp));

    % Apply the ramp to the beginning of the signal (ramp up)
    signal(1:rampSamples) = signal(1:rampSamples) .* ramp';

    % Apply the ramp to the end of the signal (ramp down)
    % signal(end-rampSamples+1:end) = signal(end-rampSamples+1:end) .* flipud(ramp');

    % Apply the condition
    switch condition_flag
        case 'forward masking'
            % If flag is 'forward masking', expect 2 additional arguments
            if nArgs < 2
                error(['The condition ''forward masking'' requires 2 additional arguments ' ...
                    'indicating the duration of the pause between noise' ...
                    'and the dB level of the noise.']);
            end

            % dB level of noise
            dB_noise = varargin{1};
            % Convert to dB SPL
            dB_noise = dB_noise + 10 * log10(800);
            fprintf('noise in dB SPL: %d\n', dB_noise);
            % Duration between noise and stimulus 
            pause_duration = varargin{2}/1000;
            % Generate noise signal (white noise)
            noise = sig_bandpassnoise(1000, fs_stim, 0.3, dB_noise, 800);
            % Initialize the full stimulus length
            pause_duration_samples = round(pause_duration * fs_stim);
            total_length_samples = length(noise) + pause_duration_samples + length(signal);
            % Initialize stimulus variable
            stim = zeros(1, total_length_samples);
            % Place noise before the tone with a pause in between
            tone_start = length(noise) + pause_duration_samples + 1;
            stim(1:length(noise)) = noise;
            stim(tone_start:tone_start + length(signal) - 1) = signal;

        case 'backward masking'
            % If flag is 'backward masking', expect 2 additional arguments
            if nArgs < 2
                error(['The condition ''backward masking'' requires 2 additional arguments ' ...
                    'indicating the duration of the pause between noise' ...
                    'and the dB level of the noise.']);
            end

            % dB level of noise
            dB_noise = varargin{1};
            % Convert to dB SPL
            dB_noise = dB_noise + (10 * log10(800));
            fprintf('noise in dB SPL: %d\n', dB_noise);
            % Duration between noise and stimulus 
            pause_duration = varargin{2}/1000;
            % Generate noise signal (white noise)
            noise = sig_bandpassnoise(1000, fs_stim, 0.3, 30, 800);
            % Initialize the full stimulus length
            pause_duration_samples = round(pause_duration * fs_stim);
            total_length_samples = length(signal) + pause_duration_samples + length(noise);
            stim = zeros(1, total_length_samples);
            % Place noise after the tone with a pause in between
            noise_start = length(signal) + pause_duration_samples + 1;
            stim(1:length(signal)) = signal;
            stim(noise_start:noise_start + length(noise) - 1) = noise;

        case 'no notch'
            % If flag is 'no notch', expect 1 additional argument
            if nArgs < 1
                error(['The condition ''no notch'' requires 1 additional argument ' ...
                    'indicating the dB level of the noise.']);
            end
    
            % dB level of noise
            dB_noise = varargin{1};
            % Convert to dB SPL
            dB_noise = dB_noise + 10 * log10(800);
            fprintf('noise in dB SPL: %d\n', dB_noise);
            % Generate noise signal (bandpass white noise)
            noise = sig_bandpassnoise(1000, fs_stim, 0.3, dB_noise, 800);
            % Ensure the tone fits within the noise duration
            if length(signal) > length(noise)
                error('The tone duration exceeds the noise duration.');
            end
            
            % Initialize the stimulus with the noise
            stim = noise;
            tone_start_index = length(stim)/2;
            % Place the tone within the noise interval
            stim(tone_start_index:tone_start_index + length(signal) - 1) = signal;
            stim = stim';

        case 'notch'
            % If flag is 'notch', expect 1 additional argument
            if nArgs < 1
                error(['The condition ''notch'' requires 1 additional argument ' ...
                    'indicating the dB level of the noise.']);
            end

            % dB level of noise
            dB_noise = varargin{1};
            % Convert to dB SPL
            dB_noise = dB_noise + (10 * log10(1200));
            fprintf('noise in dB SPL: %d\n', dB_noise);
            % Generate noise signal (bandpass white noise)
            notched_noise = sig_notchednoise(signal_freq, fs_stim, .3, dB_noise, 1.2, 0.2);
            % Ensure the tone fits within the noise duration
            if length(signal) > length(notched_noise)
                error('The tone duration exceeds the noise duration.');
            end

            % Initialize the stimulus with the filtered noise
            stim = notched_noise;
            tone_start_index = length(stim)/2;
            % Place the tone within the filtered noise interval
            stim(tone_start_index:tone_start_index + length(signal) - 1) = signal;
            stim = stim';
    
        case 'frequency discrimination'
            
            % Return tone
            stim = signal;  

        case 'vcv'
           
            % Check if the input is a string
            if ~ischar(varargin{1}) && ~isstring(varargin{1})
                error('Input must be a string or a character array.');
            end
            
            % Convert input to character array if it is a string
            inputText = char(varargin{1}); 

            % Create a temporary file path for saving the speech
            tempWavFile = tempname + ".wav";    
            % Create a COM server for SAPI.SpVoice
            speaker = actxserver('SAPI.SpVoice');    
            % Create a Speech Stream File for output
            fileStream = actxserver('SAPI.SpFileStream');
            % Set the output format and file location
            fileStream.Format.Type = 'SAFT16kHz16BitMono'; % Default speech generation rate
            invoke(fileStream, 'Open', tempWavFile, 3); % Open for writing
            speaker.AudioOutputStream = fileStream;    
            % Speak the text and save to the file
            invoke(speaker, 'Speak', inputText);    
            % Clean up
            invoke(fileStream, 'Close');
            delete(speaker);
            delete(fileStream);    
            % Read the saved WAV file
            [audioData, originalFs] = audioread(tempWavFile);  
            % Resample to the desired sampling rate (fs)
            stim = resample(audioData, fs_stim, originalFs);   
            % Adjust the amplitude 
            stim = scaletodbspl(stim, dB_signal);
            % Trim signal to 800ms max length
            stim = stim(1:fs_stim-(fs_stim*.2))';
        otherwise
            error('Invalid condition flag. Choose ''backward masking'', ''forward masking'', ''no notch'', ''notch'', ''frequency discrimination'' or ''vcv''.');
    end

    % Return length of stimulus
    length_stim = length(stim) / fs_stim;

    % If debugging flag is true, return the size of the stimulus
    if debug_flag
        varargout{1} = size(stim);
    end
end

%% Test function 
% Change condition for generation of different test stimuli
condition = 'vcv';  % condition_flag: 'backward masking', 'forward masking', 'no notch', 'notch', 'frequency discrimination', 'vcv'
% Input parameters
signal_freq = 1000; % tone_freq: frequency of the tone in Hz
signal_length_ms = 20; % signal_length_ms: length of the stimulus in (ms)
dB_signal = 75; % dB_signal: desired decibel level of the tone
% Condition dependent input parameters
db_noise = 30; % varargin1: inputs depending on the chosen condition 
gap = 50; % varargin2: inputs depending on the chosen condition
vcv_string = 'a g a';
debug_flag = false;

%% Use for condition 'frequency discrimination', 'backward masking', 'forward masking', 'no notch', 'notch'
% [stimulus, fs, length] = generate_stimulus(condition, signal_freq, signal_length_ms, dB_signal, debug_flag, db_noise, gap); % Generate stimulus
%% Use for condition 'vcv'
[stimulus, fs, length] = generate_stimulus(condition, signal_freq, signal_length_ms, dB_signal, debug_flag, vcv_string); % Generate stimulus

sound(stimulus, fs)
fprintf('stimulus length in (ms): %d\n', length);

%% define auditory pathway parameters (AN, CN, IC) ========================

species = 2; % 1 = cat; 2 = human
numCF = 10; % number of inner hair cells
CF_range = [125 14e3]; % frequency range
cf = audspace(CF_range(1),CF_range(2),numCF); % generate cfs
bm = 1; % scale basilar membrane tuning, 1 = normal, > 1 = broader
cohc = 1; % OHC function scaling factor: 1 denotes normal function
cihc = 1; % IHC function scaling factor: 1 denotes normal function
fiberType = 1; % AN fiber type, 1 = L, 2 = M, 3 = H spontaneous rate, 4 = ratio defined below
fiber_num = 20;% number of fibres
numH = 10; % high spontaneous rate fibres
numM = 5; % medium spontaneous rate fibres
numL = 5; % low spontaneous rate fibres
BMF = 100; % default CN/IN best modulation frequency


% Model auditory
[r_mean, ihc, ic_sout_BE, cn_sout_contra] = simulate_auditory_response(stimulus, ...
    fs, cf, bm, fiberType, numH, numM, numL, fiber_num, cohc, cihc, BMF);


% Plot neurograms
plot_neurograms(stimulus, r_mean, ihc, ic_sout_BE, cn_sout_contra, fs, cf);
