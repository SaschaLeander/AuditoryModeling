%% function generates stimuli for IMAP test battery =================================================================================

% written by Sascha Muhlinghaus (05/08/2024)
% contact: saschamuhlinghaus@gmail.com

function [stimulus, stimulus_fs] = generate_stimulus(condition_flag, tone_freq, signal_length_ms, dB_signal, varargin)
 % Input Parameters:
    % condition_flag: 'forward masking', 'backward masking', 'no notch',
    % 'notch', 'frequency discrimination', 'vcv'
    % tone_freq: frequency of the tone in Hz
    % signal_length_ms: length of the stimulus in (ms)
    % dB_signal: desired decibel level of the tone
    % varargin: inputs depending on the chosen condition 

    % Check the number of input arguments
    nArgs = length(varargin);
    
    % Sampling rate
    stimulus_fs = 44100; % 44.1 kHz standard sampling rate for audio

    stim_length = signal_length_ms / 1000; % Convert ms to seconds

    % Convert dB level to linear amplitude
    amp_signal = 10^(dB_signal / 20);
    
    % Generate time vector for the signal length
    t_signal = (0:round(stim_length * stimulus_fs) - 1) / stimulus_fs;
    
    % Generate tone signal (sine wave)
    signal = amp_signal * sin(2 * pi * tone_freq * t_signal);

    % Define the ramp duration in seconds
    rampDuration = 0.01; % For example, 10 ms
    
    % Convert the ramp duration to samples
    rampSamples = round(rampDuration * stimulus_fs);
    
    % Generate the cosine ramp
    t_ramp = (0:rampSamples-1)' / rampSamples;
    ramp = 0.5 * (1 - cos(pi * t_ramp));

    % Create a signal with ramp added before and after the tone
    total_length_samples = length(signal) + 2 * rampSamples;
    signal_with_ramp = zeros(1, total_length_samples);
    signal_with_ramp(rampSamples+1:rampSamples+length(signal)) = signal;

    % Apply the ramp to the beginning and end of the signal
    signal_with_ramp(1:rampSamples) = signal_with_ramp(1:rampSamples) .* ramp';
    signal_with_ramp(end-rampSamples+1:end) = signal_with_ramp(end-rampSamples+1:end) .* flip(ramp');

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
            % Duration between noise and stimulus 
            pause_duration = varargin{2}/1000;
            % Generate noise signal (white noise)
            noise = sig_bandpassnoise(1000, stimulus_fs, 0.3, dB_noise, 800);
            % Initialize the full stimulus length
            pause_duration_samples = pause_duration * stimulus_fs;
            total_length_samples = length(noise) + round(pause_duration_samples) + length(signal_with_ramp);

            stimulus = zeros(1, total_length_samples);
            % Place noise before the tone with a pause in between
            tone_start = length(noise) + pause_duration_samples + 1;
            stimulus(1:length(noise)) = noise;
            stimulus(tone_start:tone_start + length(signal_with_ramp) - 1) = signal_with_ramp;

            fprintf('Stimulus length in (s): %d\n', length(stimulus)/stimulus_fs);
          
        case 'backward masking'
            % If flag is 'backward masking', expect 2 additional arguments
            if nArgs < 2
                error(['The condition ''backward masking'' requires 2 additional arguments ' ...
                    'indicating the duration of the pause between noise' ...
                    'and the dB level of the noise.']);
            end

            % dB level of noise
            dB_noise = varargin{1};
            % Duration between noise and stimulus 
            pause_duration = varargin{2}/1000;
            % Generate noise signal (white noise)
            noise = sig_bandpassnoise(1000, stimulus_fs, 0.3, dB_noise, 800);
            % Initialize the full stimulus length
            pause_duration_samples = pause_duration * stimulus_fs;
            total_length_samples = length(signal_with_ramp) + round(pause_duration_samples) + length(noise);
            stimulus = zeros(1, total_length_samples);
            % Place noise after the tone with a pause in between
            noise_start = length(signal_with_ramp) + pause_duration_samples + 1;
            stimulus(1:length(signal_with_ramp)) = signal_with_ramp;
            stimulus(noise_start:noise_start + length(noise) - 1) = noise;

            fprintf('Stimulus length in (s): %d\n', length(stimulus)/stimulus_fs);
            
        case 'no notch'
            % If flag is 'no notch', expect 1 additional argument
            if nArgs < 1
                error(['The condition ''no notch'' requires 1 additional argument ' ...
                    'indicating the dB level of the noise.']);
            end
    
            % dB level of noise
            dB_noise = varargin{1};
            % Generate noise signal (bandpass white noise)
            noise = sig_bandpassnoise(1000, stimulus_fs, 0.3, dB_noise, 800);
            % Ensure the tone fits within the noise duration
            if length(signal_with_ramp) > length(noise)
                error('The tone duration exceeds the noise duration.');
            end
   
            % Random start position for the tone within the noise interval
            max_start_index = length(noise) - length(signal_with_ramp) + 1;
            tone_start_index = randi([1, max_start_index]);
            % Initialize the stimulus with the noise
            stimulus = noise;
            % Place the tone within the noise interval
            stimulus(tone_start_index:tone_start_index + length(signal_with_ramp) - 1) = signal_with_ramp;
            
            fprintf('Stimulus length in (s): %d\n', length(stimulus)/stimulus_fs);

        case 'notch'
            % If flag is 'notch', expect 1 additional argument
            if nArgs < 1
                error(['The condition ''notch'' requires 1 additional argument ' ...
                    'indicating the dB level of the noise.']);
            end

            % dB level of noise
            dB_noise = varargin{1};
            % Generate noise signal (bandpass white noise)
            noise = sig_bandpassnoise(1000, stimulus_fs, 0.3, dB_noise, 1200);
            % Ensure the tone fits within the noise duration
            if length(signal_with_ramp) > length(noise)
                error('The tone duration exceeds the noise duration.');
            end

            % Generate band-stop filter to remove the tone frequency from the noise
            [b, a] = butter(4, [(tone_freq - 200) / (stimulus_fs / 2), (tone_freq + 200) / (stimulus_fs / 2)], 'stop');
            noise_filtered = filter(b, a, noise);
            % Random start position for the tone within the noise interval
            max_start_index = length(noise_filtered) - length(signal_with_ramp) + 1;
            tone_start_index = randi([1, max_start_index]);
            % Initialize the stimulus with the filtered noise
            stimulus = noise_filtered;
            % Place the tone within the filtered noise interval
            stimulus(tone_start_index:tone_start_index + length(signal_with_ramp) - 1) = signal_with_ramp;

            fprintf('Stimulus length in (s): %d\n', length(stimulus)/stimulus_fs);
    
        case 'frequency discrimination'
            % If flag is 'frequency discrimination', expect 2 additional arguments
            if nArgs < 2
                error(['The condition ''frequency discrimnation'' requires 2 additional arguments ' ...
                    ', the end of the ''frequency spectrum'' (Hz) and the ''increment'' (Hz).']);
            end

            end_frequency = varargin{1};
            increment = varargin{2};

            % Number of tones
            n_samples = floor((end_frequency - tone_freq) / increment) + 1;
            % Initialize an empty array to store the tones
            tonesArray = [];
            % Define frequency for the first tone
            frequency = tone_freq;
        
            for i = 1:n_samples
                % Generate tone
                tone = sin(2 * pi * frequency * t_signal);
                % Convert dB to amplitude (assuming full scale = 0 dBFS)
                amplitude = 10^(dB_signal / 20);
                % Adjust tone amplitude
                tone = amplitude * tone;
                % Append the generated tone to the tonesArray
                tonesArray = [tonesArray; tone];
                % Update frequency for the next tone
                frequency = frequency + increment;
            end
            
            % return the tones
            stimulus = tonesArray;
            fprintf('Stimulus length in (s): %d\n', length(stimulus(1,:))/stimulus_fs);
        
        case 'vcv'
            n = 5;
            vcvArray = {};
        
            for i = 0:n-1
                % Define location of vcv stimulus
                file_name = strcat('vcv_samples/vcv_stimulus', int2str(i), '.wav');
                % Read vcv file
                [wav, Fs_wav] = audioread(file_name);
                % Concatenate audio and sampling frequency into cell array
                vcvArray = [vcvArray; {wav, Fs_wav}];
            end
            stimulus = vcvArray;

        otherwise
            error('Invalid condition flag. Choose ''backward masking'', ''forward masking'', ''no notch'', ''notch'', ''frequency discrimination'' or ''vcv''.');
    end
end
