%% Function generates stimuli for IMAP test battery =======================

% written by Sascha Muhlinghaus (05/08/2024)
% contact: saschamuhlinghaus@gmail.com

function [stimulus, stimulus_fs, N_samples, varargout] = generate_stimulus(condition_flag, tone_freq, signal_length_ms, dB_signal, debug_flag, varargin)
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
    stimulus_fs = 100e3; 
    
    % Convert ms to seconds
    stim_length = signal_length_ms / 1000; 
    
    % Generate time vector for the signal length
    t_signal = (0:round(stim_length * stimulus_fs) - 1) / stimulus_fs;
    
    % Generate tone signal (sine wave)
    signal = scaletodbspl(sin(2*pi*tone_freq*t_signal),dB_signal);

    % Define the ramp duration in seconds
    rampDuration = 0.01; % 10 ms ramp duration
    
    % Convert the ramp duration to samples
    rampSamples = round(rampDuration * stimulus_fs);
    
    % Generate cosine ramp
    t_ramp = (0:rampSamples-1)' / rampSamples;
    ramp = 0.5 * (1 - cos(pi * t_ramp));

    % Apply the ramp to the beginning of the signal
    signal(1:rampSamples) = signal(1:rampSamples) .* ramp';

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
            pause_duration_samples = round(pause_duration * stimulus_fs);
            total_length_samples = length(noise) + pause_duration_samples + length(signal);
            % Initialize stimulus variable
            stimulus = zeros(1, total_length_samples);
            % Place noise before the tone with a pause in between
            tone_start = length(noise) + pause_duration_samples + 1;
            stimulus(1:length(noise)) = noise;
            stimulus(tone_start:tone_start + length(signal) - 1) = signal;

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
            pause_duration_samples = round(pause_duration * stimulus_fs);
            total_length_samples = length(signal) + pause_duration_samples + length(noise);
            stimulus = zeros(1, total_length_samples);
            % Place noise after the tone with a pause in between
            noise_start = length(signal) + pause_duration_samples + 1;
            stimulus(1:length(signal)) = signal;
            stimulus(noise_start:noise_start + length(noise) - 1) = noise;

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
            if length(signal) > length(noise)
                error('The tone duration exceeds the noise duration.');
            end
   
            % Random start position for the tone within the noise interval
            max_start_index = length(noise) - length(signal) + 1;
            tone_start_index = randi([1, max_start_index]);
            % Initialize the stimulus with the noise
            stimulus = noise;
            % Place the tone within the noise interval
            stimulus(tone_start_index:tone_start_index + length(signal) - 1) = signal;
            stimulus = stimulus';

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
            if length(signal) > length(noise)
                error('The tone duration exceeds the noise duration.');
            end

            % Generate band-stop filter to remove the tone frequency from the noise
            [b, a] = butter(4, [(tone_freq - 200) / (stimulus_fs / 2), (tone_freq + 200) / (stimulus_fs / 2)], 'stop');
            noise_filtered = filter(b, a, noise);
            % Random start position for the tone within the noise interval
            max_start_index = length(noise_filtered) - length(signal) + 1;
            tone_start_index = randi([1, max_start_index]);
            % Initialize the stimulus with the filtered noise
            stimulus = noise_filtered;
            % Place the tone within the filtered noise interval
            stimulus(tone_start_index:tone_start_index + length(signal) - 1) = signal;
            stimulus = stimulus';
    
        case 'frequency discrimination'
            % If flag is 'frequency discrimination', expect 2 additional arguments
            if nArgs < 2
                error(['The condition ''frequency discrimnation'' requires 2 additional arguments ' ...
                    ', the end of the ''frequency spectrum'' (Hz) and the ''increment'' (Hz).']);
            end

            end_frequency = varargin{1};
            increment = varargin{2};

            % Number of tones
            n_samples = ((end_frequency - tone_freq) / increment);
            % Initialize an empty array to store the tones
            tonesArray = zeros(n_samples, length(t_signal));  % Preallocate for efficiency
            % Define frequency for the first tone
            frequency = tone_freq;
        
            for i = 1:n_samples
                % Generate tone
                tone = scaletodbspl(sin(2*pi*frequency*t_signal), dB_signal);       
                % Store the generated tone in the array
                tonesArray(i, :) = tone;
                % Update frequency for the next tone
                frequency = frequency + increment;
            end
            
            % Return the tones
            stimulus = tonesArray;

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
    stimulus = resample(audioData, stimulus_fs, originalFs);   
    % Adjust the amplitude 
    stimulus = scaletodbspl(stimulus, dB_signal);
    % Trim signal to 800ms max length
    stimulus = stimulus(1:stimulus_fs-(stimulus_fs*.2))';
        otherwise
            error('Invalid condition flag. Choose ''backward masking'', ''forward masking'', ''no notch'', ''notch'', ''frequency discrimination'' or ''vcv''.');
    end

    % Return number of stimuli
    [N_samples, cols] = size(stimulus);

    % If debugging flag is true, return the size of the stimulus
    if debug_flag
        varargout{1} = size(stimulus);
    end
end


%% Test function 
% Change condition for generation of different test stimuli
condition = 'vcv';  % condition_flag: 'backward masking', 'forward masking', 'no notch', 'notch', 'frequency discrimination', 'vcv'
% Input parameters
signal_freq = 1000; % tone_freq: frequency of the tone in Hz
signal_length_ms = 20; % signal_length_ms: length of the stimulus in (ms)
dB_signal = 80; % dB_signal: desired decibel level of the tone
% Condition dependent input parameters
db_noise = 30; % varargin1: inputs depending on the chosen condition 
gap = 50; % varargin2: inputs depending on the chosen condition
freq_limit = 1500; % varargin3: inputs depending on the chosen condition
increment = 10; % varargin4: inputs depending on the chosen condition
vcv_string = 'a b a';
debug_flag = false;

%% Use for condition 'backward masking', 'forward masking', 'no notch', 'notch'
%[stimulus, fs, N_samples] = generate_stimulus(condition, signal_freq, signal_length_ms, dB_signal, debug_flag, db_noise, gap); % Generate stimulus
%% Use for condition 'vcv'
[stimulus, fs, N_samples] = generate_stimulus(condition, signal_freq, signal_length_ms, dB_signal, debug_flag, vcv_string); % Generate stimulus
%% Use for condition 'frequency discrimination'
%[stimulus, fs, N_samples] = generate_stimulus(condition, signal_freq, signal_length_ms, dB_signal, debug_flag, freq_limit, increment);
for i = 1:length(stimulus(:,1))
    sound(stimulus(i,:),fs);
end
fprintf('number of stimuli: %d\n', N_samples);


