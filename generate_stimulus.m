%% function generates stimuli for IMAP test battery =================================================================================

% written by Sascha Muhlinghaus (05/08/2024)
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
            % Convert to dB SPL
            dB_noise = dB_noise + 10 * log10(800);
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
            dB_noise = dB_noise + 10 * log10(800);
            % Duration between noise and stimulus 
            pause_duration = varargin{2}/1000;
            % Generate noise signal (white noise)
            noise = sig_bandpassnoise(1000, fs_stim, 0.3, dB_noise, 800);
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
            % Generate noise signal (bandpass white noise)
            noise = sig_bandpassnoise(1000, fs_stim, 0.3, dB_noise, 800);
            % Ensure the tone fits within the noise duration
            if length(signal) > length(noise)
                error('The tone duration exceeds the noise duration.');
            end
   
            % Random start position for the tone within the noise interval
            max_start_index = length(noise) - length(signal) + 1;
            tone_start_index = randi([1, max_start_index]);
            % Initialize the stimulus with the noise
            stim = noise;
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
            dB_noise = dB_noise + 10 * log10(1200);
            % Generate noise signal (bandpass white noise)
            notched_noise = sig_notchednoise(1000, fs_stim, 1, dB_noise, .6, .2);
            % Ensure the tone fits within the noise duration
            if length(signal) > length(notched_noise)
                error('The tone duration exceeds the noise duration.');
            end

            % Random start position for the tone within the noise interval
            max_start_index = length(notched_noise) - length(signal) + 1;
            tone_start_index = randi([1, max_start_index]);
            % Initialize the stimulus with the filtered noise
            stim = notched_noise;
            plotfftreal(fftreal(notched_noise),fs_stim,100);
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