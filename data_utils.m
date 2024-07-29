function tones = createTone(Start_frequency, End_frequency, increment, dB, tone_duration)
    % createTone generates tones with specified start frequency, end frequency, increment, 
    % dB level, and duration.
    %
    % Returns the generated tones.
    %
    % Arguments:
    %   Start_frequency (Hz): Start frequency of the interval.
    %   End_freqeuncy (Hz): End frequency of the interval.
    %   dB (decibels): Sound pressure level of the tone.
    %   tone_duration (seconds): Duration of the tone.
    %
    % Returns:
    %   tonesArray (array): Generated tone signals.

    % Sampling rate
    fs = 44100; % Standard sampling rate for audio

    % Time vector
    t = 0:1/fs:tone_duration;

    % Number of tones
    n_samples = floor((End_frequency - Start_frequency) / increment) + 1;

    % Initialize an empty array to store the tones
    tonesArray = [];

    % Define frequency for the first tone
    frequency = Start_frequency;

    for i = 1:n_samples
        % Generate tone
        tone = sin(2 * pi * frequency * t);

        % Convert dB to amplitude (assuming full scale = 0 dBFS)
        amplitude = 10^(dB / 20);

        % Adjust tone amplitude
        tone = amplitude * tone;

        % Append the generated tone to the tonesArray
        tonesArray = [tonesArray; tone];

        % Update frequency for the next tone
        frequency = frequency + increment;
    end
    
    % return the tones
    tones = tonesArray;
end

% Example usage:

start_fs = 440;
end_fs = 500;
stepper = 5;
db = -10;
dur = .5;

tonesSamples = createTone(start_fs, end_fs, stepper, db, dur); % Generates tones from 440 Hz to 500 Hz with 10 Hz increment at -10 dB for 0.5 seconds

for i = 1:size(tonesSamples, 1)
    sound(tonesSamples(i, :), 44100);
    pause(0.5 + 0.5); % Pause to allow the tone to finish playing
end
