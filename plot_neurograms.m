%% Function plots auditory response =======================

% written by Sam Jones (30/07/2024)

function plot_neurograms(stimulus, r_mean, ihc, ic_sout_BE, cn_sout_contra, Fs, cf)

    % re-size matrices, accounting for zero-padding within the model
    [sizeN, sizeM] = size(r_mean);
    ihc = ihc(1:sizeN, 1:sizeM);
    ic_sout_BE = ic_sout_BE(1:sizeN, 1:sizeM);
    cn_sout_contra = cn_sout_contra(1:sizeN, 1:sizeM);
    
    % define common time vector
    t = (1:length(r_mean))*1/Fs;

    % plot figures
    figure 
    
    % raw stimulus
    subplot(4, 3, 1)
    t_stimulus = 1/Fs:1/Fs:length(r_mean)/Fs;
    minLength = min(length(t_stimulus), length(stimulus));
    t_stimulus = t_stimulus(1:minLength);
    stimulus = stimulus(1:minLength);
    plot(t_stimulus, stimulus)
    title('Stimulus waveform')
    xlabel('Time (s)')
    ylabel('Amplitude (Pa)')
    grid on

    % frequency spectrum subplot
    subplot(4, 3, 2)

    % calculate FFT
    n = length(stimulus);
    frequencies = (0:n-1)*(Fs/n);
    magnitude = abs(fft(stimulus))/n;
    
    % convert to dB SPL (assuming the reference pressure p0 = 20 µPa)
    referencePressure = 20e-6;
    magnitude_db_spl = 20 * log10(magnitude / referencePressure);
    
    % plot the magnitude in dB SPL by frequency
    semilogx(frequencies(1:floor(n/2)), magnitude_db_spl(1:floor(n/2)));
    xlim([125 8000]); 
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB SPL)')
    title('Stimulus frequency spectrum')
    set(gca, 'XTick', [20 50 100 200 500 1000 2000 5000 8000]);
    set(gca, 'XTickLabel', {'20', '50', '100', '200', '500', '1k', '2k', '5k', '8k'});
    grid on;

    % inner hair cell response
    subplot(4, 3, 3)
    plot(t, ihc);
    xlim([0 (length(stimulus)/Fs)*1.2]) %  scale x-axis by 1.2 * stimulus length
    title('Inner hair cell response')
    xlabel('Time (s)')
    ylabel('Volts')
    grid on

    % auditory nerve response 1
    subplot(4, 3, 4)
    imagesc(t, log10(cf/1e3), r_mean');
    axis xy;
    yticks = [0.125 0.5 2 8 16];
    set(gca, 'ytick', log10(yticks));
    set(gca, 'yticklabel', yticks);
    ylabel('CF (kHz)')
    xlabel('Time (s)')
    xlim([0 length(stimulus)/Fs])
    hcb = colorbar;
    set(get(hcb, 'ylabel'), 'string', 'Spikes / second');
    title('Auditory nerve response')
    grid on
    
    % auditory nerve response 2
    subplot(4, 3, 5)
    plot(t, r_mean)
    title('Auditory nerve response')
    xlim([0 length(stimulus)/Fs])
    xlabel('Time (s)')
    ylabel('Instantaneous mean spiking rate (sp/s)')
    grid on

    % auditory nerve mean rate by CF
    r_clipped = r_mean(1:length(stimulus),:); 
    r_mean_clipped = mean(r_clipped, 1);
    subplot(4, 3, 6);
    semilogx(cf, r_mean_clipped, 'LineWidth', 2);
    xlabel('Frequency (Hz)')
    ylabel('Mean rate (sp/s)')
    title('Auditory nerve response (mean rate)')
    xlim([0 8000]); 
    set(gca, 'XTick', [20 50 100 200 500 1000 2000 5000 8000]);
    set(gca, 'XTickLabel', {'20', '50', '100', '200', '500', '1k', '2k', '5k', '8k'});
    grid on;

    % cochlear nucleus response 1
    subplot(4, 3, 7)
    imagesc(t, log10(cf/1e3), cn_sout_contra');
    axis xy;
    yticks = [0.125 0.5 2 8 16];
    set(gca, 'ytick', log10(yticks));
    set(gca, 'yticklabel', yticks);
    ylabel('CF (kHz)')
    xlabel('Time (s)')
    xlim([0 length(stimulus)/Fs])
    hcb = colorbar;
    set(get(hcb, 'ylabel'), 'string', 'Spikes / second');
    title('Cochlear nucleus response')
    grid on

    % cochlear nucleus response 2
    subplot(4, 3, 8)
    plot(t, cn_sout_contra)
    title('Cochlear nucleus response')
    xlim([0 length(stimulus)/Fs])
    xlabel('Time (s)')
    ylabel('Instantaneous mean spiking rate (sp/s)')
    grid on

    % cochlear nucleus response mean rate by CF
    cn_clipped = cn_sout_contra(1:length(stimulus),:);
    cn_mean_clipped = mean(cn_clipped, 1);
    subplot(4, 3, 9);
    semilogx(cf, cn_mean_clipped, 'LineWidth', 2); 
    xlabel('Frequency (Hz)')
    ylabel('Mean rate (sp/s)')
    title('Cochlear nucleus response (mean rate)')
    xlim([0 8000]); 
    set(gca, 'XTick', [20 50 100 200 500 1000 2000 5000 8000]);
    set(gca, 'XTickLabel', {'20', '50', '100', '200', '500', '1k', '2k', '5k', '8k'});
    grid on;

    % inferior colliculus response 1
    subplot(4, 3, 10)
    imagesc(t, log10(cf/1e3), ic_sout_BE');
    axis xy;
    yticks = [0.125 0.5 2 8 16];
    set(gca, 'ytick', log10(yticks));
    set(gca, 'yticklabel', yticks);
    ylabel('CF (kHz)')
    xlabel('Time (s)')
    xlim([0 length(stimulus)/Fs])
    hcb = colorbar;
    set(get(hcb, 'ylabel'), 'string', 'Spikes / second');
    title('Inferior colliculus response')
    grid on

    % inferior colliculus response 2
    subplot(4, 3, 11)
    plot(t, ic_sout_BE)
    xlim([0 length(stimulus)/Fs])
    xlabel('Time (s)')
    ylabel('Instantaneous mean spiking rate (sp/s)')
    title('Inferior colliculus response')
    grid on

    % inferior colliculus response mean rate by CF
    ic_clipped = ic_sout_BE(1:length(stimulus),:);
    ic_mean_clipped = mean(ic_clipped, 1); 
    subplot(4, 3, 12);
    semilogx(cf, ic_mean_clipped, 'LineWidth', 2); 
    xlabel('Frequency (Hz)')
    ylabel('Mean rate (sp/s)')
    title('Inferior colliculus response (mean rate)')
    xlim([0 8000]); 
    set(gca, 'XTick', [20 50 100 200 500 1000 2000 5000 8000]);
    set(gca, 'XTickLabel', {'20', '50', '100', '200', '500', '1k', '2k', '5k', '8k'});
    grid on;

end
