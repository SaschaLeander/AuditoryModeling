% written by Sam Jones and Sascha Muhlinghaus(30/08/2024)
% contact: saschamuhlinghaus@gmail.com

%% Function plots two auditory stimuli as a comparison=======================
function plot_neurogram_comparison(stimulus, r_mean, ihc, ic_sout_BE, cn_sout_contra, Fs, cf, stimulus2, r_mean2, ihc2, ic_sout_BE2, cn_sout_contra2, varArray, nsimArray, ssimArray)

    
    % PLOT_NEUROGRAM_COMPARISON plots auditory modelling output as stimulus comparison
    % 
    %   Usage: plot_neurogram_comparison(stimulus, r_mean, ihc, ic_sout_BE, cn_sout_contra, Fs, cf, stimulus2, r_mean2, ihc2, ic_sout_BE2, cn_sout_contra2, varArray, nsimArray, ssimArray)
    %
    %   Input parameters: (stimulus parameters must be provided for both stimuli)
    %       stimulus:           generated audio stimulus
    %       r_mean:             mean firing rates of AN, CN, IC
    %       ihc:                inner hair cell firing
    %       ic_sout_BE:         inferior colliculus output
    %       cn_sout_contra:     cochlear nucleus output 
    %
    %       Fs:                 sampling frequency
    %       cf:                 center frequencies
    %       varArray:           dependent variable (freqeuncy/dB SPL)
    %       nsimArray:          neurogram similarity index
    %       ssimArray:          structural similarity index  

    % re-size matrices, accounting for zero-padding within the model
    [sizeN, sizeM] = size(r_mean);
    ihc = ihc(1:sizeN, 1:sizeM);
    ic_sout_BE = ic_sout_BE(1:sizeN, 1:sizeM);
    cn_sout_contra = cn_sout_contra(1:sizeN, 1:sizeM);

    [sizeN2, sizeM2] = size(r_mean2);
    ihc2 = ihc2(1:sizeN2, 1:sizeM2);
    ic_sout_BE2 = ic_sout_BE2(1:sizeN2, 1:sizeM2);
    cn_sout_contra2 = cn_sout_contra2(1:sizeN2, 1:sizeM2);
    
    % define common time vector
    t = (1:length(r_mean))*1/Fs;
    t2 = (1:length(r_mean2))*1/Fs;

    % plot figures
    figure;
    
    % raw stimulus
    subplot(5, 5, 1)
    t_stimulus = 1/Fs:1/Fs:length(r_mean)/Fs;
    minLength = min(length(t_stimulus), length(stimulus));
    t_stimulus = t_stimulus(1:minLength);
    stimulus = stimulus(1:minLength);
    plot(t_stimulus, stimulus);
    title('Stimulus waveform (reference)');
    xlabel('Time (s)');
    ylabel('Amplitude (Pa)');
    grid on;

    % raw stimulus
    subplot(5, 5, 2)
    t_stimulus = 1/Fs:1/Fs:length(r_mean2)/Fs;
    minLength = min(length(t_stimulus), length(stimulus2));
    t_stimulus = t_stimulus(1:minLength);
    stimulus2 = stimulus2(1:minLength);
    plot(t_stimulus, stimulus2);
    title('Stimulus waveform (altered)');
    xlabel('Time (s)');
    ylabel('Amplitude (Pa)');
    grid on;

    % inner hair cell response
    subplot(5, 5, 3)
    plot(t, ihc);
    xlim([0 (length(stimulus)/Fs)*1.2]) %  scale x-axis by 1.2 * stimulus length
    title('Inner hair cell response (reference)')
    xlabel('Time (s)')
    ylabel('Volts')
    grid on

    % inner hair cell response 2
    subplot(5, 5, 4)
    plot(t2, ihc2);
    xlim([0 (length(stimulus)/Fs)*1.2]) %  scale x-axis by 1.2 * stimulus length
    title('Inner hair cell response (altered)')
    xlabel('Time (s)')
    ylabel('Volts')
    grid on

    % auditory nerve mean rate by CF
    r_clipped = r_mean(1:length(stimulus),:); 
    r_mean_clipped = mean(r_clipped, 1);

    r_clipped2 = r_mean2(1:length(stimulus2),:); 
    r_mean_clipped2 = mean(r_clipped2, 1);

    subplot(5, 5, 5);
    semilogx(cf, r_mean_clipped, 'b-', 'LineWidth', 2);
    hold on; % Hold on to add the second line to the same plot
    semilogx(cf, r_mean_clipped2, 'r-', 'LineWidth', 2);
    hold off;
    xlabel('Frequency (Hz)')
    ylabel('Mean rate (sp/s)')
    title('Auditory nerve response (mean rate)')
    xlim([0 8000]); 
    legend('reference', 'altered');
    set(gca, 'XTick', [20 50 100 200 500 1000 2000 5000 8000]);
    set(gca, 'XTickLabel', {'20', '50', '100', '200', '500', '1k', '2k', '5k', '8k'});
    grid on;

    % frequency spectrum subplot
    subplot(5, 5, 6)

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
    title('Stimulus frequency spectrum (reference)')
    set(gca, 'XTick', [20 50 100 200 500 1000 2000 5000 8000]);
    set(gca, 'XTickLabel', {'20', '50', '100', '200', '500', '1k', '2k', '5k', '8k'});
    grid on;

    % frequency spectrum subplot 2
    subplot(5, 5, 7)

    % calculate FFT
    n = length(stimulus2);
    frequencies = (0:n-1)*(Fs/n);
    magnitude = abs(fft(stimulus2))/n;
    
    % convert to dB SPL (assuming the reference pressure p0 = 20 µPa)
    referencePressure = 20e-6;
    magnitude_db_spl = 20 * log10(magnitude / referencePressure);
    
    % plot the magnitude in dB SPL by frequency
    semilogx(frequencies(1:floor(n/2)), magnitude_db_spl(1:floor(n/2)));
    xlim([125 8000]); 
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB SPL)')
    title('Stimulus frequency spectrum (altered)')
    set(gca, 'XTick', [20 50 100 200 500 1000 2000 5000 8000]);
    set(gca, 'XTickLabel', {'20', '50', '100', '200', '500', '1k', '2k', '5k', '8k'});
    grid on;

    % auditory nerve response  1
    subplot(5, 5, 8)
    plot(t, r_mean)
    title('Auditory nerve response')
    xlim([0 length(stimulus)/Fs])
    xlabel('Time (s)')
    ylabel('Instantaneous rate (sp/s) (reference)')
    grid on;

    % auditory nerve response  2
    subplot(5, 5, 9)
    plot(t, r_mean2)
    title('Auditory nerve response')
    xlim([0 length(stimulus)/Fs])
    xlabel('Time (s)')
    ylabel('Instantaneous rate (sp/s) (altered)')
    grid on;
    
    % cochlear nucleus response mean rate by CF
    cn_clipped = cn_sout_contra(1:length(stimulus),:);
    cn_mean_clipped = mean(cn_clipped, 1);

    cn_clipped2 = cn_sout_contra2(1:length(stimulus2),:);
    cn_mean_clipped2 = mean(cn_clipped2, 1);
    subplot(5, 5, 10);
    semilogx(cf, cn_mean_clipped, 'b-', 'LineWidth', 2); 
    hold on;
    semilogx(cf, cn_mean_clipped2, 'r-', 'LineWidth', 2); 
    hold off;
    xlabel('Frequency (Hz)')
    ylabel('Mean rate (sp/s)')
    title('Cochlear nucleus response (mean rate)')
    legend('reference', 'altered');
    xlim([0 8000]); 
    set(gca, 'XTick', [20 50 100 200 500 1000 2000 5000 8000]);
    set(gca, 'XTickLabel', {'20', '50', '100', '200', '500', '1k', '2k', '5k', '8k'});
    grid on;

    % auditory nerve response CF 1
    subplot(5, 5, 11)
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
    title('Auditory nerve response (reference)')
    grid on

    % auditory nerve response CF 2
    subplot(5, 5, 12)
    imagesc(t2, log10(cf/1e3), r_mean2');
    axis xy;
    yticks = [0.125 0.5 2 8 16];
    set(gca, 'ytick', log10(yticks));
    set(gca, 'yticklabel', yticks);
    ylabel('CF (kHz)')
    xlabel('Time (s)')
    xlim([0 length(stimulus)/Fs])
    hcb = colorbar;
    set(get(hcb, 'ylabel'), 'string', 'Spikes / second');
    title('Auditory nerve response (altered)')
    grid on;

    % cochlear nucleus response 1
    subplot(5, 5, 13)
    plot(t, cn_sout_contra)
    title('Cochlear nucleus response')
    xlim([0 length(stimulus)/Fs])
    xlabel('Time (s)')
    ylabel('Instantaneous rate (sp/s) (reference)')
    grid on;

    % cochlear nucleus response 2
    subplot(5, 5, 14)
    plot(t, cn_sout_contra2)
    title('Cochlear nucleus response')
    xlim([0 length(stimulus)/Fs])
    xlabel('Time (s)')
    ylabel('Instantaneous rate (sp/s) (altered)')
    grid on;

    % inferior colliculus response mean rate by CF
    ic_clipped = ic_sout_BE(1:length(stimulus),:);
    ic_mean_clipped = mean(ic_clipped, 1); 

    ic_clipped2 = ic_sout_BE2(1:length(stimulus2),:);
    ic_mean_clipped2 = mean(ic_clipped2, 1);
    subplot(5, 5, 15);
    semilogx(cf, ic_mean_clipped, 'b-', 'LineWidth', 2); 
    hold on;
    semilogx(cf, ic_mean_clipped2, 'r-', 'LineWidth', 2); 
    hold off;
    xlabel('Frequency (Hz)');
    ylabel('Mean rate (sp/s)');
    legend('reference', 'altered');
    title('Inferior colliculus response (mean rate)');
    xlim([0 8000]); 
    set(gca, 'XTick', [20 50 100 200 500 1000 2000 5000 8000]);
    set(gca, 'XTickLabel', {'20', '50', '100', '200', '500', '1k', '2k', '5k', '8k'});
    grid on;

    % cochlear nucleus response CF 1
    subplot(5, 5, 16)
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
    title('Cochlear nucleus response (reference)')
    grid on;

    % cochlear nucleus response CF 2
    subplot(5, 5, 17)
    imagesc(t2, log10(cf/1e3), cn_sout_contra2');
    axis xy;
    yticks = [0.125 0.5 2 8 16];
    set(gca, 'ytick', log10(yticks));
    set(gca, 'yticklabel', yticks);
    ylabel('CF (kHz)')
    xlabel('Time (s)')
    xlim([0 length(stimulus2)/Fs])
    hcb = colorbar;
    set(get(hcb, 'ylabel'), 'string', 'Spikes / second');
    title('Cochlear nucleus response (altered)')
    grid on;

    % inferior colliculus response 1
    subplot(5, 5, 18)
    plot(t, ic_sout_BE)
    xlim([0 length(stimulus)/Fs])
    xlabel('Time (s)')
    ylabel('Instantaneous rate (sp/s)')
    title('Inferior colliculus response (reference)')
    grid on;

    % inferior colliculus response 2
    subplot(5, 5, 19)
    plot(t2, ic_sout_BE2)
    xlim([0 length(stimulus2)/Fs])
    xlabel('Time (s)')
    ylabel('Instantaneous rate (sp/s)')
    title('Inferior colliculus response (altered)')
    grid on;
    
    subplot(5, 5, 20);
    plot(varArray, nsimArray, 'r-', 'LineWidth', 2); % Plot the first line in red
    hold on; % Hold on to add the second line to the same plot
    plot(varArray, ssimArray, 'b-', 'LineWidth', 2); % Plot the first line in green
    hold off;  
    % Add labels and title
    xlabel('Signal Frequency/Amplitude [Hz/dB SPL]');
    ylabel('Similarity Metric [0 - 1]');
    legend('reference', 'altered');
    title('Neurogram Similarity');
    % Add legend
    legend('NSIM','SSIM'); 
    % Show grid
    grid on;

    % inferior colliculus response CF 1
    subplot(5, 5, 21)
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
    title('Inferior colliculus response (reference)')
    grid on;

    % inferior colliculus response CF 2
    subplot(5, 5, 22)
    imagesc(t2, log10(cf/1e3), ic_sout_BE2');
    axis xy;
    yticks = [0.125 0.5 2 8 16];
    set(gca, 'ytick', log10(yticks));
    set(gca, 'yticklabel', yticks);
    ylabel('CF (kHz)')
    xlabel('Time (s)')
    xlim([0 length(stimulus2)/Fs])
    hcb = colorbar;
    set(get(hcb, 'ylabel'), 'string', 'Spikes / second');
    title('Inferior colliculus response (altered)')
    grid on;

end
