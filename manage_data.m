%% Function stores results and stimulus parameters from neurogram similarity comparison =======================

% Author: Sascha Muhlinghaus 
% contact: saschamuhlinghaus@gmail.com

function manage_data(index, numCF, intensity, frequency, cihc, cohc, bm, task, competitor, nsim, ssim)

    % MANAGE_DATA saves auditory modeling parameters and results to a csv file
    %
    %   Usage: manage_data(index, intensity, frequency, ihc, ohc, bm, task, competitor, nsim, ssim)
    %
    %   Input parameters:
    %       index:      stimulus sample index
    %       intensity:  stimulus level (dB SPL)
    %       frequency:  stimulus frequency (Hz)
    %       cihc:       innerhair cell function scaling factor: 1 denotes normal function 
    %       cohc:       outerhair cell function scaling factor: 1 denotes normal function
    %       bm:         basilar membrane tuning, 1 = normal, > 1 = broader 
    %       task:       IMAP task 
    %       competitor: reference stimulus with baseline dep. var. value
    %       nsim:       Neurogram Similarity Index 
    %       ssim:       Structural Similarity Index

    % Define the filename
    filename = 'Test_AuditoryData.csv';

    % Check if the file exists
    if isfile(filename)
        % Load existing dataset
        opts = detectImportOptions(filename);
        data = readtable(filename, opts);
        
        % Prepare new row to add, ensuring numeric values are kept as they are
        newRow = {index, numCF, intensity, frequency, ihc, ohc, bm, task, competitor, nsim, ssim};
        
        % Append the new row to the dataset
        data = [data; newRow];
    else
        % Define headers
        headers = {'Stimulus Index', 'number of fibres', 'dB SPL', 'Frequency (Hz)', 'cihc', 'cohc', 'bm', 'IMAP task', ...
                   'Competitor Stimulus (Hz/dB SPL)', 'NSIM', 'SSIM'};
        
        % Prepare the first row of data
        firstRow = {index, numCF, intensity, frequency, cihc, cohc, bm, task, competitor, nsim, ssim};
        
        % Create a table with headers and the first row
        data = cell2table(firstRow, 'VariableNames', headers);
    end
    
    % Write the table to the file
    writetable(data, filename);
end

