%% Function stores results and stimulus parameters from neurogram similarity comparison =======================

% written by Sascha Muhlinghaus (16/08/2024)
% contact: saschamuhlinghaus@gmail.com

function manage_data(index, intensity, frequency, ihc, ohc, bm, task, competitor, nsim, ssim)
    % Define the filename
    filename = 'Test_AuditoryData.csv';

    % Check if the file exists
    if isfile(filename)
        % Load existing dataset
        opts = detectImportOptions(filename);
        data = readtable(filename, opts);
        
        % Prepare new row to add, ensuring numeric values are kept as they are
        newRow = {index, intensity, frequency, ihc, ohc, bm, task, competitor, nsim, ssim};
        
        % Append the new row to the dataset
        data = [data; newRow];
    else
        % Define headers
        headers = {'Stimulus Index', 'dB SPL', 'Frequency (Hz)', 'ihc', 'ohc', 'bm', 'IMAP task', ...
                   'Competitor Stimulus (Hz/dB SPL)', 'NSIM', 'SSIM'};
        
        % Prepare the first row of data
        firstRow = {index, intensity, frequency, ihc, ohc, bm, task, competitor, nsim, ssim};
        
        % Create a table with headers and the first row
        data = cell2table(firstRow, 'VariableNames', headers);
    end
    
    % Write the table to the file
    writetable(data, filename);
end

