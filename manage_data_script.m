function manage_data(index, stimulus, ihc, ohc, bm, task, target, competitor_index, nsim, ssim)
    % Define the filename
    filename = 'Test_AuditoryData.csv';
    
    % Check if the file exists
    if isfile(filename)
        % Load existing dataset
        opts = detectImportOptions(filename);
        data = readtable(filename, opts);
        
        % Prepare new row to add, ensuring numeric values are kept as they are
        newRow = {index, {stimulus}, ihc, ohc, bm, task, target, competitor_index, nsim, ssim};
        
        % Append the new row to the dataset
        data = [data; newRow];
    else
        % Define headers
        headers = {'index', 'stimulus', 'ihc', 'ohc', 'bm', 'task', ...
                   'target', 'competitor_index', 'nsim', 'ssim'};
        
        % Prepare the first row of data
        firstRow = {index, {stimulus}, ihc, ohc, bm, task, target, competitor_index, nsim, ssim};
        
        % Create a table with headers and the first row
        data = cell2table(firstRow, 'VariableNames', headers);
    end
    
    % Write the table to the file
    writetable(data, filename);
end


% Call the function to manage the dataset

index = 1;
stimulus = stimulus2;
ihc = 1;
ohc = 2;
bm = 3;
task = 'classification';
target = 1;
competitor_index = 2;
nsim = 0.89;
ssim = 0.9;

manage_data(index, stimulus, ihc, ohc, bm, task, target, competitor_index, nsim, ssim);

