function [r_mean, ihc, ic_sout_BE, cn_sout_contra] = simulate_auditory_response(stimulus, Fs, cf, bm, fiberType, numH, numM, numL, fiber_num, cohc, cihc, BMF, filename)
   %% run models =============================================================

    % inner hair cells and synapse out 
    
    [r_mean, psth, ihc] = zilany2014(stimulus, Fs, cf, 'bm', bm, 'fiberType', fiberType, ...
        'numH', numH, 'numM', numM, 'numL', numL, 'nrep', fiber_num, ...
        'cohc', cohc, 'cihc', cihc);
    
    % cochlear nucleus and inferior colliculus out
    
    [ic_sout_BE, ic_sout_BS, cn_sout_contra] = carney2015(r_mean, BMF, Fs);
    
    %% save output ============================================================

    function manage_dataset(stimulus, ic_sout_BE, filename)
        % Check if the file exists
        if isfile(filename)
            % Load the existing dataset
            load(filename, 'dataTable');
            
            % Create a new row with the current variables as cell arrays
            newRow = {stimulus, ic_sout_BE};
            
            % Append the new row to the existing table
            newRowTable = cell2table(newRow, 'VariableNames', {'stimulus', 'ic_sout_BE'});
            dataTable = [dataTable; newRowTable];
        else
            % Create a new table with the given variables as cell arrays
            dataTable = table({stimulus}, {ic_sout_BE}, 'VariableNames', {'stimulus', 'ic_sout_BE'});
        end
        
        % Save the updated table to the file
        save(filename, 'dataTable');
    end

    % Call the function to manage the dataset
    %manage_dataset(stimulus, ic_sout_BE , filename);
    
end

