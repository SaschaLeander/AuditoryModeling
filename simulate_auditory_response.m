function [r_mean, ihc, ic_sout_BE, cn_sout_contra] = simulate_auditory_response(stimulus, Fs, cf, bm, fiberType, numH, numM, numL, fiber_num, cohc, cihc, BMF)
   %% run models =============================================================

    % inner hair cells and synapse out 
    
    [r_mean, psth, ihc] = zilany2014(stimulus, Fs, cf, 'bm', bm, 'fiberType', fiberType, ...
        'numH', numH, 'numM', numM, 'numL', numL, 'nrep', fiber_num, ...
        'cohc', cohc, 'cihc', cihc);
    
    % cochlear nucleus and inferior colliculus out
    
    [ic_sout_BE, ic_sout_BS, cn_sout_contra] = carney2015(r_mean, BMF, Fs);
    
end

