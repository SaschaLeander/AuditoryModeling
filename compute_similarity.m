%% Function computes similarity between two neurograms =======================

% written by Sam Jones & Sascha Muhlinghaus (12/08/2024)
% contact: saschamuhlinghaus@gmail.com

function [nsi, ssi] = compute_similarity(neurogram_ref, neurogram_degraded, flag, C1, C2, C3)

    % COMPUTE_SIMILARITY returns the similarity measures SSIM and NSIM for
    % two neurograms
        %
        %   Usage: compute_similarity(neurogram_ref, neurogram_degraded, flag, C1, C2, C3)
        %
        %   Input parameters:
        %       neurogram_ref:          neurogram output with normal cochelear parameters         
        %       neurogram_degraded:     neurogram output with degraded cochelear parameters
        %       flag:                   return 'ssi', 'nsi' or 'both' (string)
        %       C1:                     parameters C1 for luminance (float)
        %       C2:                     parameter C2 for contrast (float)
        %       C3:                     parameter C3 for structure (float)
        %   
        %   Output parameters:
        %       nsi:    Neurogram Similarity Index (0-1)
        %       ssi:    Structural Similarity Index (0-1)

        
    % Compute the neurogram similarity index measure (NSIM).
    if nargin < 4
        C1 = 0.01;
    end
    if nargin < 5
        C2 = 0.03;
    end
    if nargin < 6
        C3 = 0.03^(2);
    end
    neurogram_ref = normalize(neurogram_ref);
    neurogram_degraded = normalize(neurogram_degraded);

    mu_r = mean(neurogram_ref(:));
    mu_d = mean(neurogram_degraded(:));
    sigma_r = var(neurogram_ref(:)) + 1e-10;
    sigma_d = var(neurogram_degraded(:)) + 1e-10;
    sigma_rd = cov(neurogram_ref(:), neurogram_degraded(:));
    sigma_rd = sigma_rd(1, 2);

    alpha = 1;
    beta = 1;
    gamma = 1;

    % Edge case: If the variance is zero for both, return perfect similarity
    if sigma_r <= 1e-10 && sigma_d <= 1e-10
        nsi = 1; % NSIM should be 1 if both images are identical
        ssi = 1; % SSIM should be 1 if both images are identical
        return;
    end

    % Compute luminance term
    luminance = (2 * mu_r * mu_d + C1) / (mu_r^2 + mu_d^2 + C1);

    % Compute contrast term
    contrast = (sigma_rd + C2) / (sqrt(sigma_r * sigma_d) + C2);
    
    % Compute structure term 
    structure = (2 * sigma_r * sigma_d + C3) / (sigma_r^2 + sigma_d^2 + C3);

    switch flag
        case 'nsi'
            nsi = luminance * contrast;
            % Clamp the NSIM score to the range [0, 1] to avoid slight deviations above 1
            nsi = min(max(nsi, 0), 1);
            ssi = NaN;  % Return NaN as placeholder since only NSIM is computed

        case 'ssi'
            ssi = luminance^(alpha) * structure^(beta) * contrast^(gamma);
            % Clamp the SSIM score to the range [0, 1] to avoid slight deviations above 1
            ssi = min(max(ssi, 0), 1);
            nsi = NaN;  % Return NaN as placeholder since only SSIM is computed
        
        case 'both'
            % Compute and return both metrics 
            nsi = luminance * contrast;
            nsi = min(max(nsi, 0), 1);

            ssi = luminance^(alpha) * structure^(beta) * contrast^(gamma);
            ssi = min(max(ssi, 0), 1);
    end  
end

function normalized_neurogram = normalize(neurogram)
    % Normalize the neurogram to the range [0, 1].
    min_val = min(neurogram(:));
    max_val = max(neurogram(:));
    if max_val - min_val == 0
        normalized_neurogram = zeros(size(neurogram));
    else
        normalized_neurogram = (neurogram - min_val) / (max_val - min_val);
    end
end
