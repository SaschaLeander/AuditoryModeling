function return_metric = compute_similarity(neurogram_ref, neurogram_degraded, flag, C1, C2, C3)
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

    % Check if variances are zero and handle the edge case
    if sigma_r <= 1e-10 && sigma_d <= 1e-10
        return_metric = 0;
        return;
    end

    % Compute luminance term
    luminance = (2 * mu_r * mu_d + C1) / (mu_r^2 + mu_d^2 + C1);

    % Compute contrast term
    contrast = (sigma_rd + C2) / (sqrt(sigma_r * sigma_d) + C2);
    
    % Compute structure term 
    structure = (2*sigma_r * sigma_d + C3) / (sigma_r^(2) + sigma_d^(2) + C3);

    switch flag
        case 'nsi'
            
            return_metric = luminance * contrast;
            % Clamp the NSIM score to the range [0, 1] to avoid slight deviations above 1
            return_metric = min(max(return_metric, 0), 1);

        case 'ssi'

            alpha = 1;
            beta = 1;
            gamma = 1;
    
            return_metric = luminance^(alpha) * structure^(beta) * contrast^(gamma);
            % Clamp the SSIM score to the range [0, 1] to avoid slight deviations above 1
            return_metric = min(max(return_metric, 0), 1);
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
