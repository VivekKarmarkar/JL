% Waviness parameters
waviness_period_start_muM = 2;
waviness_period_end_muM = 9;
waviness_period_array_muM = linspace(waviness_period_start_muM, waviness_period_end_muM, 200);
waviness_period_array = (10^(-6)) * waviness_period_array_muM;

% Grid
grid_size = 1000;
phi = linspace(0, pi, grid_size);
phi_deg = rad2deg(phi);

% Intensity Plot vs angle Phi animation with respect to waviness period
vidObj = VideoWriter('intensity_circlePhiPlot_animation.mp4', 'MPEG-4');
open(vidObj);
figure('WindowState','maximized');
for iter = 1:length(waviness_period_array)
    waviness_period = waviness_period_array(iter);
    [F, theta_deg] = WavyFormFactor(grid_size, phi, waviness_period);
    generate_plot(waviness_period, theta_deg, phi_deg, F);
    pause(0.5);
    frame = getframe(gcf);
    writeVideo(vidObj, frame);
    clf;
end
close(vidObj);
close(gcf);

% Plot Function
function generate_plot(waviness_period, theta_deg, phi_deg, F)
    plot(phi_deg, F.^2, 'r', 'LineWidth', 1);
    xlabel('\phi (degrees)', 'FontWeight','bold');
    ylabel('I (\phi)', 'FontWeight','bold');
    grid;
    title("Single Fiber Intensity as (Form Factor)^2 on a circle at fixed \theta = " + num2str(theta_deg) + " degrees , \Lambda = " + num2str(waviness_period/(10 ^ (-6))) + " \mum");
    ylim([0 1.2]);
end

% Form Factor Function
function [F, theta_deg] = WavyFormFactor(grid_size, phi, waviness_period)
    % Optical parameters
    lambda = 450 * 10^(-9);
    
    % Cylinder dimensions
    a = 20 * 10^(-9);
    L = 300 * 10^(-6);
    
    % Waviness parameters
    waviness_amplitude = 3 * 10^(-6);

    % Theta
    theta_deg = 0.5;
    theta = deg2rad(theta_deg);
    
    % Wavenumber related computations
    k = (2 * pi) / lambda;
    P = double(int32(L / (2 * waviness_period)));
    normalization_factor = waviness_period / L;
    wavelength_ratio = waviness_period / lambda;
    
    % Scaling facors
    scalingFactor_transverse = k*a;
    scalingFactor_longitudinal = k*L;
    scalingFactor_waviness = k*waviness_amplitude;
    
    % Vectorized computation of form factor building blocks
    A = sin(theta) * cos(phi);
    B = sin(theta) * sin(phi);
    C = cos(theta) - 1;
    M = sqrt(B.^2 + C^2);
    
    % Transverse Form Factor
    FT = besselj(1, scalingFactor_transverse * M) ./ (scalingFactor_transverse * M);
    FT(:,1) = 0.5;
    
    % Computing nu
    nu = wavelength_ratio * A;
    
    % Anger Function non-vectorized computations
    FA = nan(1, grid_size);
    FA_arg = scalingFactor_waviness * (cos(theta) - 1);
    for m_iter = 1:grid_size
        FA(1, m_iter) = AngerFunc(nu(1, m_iter), FA_arg);
        disp(m_iter);
    end
    
    % Wavy Cylinder form factor from paper
    FR = (sin(2 * pi * P * nu) ./ sin(pi * nu));
    FRN = normalization_factor * FR;
    FL = FA .* FRN;
    F = 2 * (FT .* FL);

end

% Anger Function Definition
function j = AngerFunc(nu, z)
    j = (1/pi) * integral(@(t) cos(nu .* t - z .* sin(t)), 0, pi);
end