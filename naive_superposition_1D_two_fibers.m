% Separation between fibers
separatation_small = true;

% Grid
theta_deg = 20;
theta_start = 0;
theta_end = deg2rad(theta_deg);
grid_size = 10000;
theta = linspace(theta_start, theta_end, grid_size);

% Waviness Inputs
waviness_period_one = 5 * 10^(-6);
if separatation_small
    separatation_str = "small";
    waviness_period_two = 5.1 * 10^(-6);
else
    separatation_str = "large";
    waviness_period_two = 9 * 10^(-6);
end

% Coordinate values
x = sin(theta);
y = zeros(size(x));
z = cos(theta);

% Grid struct
physical_grid = struct;
physical_grid.theta_val = theta;
physical_grid.theta_deg_val = theta_deg;
physical_grid.x_val = x;
physical_grid.y_val = y;
physical_grid.z_val = z;

% Form Factor computation
F1 = compute_1D_form_factor(physical_grid, waviness_period_one);
F2 = compute_1D_form_factor(physical_grid, waviness_period_two);

% Fiber Info struct
fiber_info = struct;
fiber_info.separation_str_val = separatation_str;
fiber_info.waviness_period_one_val = waviness_period_one;
fiber_info.waviness_period_two_val = waviness_period_two;
fiber_info.F1_val = F1;
fiber_info.F2_val = F2;

% Visualize Superposition Results
visualize_superposition_results(physical_grid, fiber_info);

% Visualizing Results Function
function visualize_superposition_results(physical_grid, fiber_info)
    
    x = physical_grid.x_val;
    theta_deg = physical_grid.theta_deg_val;

    waviness_period_one = fiber_info.waviness_period_one_val;
    waviness_period_two = fiber_info.waviness_period_two_val;

    F1 = fiber_info.F1_val;
    F2 = fiber_info.F2_val;
    F = F1 + F2;

    separatation_str = fiber_info.separation_str_val;

    figure(WindowState="maximized");
    
    subplot(2,2,1);
    plot(x, F1.^2, 'r', 'DisplayName', "\Lambda_{1} = " + num2str(waviness_period_one / (10^(-6))) + " \mum");
    hold on;
    plot(x, F2.^2, 'b', 'DisplayName', "\Lambda_{2} = " + num2str(waviness_period_two / (10^(-6))) + " \mum");
    plot(x, F.^2, 'k', 'DisplayName', '\Lambda_{1} + \Lambda_{2}');
    grid;
    title("(Form Factor)^2");
    
    subplot(2,2,3);
    plot(x, F1.^2, 'r', 'DisplayName', "\Lambda_{1} = " + num2str(waviness_period_one / (10^(-6))) + " \mum");
    hold on;
    plot(x, F2.^2, 'b', 'DisplayName', "\Lambda_{2} = " + num2str(waviness_period_two / (10^(-6))) + " \mum");
    plot(x, F.^2, 'k', 'DisplayName', '\Lambda_{1} + \Lambda_{2}');
    grid;
    xlim([0 0.01]);
    xlabel('c_{||}', 'FontWeight', 'bold');
    title("(Form Factor)^2, Zeroth Peak, Zoomed In");
    
    subplot(2,2,2);
    plot(x, F1.^2, 'r', 'DisplayName', "\Lambda_{1} = " + num2str(waviness_period_one / (10^(-6))) + " \mum");
    hold on;
    plot(x, F2.^2, 'b', 'DisplayName', "\Lambda_{2} = " + num2str(waviness_period_two / (10^(-6))) + " \mum");
    plot(x, F.^2, 'k', 'DisplayName', '\Lambda_{1} + \Lambda_{2}');
    legend;
    grid;
    if separatation_str == "small"
        xlim([0.081 0.28]);
    else
        xlim([0.045 0.28]);
    end
    title("(Form Factor)^2 , Zoomed In");
    
    subplot(2,2,4);
    plot(x, F1.^2, 'r', 'DisplayName', "\Lambda_{1} = " + num2str(waviness_period_one / (10^(-6))) + " \mum");
    hold on;
    plot(x, F2.^2, 'b', 'DisplayName', "\Lambda_{2} = " + num2str(waviness_period_two / (10^(-6))) + " \mum");
    plot(x, F.^2, 'k', 'DisplayName', '\Lambda_{1} + \Lambda_{2}');
    grid;
    if separatation_str == "small"
        xlim([0.081 0.1]);
    else
        xlim([0.045 0.1]);
    end
    title("(Form Factor)^2 , First Peak , Zoomed In");

    sgtitle_str = "Naively Developing Intuition for 1D Superposition of Two Fibers (with " + separatation_str + " \Lambda separation) at:  \phi = 0 and 0 \leq \theta \leq " + num2str(theta_deg) + char(176);
    sgtitle(sgtitle_str);

end
 
% 1D Form Factor Computation Function
function form_factor = compute_1D_form_factor(physical_grid, waviness_period)
    % Optical parameters
    lambda = 450 * 10^(-9);
    
    % Cylinder dimensions
    a = 20 * 10^(-9);
    L = 300 * 10^(-6);
    
    % Waviness parameters
    waviness_amplitude = 3 * 10^(-6);

    % Physical Grid
    theta = physical_grid.theta_val;
    x = physical_grid.x_val;

    % Wavenumber related computations
    k = (2 * pi) / lambda;
    P = double(int32(L / (2 * waviness_period)));
    normalization_factor = waviness_period / L;
    wavelength_ratio = waviness_period / lambda;

     % Scaling facors
    scalingFactor_transverse = k*a;
    scalingFactor_waviness = k*waviness_amplitude;
    
    % Computing nu by scaling x
    nu = wavelength_ratio * x;

    % Anger Function Factor
    FA_arg = scalingFactor_waviness * (cos(theta) - 1);
    FA = nan(size(nu));
    for n_iter = 1: length(nu)
        FA(1, n_iter) = AngerFunc(nu(1, n_iter), FA_arg(1, n_iter));
    end

    % Vectorized computation of form factor building blocks
    A = sin(theta);
    B = zeros(size(A));
    C = cos(theta) - 1;
    M = sqrt(B.^2 + C.^2);
    
    % Transverse Form Factor
    FT = besselj(1, scalingFactor_transverse * M) ./ (scalingFactor_transverse * M);
    FT(:,1) = 0.5;

    % Wavy Cylinder form factor from paper
    FR = (sin(2 * pi * P * nu) ./ sin(pi * nu));
    FR(:,1) = 1;
    FRN = normalization_factor * FR;
    FRN(:,1) = 1;
    FL = FA .* FRN;
    F = 2 * (FT .* FL);

    % Assign form factor to return variable
    form_factor = F;

end

% Anger Function Definition
function j = AngerFunc(nu, z)
    j = (1/pi) * integral(@(t) cos(nu .* t - z .* sin(t)), 0, pi);
end