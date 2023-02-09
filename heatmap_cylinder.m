% Visualization flags for straight cylinder
form_factor_flag = false;
single_plot_flag = false;
animation_flag = false;

% Generate meshgrid of theta and phi values
mesh_size = 100;
theta = linspace(0, pi, mesh_size);
phi = linspace(0, 2 * pi, mesh_size);
[theta_mesh, phi_mesh] = meshgrid(theta, phi);

% Physical parameter input
lambda = 450 * 10^(-9);
L = 300 * 10^(-9);
a = 20 * 10^(-9);
waviness_period = 15 * 10^(-6);
waviness_amplitude = 3 * 10^(-6);

% Physical parameter computation
k = (2*pi)/lambda;

% Evaluate the heat function at each (theta, phi) point
h_straight = heat_func_straight(theta_mesh, phi_mesh, k, L, a);

% Convert spherical coordinates to Cartesian coordinates
x = sin(theta_mesh) .* cos(phi_mesh);
y = sin(theta_mesh) .* sin(phi_mesh);
z = cos(theta_mesh);

% Visualize building blocks of the form factor
if form_factor_flag
    figure('WindowState', 'maximized');
    visualize_form_factor_blocks(theta_mesh, phi_mesh, k, L, a);
end

% Plot the heat map
if single_plot_flag
    figure('WindowState', 'maximized');
    generate_heat_map_sphere(x, y, z, h_straight);
    title("(Form Factor)^2 , L / \lambda = " + num2str(round(L/lambda, 2)));
end

% Struct initialization and assignment for animation function input
mesh_parameters = struct;
mesh_parameters.theta_mesh_val = theta_mesh;
mesh_parameters.phi_mesh_val = phi_mesh;

physical_parameters = struct;
physical_parameters.lambda = lambda;
physical_parameters.k_val = k;
physical_parameters.L0_val = 0.5 * lambda;
physical_parameters.Lm_val = 10 * lambda;
physical_parameters.Lf_val = 150 * lambda;
physical_parameters.a_val = a;

animation_inputs = struct;
animation_inputs.mesh_parameters_val = mesh_parameters;
animation_inputs.physical_parameters_val = physical_parameters;
animation_inputs.animation_flag_val = animation_flag;

% Create animation of heat map
generate_heat_map_animation(animation_inputs);

% Define function for visualizing form factor building blocks
function visualize_form_factor_blocks(theta, phi, k, L, a)
    A = sin(theta) .* cos(phi);
    B = sin(theta) .* sin(phi);
    C = cos(theta) - 1;
    M = sqrt(B.^2 + C.^2);
    
    scalingFactor_transverse = k*a;
    scalingFactor_longitudinal = k*L;
    
    F_transverse = besselj(1, scalingFactor_transverse * M) ./ (scalingFactor_transverse * M);
    F_transverse(:,1) = 0.5;
    F_longitudinal = sinc(scalingFactor_longitudinal * A);
    F = 2 * (F_transverse .* F_longitudinal);
    
    x = sin(theta) .* cos(phi);
    y = sin(theta) .* sin(phi);
    z = cos(theta);
    
    subplot(3, 3, 1);
    surf(x, y, z, A.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    title("A^2 = (sin\theta cos\phi)^2")
    hold on;
    colormap hot;
    colorbar;
    clim([0, 1]);
    view(3);
    axis equal
    axes_plotter();
    
    subplot(3, 3, 2);
    surf(x, y, z, sinc(A).^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    title("sinc(A)^2")
    hold on;
    colormap hot;
    colorbar;
    clim([0, 1]);
    view(3);
    axis equal
    axes_plotter();
    
    subplot(3, 3, 3);
    surf(x, y, z, F_longitudinal, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    title("(Longitudinal Form Factor F_L)^2 = sinc(kLA)^2 , kL = " + num2str(round(k*L, 2)));
    hold on;
    colormap hot;
    colorbar;
    clim([0, 1]);
    view(3);
    axis equal
    axes_plotter();
    
    subplot(3, 3, 4);
    surf(x, y, z, B.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    title("B^2 = (sin\theta sin\phi)^2")
    hold on;
    colormap hot;
    colorbar;
    clim([0, 5]);
    view(3);
    axis equal
    axes_plotter();
    
    subplot(3, 3, 5);
    surf(x, y, z, C.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    title("C^2 = (cos\theta - 1)^2")
    hold on;
    colormap hot;
    colorbar;
    clim([0, 5]);
    view(3);
    axis equal
    axes_plotter();
    
    subplot(3, 3, 6);
    surf(x, y, z, M.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    title("M^2 = B^2 + C^2")
    hold on;
    colormap hot;
    colorbar;
    clim([0, 5]);
    view(3);
    axis equal
    axes_plotter();
    
    subplot(3, 3, 7);
    surf(x, y, z, (besselj(1, M)./M).^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    title("(bessel(1, M) / M)^2")
    hold on;
    colormap hot;
    colorbar;
    clim([0, 0.3]);
    view(3);
    axis equal
    axes_plotter();
    
    subplot(3, 3, 8);
    surf(x, y, z, F_transverse.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    title("(Transverse Form Factor F_T)^2 = (bessel(1, kaM) / kaM)^2 , ka = " + num2str(round(k*a, 2)));
    hold on;
    colormap hot;
    colorbar;
    clim([0, 0.3]);
    view(3);
    axis equal
    axes_plotter();
    
    subplot(3, 3, 9);
    surf(x, y, z, F.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    title("(Form Factor F)^2 = 4 * F_T * F_L");
    hold on;
    colormap hot;
    colorbar;
    clim([0, 1]);
    view(3);
    axis equal
    axes_plotter();

    sgtitle("Developing Intuition for Straight Cylinder Form factor");
  
end

% Define animation generating function
function generate_heat_map_animation(animation_inputs)
    mesh_parameters = animation_inputs.mesh_parameters_val;
    physical_parameters = animation_inputs.physical_parameters_val;
    animation_flag = animation_inputs.animation_flag_val;

    theta_mesh = mesh_parameters.theta_mesh_val;
    phi_mesh = mesh_parameters.phi_mesh_val;

    lambda = physical_parameters.lambda;
    k = physical_parameters.k_val;
    L0 = physical_parameters.L0_val;
    Lf = physical_parameters.Lf_val;
    Lm = physical_parameters.Lm_val;
    a = physical_parameters.a_val;

    x = sin(theta_mesh) .* cos(phi_mesh);
    y = sin(theta_mesh) .* sin(phi_mesh);
    z = cos(theta_mesh);

    if animation_flag
        vidObj = VideoWriter('heat_map_animation.mp4', 'MPEG-4');
        open(vidObj);
        figure('WindowState', 'maximized');
        numPts_array_start = 100;
        numPts_array_finish = 500;
        length_array_start = linspace(L0, Lm, numPts_array_start);
        length_array_finish = linspace(Lm, Lf, numPts_array_finish);
        length_array = horzcat(length_array_start, length_array_finish);
        for iter = 1:length(length_array)
            L = length_array(iter);
            h_straight = heat_func_straight(theta_mesh, phi_mesh, k, L, a);
            generate_heat_map_sphere(x, y, z, h_straight);
            title("(Form Factor)^2 , L / \lambda = " + num2str(round(L/lambda, 2)));
            pause(0.5);
            frame = getframe(gcf);
            writeVideo(vidObj, frame);
            clf;
        end
        close(vidObj);
        close(gcf);
    end

end

% Define heat map generating function on sphere
function generate_heat_map_sphere(x, y, z, h)    
    surf(x, y, z, h.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    hold on;
    colormap hot;
    colorbar;
    clim([0, 1]);
    view(3);
    axis equal
    
    axes_plotter();
    arrows_plotter();
    cylinder_plotter();
end

% Define the heat function for a straight cylinder
function h_straight = heat_func_straight(theta, phi, k, L, a)
    A = sin(theta) .* cos(phi);
    B = sin(theta) .* sin(phi);
    C = cos(theta) - 1;
    M = sqrt(B.^2 + C.^2);
    
    scalingFactor_transverse = k*a;
    scalingFactor_longitudinal = k*L;
    
    F_transverse = besselj(1, scalingFactor_transverse * M) ./ (scalingFactor_transverse * M);
    F_transverse(:,1) = 0.5;
    
    F_longitudinal = sinc(scalingFactor_longitudinal * A);
    
    F = 2 * (F_transverse .* F_longitudinal);
    h_straight = F;
end

% Define axes plotter
function axes_plotter
    line([0, 1.5], [0, 0], [0, 0], 'Color', 'red', 'LineWidth', 2);
    line([0, 0], [0, 1.5], [0, 0], 'Color', 'green', 'LineWidth', 2);
    line([0, 0], [0, 0], [0, 1.5], 'Color', 'blue', 'LineWidth', 2);
end

% Define cylinder plotter
function cylinder_plotter
    [X, Y, Z] = cylinder(0.1, 100);
    cyl_pltX = 0.4*(Z-0.5);
    cyl_pltY = X;
    cyl_pltZ = Y;
    c = surf(cyl_pltX, cyl_pltY, cyl_pltZ, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 1);
end

% Define arrows plotter
function arrows_plotter
    num_arrows = 6;
    arrow_length = 0.35;
    arrow_base = repmat([0, 0, -1.5], num_arrows, 1);
    arrow_base(:,1) = linspace(-0.2, 0.2, num_arrows)';
    arrow_direction = repmat([0, 0, arrow_length], num_arrows, 1);
    quiver3(arrow_base(:,1), arrow_base(:,2), arrow_base(:,3), arrow_direction(:,1), arrow_direction(:,2), arrow_direction(:,3), 0, 'Color', 'black', 'LineWidth', 0.5);
end