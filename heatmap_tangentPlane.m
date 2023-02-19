% Optical parameters
lambda = 450 * 10^(-9);

% Cylinder dimensions
a = 20 * 10^(-9);
L = 300 * 10^(-6);

% Waviness parameters
waviness_amplitude = 3 * 10^(-6);
waviness_period = 5 * 10^(-6);

% Grid
plane_size = 0.35;
grid_size = 1000;
x = linspace(-plane_size, plane_size, grid_size)';
y = linspace(-plane_size, plane_size, grid_size);
z = ones(grid_size, grid_size);
[x_mesh, y_mesh] = meshgrid(x, y);
z_mesh = z;

% Wavenumber related computations
k = (2 * pi) / lambda;
P = double(int32(L / (2 * waviness_period)));
normalization_factor = waviness_period / L;
wavelength_ratio = waviness_period / lambda;

% Scaling facors
scalingFactor_transverse = k*a;
scalingFactor_longitudinal = k*L;
scalingFactor_waviness = k*waviness_amplitude;

% Create Cartesian mesh struct
cartesian_mesh = struct;
cartesian_mesh.x_mesh_val = x_mesh;
cartesian_mesh.y_mesh_val = y_mesh;
cartesian_mesh.z_mesh_val = z_mesh;
cartesian_mesh.grid_size_val = grid_size;

% Extract spherical coordinate angles in radians from Cartesian coordinates
[theta, phi] = custom_CartToSph(cartesian_mesh);

% Converting the angles to degrees
theta_deg = rad2deg(theta);
phi_deg = rad2deg(phi);

% Vectorized computation of form factor building blocks
A = sin(theta) .* cos(phi);
B = sin(theta) .* sin(phi);
C = cos(theta) - 1;
M = sqrt(B.^2 + C.^2);

% Transverse Form Factor
FT = besselj(1, scalingFactor_transverse * M) ./ (scalingFactor_transverse * M);
FT(:,1) = 0.5;

% Computing nu
nu = wavelength_ratio * A;

% Anger Function non-vectorized computations
FA = nan(grid_size, grid_size);
FA_arg = scalingFactor_waviness * (cos(theta) - 1);
for n_iter = 1:grid_size
    disp(n_iter);
    for m_iter = 1:grid_size
        FA(n_iter, m_iter) = AngerFunc(nu(n_iter, m_iter), FA_arg(n_iter, m_iter));
        disp(m_iter);
    end
end

% Wavy Cylinder form factor from paper
FR = (sin(2 * pi * P * nu) ./ sin(pi * nu));
FRN = normalization_factor * FR;
FL = FA .* FRN;
F = 2 * (FT .* FL);

% Form Factor Visualization
figure(WindowState="maximized");

subplot(4,2,1);
surf(x_mesh, y_mesh, z_mesh, FR, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
colormap hot;
colorbar;
view(2);
title("Ratio of Sines Factor F_R");

subplot(4,2,3);
surf(x_mesh, y_mesh, z_mesh, FRN.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
colormap hot;
colorbar;
clim([0 1]);
view(2);
title("Ratio of Sines Factor (normalized) , squared  (F_{RN})^2");

subplot(4,2,5);
surf(x_mesh, y_mesh, z_mesh, FA.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
colormap hot;
colorbar;
clim([0 1]);
view(2);
title("Anger function , squared  (F_A)^2");

subplot(4,2,7);
surf(x_mesh, y_mesh, z_mesh, FL.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
colormap hot;
colorbar;
clim([0 1]);
view(2);
title("Longitudinal Form Factor , squared  (F_{L})^2");
ylabel('Lc_{\perp}', 'FontWeight','bold');
xlabel('c_{||}', 'FontWeight','bold');

subplot(4,2,2);
surf(x_mesh, y_mesh, z_mesh, FL.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
colormap hot;
colorbar;
clim([0 1]);
view(2);
title("Longitudinal Form Factor , squared  (F_{L})^2");

subplot(4,2,4);
surf(x_mesh, y_mesh, z_mesh, FT.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
colormap hot;
colorbar;
clim([0 0.3]);
view(2);
title("Transverse Form Factor (from paper), squared  (F_T)^2");

subplot(4,2,6);
surf(x_mesh, y_mesh, z_mesh, F.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
colormap hot;
colorbar;
clim([0 1]);
view(2);
title("Overall Form Factor (from paper), squared  (F)^2");

subplot(4,2,8);
surf(x_mesh, y_mesh, z_mesh, F.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
colormap hot;
colorbar;
clim([0 8 * 10^(-3)]);
view(2);
title("Overall Form Factor (from paper), squared, (F_P)^2, Intensity Scale Zoom In");

sgtitle("Developing Intuition for Overall Form factor (from paper) for wavy cylinder on the top Tangent Plane at \theta = " + num2str(0));

% Form Factor Comparison Plot
figure(WindowState="maximized");

subplot(2,1,1);
surf(x_mesh, y_mesh, z_mesh, F.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
colormap hot;
colorbar;
clim([0 1]);
view(2);
title("Overall Form Factor (from paper), squared  (F)^2");

subplot(2,1,2);
surf(x_mesh, y_mesh, z_mesh, F.^2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
colormap hot;
colorbar;
clim([0 8 * 10^(-3)]);
view(2);
title("Overall Form Factor (from paper), squared, (F_P)^2, Intensity Scale Zoom In");
ylabel('Lc_{\perp}', 'FontWeight','bold');
xlabel('c_{||}', 'FontWeight','bold');

sgtitle("Overall Form factor (from paper) for wavy cylinder on the top Tangent Plane at \theta = " + num2str(0));


% Anger Function Definition
function j = AngerFunc(nu, z)
    j = (1/pi) * integral(@(t) cos(nu .* t - z .* sin(t)), 0, pi);
end

% Custom Cartesian to Spherical mapping
function [theta_rad, phi_rad] = custom_CartToSph(cartsian_mesh)
    x_mesh = cartsian_mesh.x_mesh_val;
    y_mesh = cartsian_mesh.y_mesh_val;
    z_mesh = cartsian_mesh.z_mesh_val;
    grid_size = cartsian_mesh.grid_size_val;

    theta_rad = nan(grid_size, grid_size);
    phi_rad = nan(grid_size, grid_size);

    for i_iter = 1:grid_size
        
        for j_iter = 1:grid_size
            
            x_iter = x_mesh(i_iter, j_iter);
            y_iter = y_mesh(i_iter, j_iter);
            z_iter = z_mesh(i_iter, j_iter);
    
            [azimuth_rad, elevation_rad, ~] = cart2sph(x_iter, y_iter, z_iter);
            theta_rad_iter = pi/2 - elevation_rad;
            if azimuth_rad >= 0
                phi_rad_iter = azimuth_rad;
            else
                phi_rad_iter = 2*pi + azimuth_rad;
            end
    
            theta_rad(i_iter, j_iter) = theta_rad_iter;
            phi_rad(i_iter, j_iter) = phi_rad_iter;

        end

    end

end