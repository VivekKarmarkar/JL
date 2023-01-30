% Generate meshgrid of theta and phi values
theta = linspace(0, pi, 100);
phi = linspace(0, 2 * pi, 100);
[theta_mesh, phi_mesh] = meshgrid(theta, phi);

% Physical parameter input
lambda = 450 * 10^(-9);
L = 300 * 10^(-9);
a = 20 * 10^(-9);

% Physical parameter computation
k = (2*pi)/lambda;

% Evaluate the heat function at each (theta, phi) point
h = heat_func(theta_mesh, phi_mesh, k, L, a);

% Convert spherical coordinates to Cartesian coordinates
x = sin(theta_mesh) .* cos(phi_mesh);
y = sin(theta_mesh) .* sin(phi_mesh);
z = cos(theta_mesh);

% Plot the heat map
figure;
%title("Scattering from a cylinder");
surf(x, y, z, h, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
colormap hot;
colorbar;
view(3);
axis equal

% Plot the axes
line([0, 1.5], [0, 0], [0, 0], 'Color', 'red', 'LineWidth', 2);
line([0, 0], [0, 1.5], [0, 0], 'Color', 'green', 'LineWidth', 2);
line([0, 0], [0, 0], [0, 1.5], 'Color', 'blue', 'LineWidth', 2);

% Plotting black arrows parallel to the z-axis
num_arrows = 6;
arrow_length = 0.35;
arrow_base = repmat([0, 0, -1.5], num_arrows, 1);
arrow_base(:,1) = linspace(-0.2, 0.2, num_arrows)';
arrow_direction = repmat([0, 0, arrow_length], num_arrows, 1);
quiver3(arrow_base(:,1), arrow_base(:,2), arrow_base(:,3), arrow_direction(:,1), arrow_direction(:,2), arrow_direction(:,3), 0, 'Color', 'black', 'LineWidth', 0.5);

% Plot the cylinder
[X, Y, Z] = cylinder(0.1, 100);
cyl_pltX = 0.4*(Z-0.5);
cyl_pltY = X;
cyl_pltZ = Y;
c = surf(cyl_pltX, cyl_pltY, cyl_pltZ, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 1);
hold off;

% Define the heat function
function h = heat_func(theta, phi, k, L, a)
  A = sin(theta) .* cos(phi);
  B = sin(theta) .* sin(phi);
  C = cos(theta) - 1;
  M = sqrt(B^2 + C^2);
  F1 = sinc(k*L*A);
  F2 = besselj(1, k*a*M);
  F = 2 * F1 * F2;
  h = F;
end