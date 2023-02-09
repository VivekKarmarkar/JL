% Optical parameters
lambda = 450 * 10^(-9);

% Cylinder dimensions
a = 20 * 10^(-9);
L = 300 * 10^(-6);

% Waviness parameters
waviness_amplitude = 3 * 10^(-6);
waviness_period = 5 * 10^(-6);

% Grid
theta_deg = 20;
theta_start = 0;
theta_end = deg2rad(theta_deg);
grid_size = 10000;
theta = linspace(theta_start, theta_end, grid_size);

% Coordinate values
x = sin(theta);
y = zeros(size(x));
z = cos(theta);

% Wavenumber related computations
k = (2 * pi) / lambda;
FA_arg = k * waviness_amplitude * (cos(theta) - 1);

% Computing nu by scaling x
nu = (waviness_period / lambda) * x;

% Sinc and Inverse Sinc Factors
J0 = sinc(nu);
J0_reciprocal = ones(size(J0)) ./ J0;

% Anger Function Factor
FA = nan(size(J0));
for n_iter = 1: length(nu)
    FA(1, n_iter) = AngerFunc(nu(1, n_iter), FA_arg(1, n_iter));
end

% Waviness Form Factor
FW = FA .* J0_reciprocal;

% Straight Cylinder Form Factor
A = sin(theta);
B = zeros(size(A));
C = cos(theta) - 1;
M = sqrt(B.^2 + C.^2);

scalingFactor_transverse = k*a;
scalingFactor_longitudinal = k*L;

F_transverse = besselj(1, scalingFactor_transverse * M) ./ (scalingFactor_transverse * M);
F_transverse(:,1) = 0.5;
F_longitudinal = sinc(scalingFactor_longitudinal * (A / pi));
FS = 2 * (F_transverse .* F_longitudinal);

% Overall Form Factor
F = FS .* FW;

% Plots for visualizing wavy form factor blocks using factorization
figure(WindowState="maximized");

subplot(4,2,1);
plot(x, J0, 'r', 'DisplayName', 'J0');
hold on;
plot(x, FA, 'b', 'DisplayName', 'Anger Funtion');
yline(0, 'k', 'DisplayName', 'X axis');
title("J0 and Anger Function");
grid;
legend;

subplot(4,2,3);
plot(x, J0_reciprocal.^2, 'r');
hold on;
yline(0, 'k');
title("(1 / J0)^2");
grid;

subplot(4,2,5);
plot(x, FA.^2, 'b');
hold on;
yline(0, 'k');
title("(Anger Function)^2");
grid;

subplot(4,2,7);
plot(x, FW.^2, 'r');
hold on;
yline(0, 'k');
xlabel('coordinate parallel to cylinder axis');
title("(Wavy Form Factor)^2");
grid;

subplot(4,2,2);
plot(x, FW.^2, 'r');
hold on;
yline(0, 'k');
title("(Wavy Form Factor)^2");
grid;

subplot(4,2,4);
plot(x, FS.^2, 'r');
hold on;
yline(0, 'k');
title("(Straight Form Factor)^2");
grid;

subplot(4,2,6);
plot(x, F.^2, 'r');
hold on;
yline(0, 'k');
title("(Overall Form Factor)^2");
grid;

subplot(4,2,8);
plot(x, F.^2, 'r');
hold on;
yline(0, 'k');
xlim([0.081 0.28]);
xlabel('coordinate parallel to cylinder axis');
title("(Overall Form Factor)^2 , Zoomed In");
grid;

sgtitle("Developing Intuition for Wavy Form factor at:  \phi = 0 and 0 \leq \theta \leq " + num2str(theta_deg));

% Wavy Cylinder form factor from paper
P = double(int32(L / (2 * waviness_period)));
normalization_factor = waviness_period / L;
FR = (sin(2 * pi * P * nu) ./ sin(pi * nu));
FR(:,1) = 1;
FRN = normalization_factor * FR;
FPL = (FA) .* FRN;
FP = 2 * (F_transverse .* FPL);

% Plots for visualizing wavy form factor from paper
figure(WindowState="maximized");

subplot(4,2,1);
plot(x, FR, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
ylim([min(FR) - 5, max(FR) + 5]);
title("Ratio of Sines Factor F_R");
grid;

subplot(4,2,3);
plot(x, FRN.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
ylim([0, 1.2]);
title("Ratio of Sines Factor (normalized) , squared  (F_{RN})^2");
grid;

subplot(4,2,5);
plot(x, FA.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Anger function , squared  (F_A)^2");
grid;

subplot(4,2,7);
plot(x, FPL.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Longitudinal Form Factor , squared  (F_{PL})^2");
grid;
xlabel('coordinate parallel to cylinder axis');

subplot(4,2,2);
plot(x, FPL.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Longitudinal Form Factor (from paper), squared  (F_{PL})^2");
grid;

subplot(4,2,4);
plot(x, F_transverse.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Transverse Form Factor (from paper), squared  (F_T)^2");
grid;

subplot(4,2,6);
plot(x, FP.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Overall Form Factor (from paper), squared  (F_P)^2");
grid;

subplot(4,2,8);
plot(x, FP.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Overall Form Factor (from paper), squared, Zoomed In  (F_P)^2");
xlim([0.081 0.28]);
grid;
xlabel('coordinate parallel to cylinder axis');

sgtitle("Developing Intuition for Overall Form factor (from paper) at:  \phi = 0 and 0 \leq \theta \leq " + num2str(theta_deg));

% Plots for visualizing the comparison between overall form factors
figure(WindowState="maximized");

subplot(2,1,1);
plot(x, F.^2, 'r', 'DisplayName', 'Factorization based');
hold on;
plot(x, FP.^2, 'b', 'DisplayName', 'Formula from paper');
yline(0, 'k', 'DisplayName', 'X axis');
title("(Overall Form Factor)^2");
grid;
legend;

subplot(2,1,2);
plot(x, F.^2, 'r', 'DisplayName', 'Factorization based');
hold on;
plot(x, FP.^2, 'b', 'DisplayName', 'Formula from paper');
yline(0, 'k', 'DisplayName', 'X axis');
title("(Overall Form Factor)^2 , Zoomed In");
xlim([0.081 0.28]);
grid;
xlabel('coordinate parallel to cylinder axis');

sgtitle("Comparing Overall Form factors at:  \phi = 0 and 0 \leq \theta \leq " + num2str(theta_deg));

% All plots repeated with nu on plot x-axis

% nu Plots for visualizing wavy form factor blocks using factorization
figure(WindowState="maximized");

subplot(4,2,1);
plot(nu, J0, 'r', 'DisplayName', 'J0');
hold on;
plot(nu, FA, 'b', 'DisplayName', 'Anger Funtion');
yline(0, 'k', 'DisplayName', 'X axis');
title("J0 and Anger Function");
grid;
legend;

subplot(4,2,3);
plot(nu, J0_reciprocal.^2, 'r');
hold on;
yline(0, 'k');
title("(1 / J0)^2");
grid;

subplot(4,2,5);
plot(nu, FA.^2, 'b');
hold on;
yline(0, 'k');
title("(Anger Function)^2");
grid;

subplot(4,2,7);
plot(nu, FW.^2, 'r');
hold on;
yline(0, 'k');
xlabel('nu');
title("(Wavy Form Factor)^2");
grid;

subplot(4,2,2);
plot(nu, FW.^2, 'r');
hold on;
yline(0, 'k');
title("(Wavy Form Factor)^2");
grid;

subplot(4,2,4);
plot(nu, FS.^2, 'r');
hold on;
yline(0, 'k');
title("(Straight Form Factor)^2");
grid;

subplot(4,2,6);
plot(nu, F.^2, 'r');
hold on;
yline(0, 'k');
title("(Overall Form Factor)^2");
grid;

subplot(4,2,8);
plot(nu, F.^2, 'r');
hold on;
yline(0, 'k');
xlim([0.9 3.11]);
xlabel('nu');
title("(Overall Form Factor)^2 , Zoomed In");
grid;

sgtitle("nu Plots: Developing Intuition for Wavy Form factor at:  \phi = 0 and 0 \leq \theta \leq " + num2str(theta_deg));

% nu Plots for visualizing wavy form factor from paper
figure(WindowState="maximized");

subplot(4,2,1);
plot(nu, FR, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
ylim([min(FR) - 5, max(FR) + 5]);
title("Ratio of Sines Factor F_R");
grid;

subplot(4,2,3);
plot(nu, FRN.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
ylim([0, 1.2]);
title("Ratio of Sines Factor (normalized) , squared  (F_{RN})^2");
grid;

subplot(4,2,5);
plot(nu, FA.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Anger function , squared  (F_A)^2");
grid;

subplot(4,2,7);
plot(nu, FPL.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Longitudinal Form Factor , squared  (F_{PL})^2");
grid;
xlabel('nu');

subplot(4,2,2);
plot(nu, FPL.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Longitudinal Form Factor (from paper), squared  (F_{PL})^2");
grid;

subplot(4,2,4);
plot(nu, F_transverse.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Transverse Form Factor (from paper), squared  (F_T)^2");
grid;

subplot(4,2,6);
plot(nu, FP.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Overall Form Factor (from paper), squared  (F_P)^2");
grid;

subplot(4,2,8);
plot(nu, FP.^2, 'r');
hold on;
yline(0, 'k', 'DisplayName', 'X axis');
title("Overall Form Factor (from paper), squared, Zoomed In  (F_P)^2");
xlim([0.9 3.11]);
grid;
xlabel('nu');

sgtitle("nu Plots: Developing Intuition for Overall Form factor (from paper) at:  \phi = 0 and 0 \leq \theta \leq " + num2str(theta_deg));

% Plots for visualizing the comparison between overall form factors
figure(WindowState="maximized");

subplot(2,1,1);
plot(nu, F.^2, 'r', 'DisplayName', 'Factorization based');
hold on;
plot(nu, FP.^2, 'b', 'DisplayName', 'Formula from paper');
yline(0, 'k', 'DisplayName', 'X axis');
title("(Overall Form Factor)^2");
grid;
legend;

subplot(2,1,2);
plot(nu, F.^2, 'r', 'DisplayName', 'Factorization based');
hold on;
plot(nu, FP.^2, 'b', 'DisplayName', 'Formula from paper');
yline(0, 'k', 'DisplayName', 'X axis');
title("(Overall Form Factor)^2 , Zoomed In");
xlim([0.9 3.11]);
grid;
xlabel('nu');

sgtitle("nu Plots: Comparing Overall Form factors at:  \phi = 0 and 0 \leq \theta \leq " + num2str(theta_deg));

% Anger Function Definition
function j = AngerFunc(nu, z)
    j = (1/pi) * integral(@(t) cos(nu .* t - z .* sin(t)), 0, pi);
end