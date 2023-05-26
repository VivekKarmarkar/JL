fileName = "Collagen_stack_1_uncut_X20_z20_processed_1.jpg";

I = imread(fileName);
imshow(I);

nameSplit = strsplit(fileName, ".");
nameSplitFirst = nameSplit(1);
img_idx_str = extract(nameSplitFirst,strlength(nameSplitFirst));
img_idx = str2double(img_idx_str);

if img_idx < 3
    beep;
    disp("grab calibration data");
    calibration_data = ginput(2);
    x = calibration_data(:,1);
    y = calibration_data(:,2);
    calibration_length = 20 * (10^(-6));
    calibration_measurement = diff(x);
    calibration_factor = calibration_length / calibration_measurement;
else
    calibration_factor = 2.8590e-07;
end

beep;
disp("grab waviness data");
[waviness_x, waviness_y] = ginput(20);

beep;
disp("grab rotation data");
rotation_data = ginput(2);
theta = atan((rotation_data(2,2) - rotation_data(1,2))/(rotation_data(2,1) - rotation_data(1,1)));
theta_deg = rad2deg(theta);
R = rotz(-theta_deg);

waviness_x_rot = nan(20,1);
waviness_y_rot = nan(20,1);
for j=1:20
    waviness_iter = [waviness_x(j); waviness_y(j); 0];
    waviness_iter_rot = R*waviness_iter;
    waviness_x_rot(j) = waviness_iter_rot(1);
    waviness_y_rot(j) = waviness_iter_rot(2);
end
figure
plot(waviness_x_rot, waviness_y_rot);
hold on;
f = fit(waviness_x_rot, waviness_y_rot, 'poly9');
plot(f, waviness_x_rot, waviness_y_rot);

domain_pts = linspace(waviness_x_rot(1), waviness_x_rot(end));
func_domain_pts = feval(f, domain_pts);
func_domain_pts_shifted = func_domain_pts - mean(func_domain_pts);
figure
plot(domain_pts, func_domain_pts_shifted, 'r');
hold on;
yline(0, 'k');
disp("grab roots data");
[roots_x, roots_y] = ginput(3);
half_period_uncalibrated = mean(diff(roots_x));
half_period = half_period_uncalibrated * calibration_factor;
waviness_period = half_period * 2;
waviness_period_microns = waviness_period * (10^(6));
disp("waviness period = " + num2str(waviness_period_microns) + " microns");
