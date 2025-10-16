%% SIMPLE RADIAL VELOCITY vs DISTANCE ANALYSIS
% Manual region selection + radial velocity calculation
% Similar to Ramírez-Ramírez et al. analysis

clear all; clc; close all;

%% ========================================================================
%  PARAMETERS
%% ========================================================================
MIN_TRACK_LENGTH = 50;           % Minimum track length to include

%% ========================================================================
%  LOAD DATA
%% ========================================================================
[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), return; end

fprintf('Loading data...\n');
data = readtable(fullfile(path_name, file_name));
data = data(4:end,:);  % Remove TrackMate headers

% Convert to numeric
if iscell(data.TRACK_ID), data.TRACK_ID = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X), data.POSITION_X = cellfun(@str2double, data.POSITION_X); end
if iscell(data.POSITION_Y), data.POSITION_Y = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T), data.POSITION_T = cellfun(@str2double, data.POSITION_T); end

sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

%% ========================================================================
%  EXTRACT TRAJECTORIES
%% ========================================================================
fprintf('Extracting trajectories...\n');
X = {}; Y = {}; T = {};
count = 0;

for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    if height(track_data) >= MIN_TRACK_LENGTH
        count = count + 1;
        X{count} = track_data.POSITION_X;
        Y{count} = track_data.POSITION_Y;
        T{count} = track_data.POSITION_T;
    end
end

fprintf('Valid trajectories: %d\n', count);

%% ========================================================================
%  STEP 1: SHOW ALL TRAJECTORIES - USER SELECTS CENTER
%% ========================================================================
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

figure('Position', [100, 100, 900, 700], 'Color', 'w');
hold on;

% Plot all trajectories in gray
for i = 1:length(X)
    plot(X{i}, Y{i}, '-', 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
end

% Plot final positions
final_x = cellfun(@(x) x(end), X);
final_y = cellfun(@(y) y(end), Y);
scatter(final_x, final_y, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);

xlabel('$x$ ($\mu$m)', 'FontSize', 14);
ylabel('$y$ ($\mu$m)', 'FontSize', 14);
title('Click on the CENTER of the light spot', 'FontSize', 16);
grid on; box on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'FontSize', 12);
axis equal;

fprintf('\n=== CLICK ON CENTER ===\n');
[x_center, y_center] = ginput(1);

% Mark center
scatter(x_center, y_center, 200, 'r', 'x', 'LineWidth', 4);
scatter(x_center, y_center, 100, 'r', 'o', 'LineWidth', 3);

fprintf('Center: (%.2f, %.2f) μm\n', x_center, y_center);

%% ========================================================================
%  STEP 2: SELECT RADIUS OF INTEREST
%% ========================================================================
title('Draw a CIRCLE around the region of interest (double-click to finish)', 'FontSize', 14);

% Let user draw ROI
h = drawcircle('Center', [x_center, y_center], 'Radius', 20);
wait(h);  % Wait for user to finish

roi_center = h.Center;
roi_radius = h.Radius;

fprintf('ROI radius: %.2f μm\n', roi_radius);

%% ========================================================================
%  STEP 3: CALCULATE RADIAL VELOCITIES IN ROI
%% ========================================================================
fprintf('\n=== CALCULATING VELOCITIES ===\n');

all_r = [];           % Distance from center
all_vr = [];          % Radial velocity
all_v_mag = [];       % Total velocity magnitude

for i = 1:length(X)
    % For each segment in trajectory
    for j = 1:length(X{i})-1
        % Midpoint position
        x_mid = (X{i}(j) + X{i}(j+1)) / 2;
        y_mid = (Y{i}(j) + Y{i}(j+1)) / 2;
        
        % Distance from center
        r = sqrt((x_mid - x_center)^2 + (y_mid - y_center)^2);
        
        % Only include if inside ROI
        if r <= roi_radius
            % Calculate velocity
            dt = T{i}(j+1) - T{i}(j);
            vx = (X{i}(j+1) - X{i}(j)) / dt;
            vy = (Y{i}(j+1) - Y{i}(j)) / dt;
            
            % Radial component
            dx = x_mid - x_center;
            dy = y_mid - y_center;
            v_radial = (vx * dx + vy * dy) / r;
            
            % Velocity magnitude
            v_mag = sqrt(vx^2 + vy^2);
            
            % Store
            all_r = [all_r; r];
            all_vr = [all_vr; v_radial];
            all_v_mag = [all_v_mag; v_mag];
        end
    end
end

fprintf('Data points in ROI: %d\n', length(all_r));

%% ========================================================================
%  STEP 4: BIN AND AVERAGE
%% ========================================================================
% Create distance bins
dr = 2;  % bin width in μm
r_bins = 0:dr:roi_radius;
n_bins = length(r_bins) - 1;

r_centers = zeros(n_bins, 1);
vr_mean = zeros(n_bins, 1);
vr_std = zeros(n_bins, 1);
vmag_mean = zeros(n_bins, 1);

for i = 1:n_bins
    idx = all_r >= r_bins(i) & all_r < r_bins(i+1);
    
    if sum(idx) > 5  % At least 5 points
        r_centers(i) = (r_bins(i) + r_bins(i+1)) / 2;
        vr_mean(i) = mean(all_vr(idx));
        vr_std(i) = std(all_vr(idx));
        vmag_mean(i) = mean(all_v_mag(idx));
    else
        r_centers(i) = NaN;
    end
end

% Remove NaN bins
valid = ~isnan(r_centers);
r_centers = r_centers(valid);
vr_mean = vr_mean(valid);
vr_std = vr_std(valid);
vmag_mean = vmag_mean(valid);

%% ========================================================================
%  STEP 5: PLOT RESULTS
%% ========================================================================

%% Figure 1: Main result - v_r vs r (like the paper)
fig1 = figure('Position', [150, 150, 900, 700], 'Color', 'w');

% Scatter plot of all data
scatter(all_r, all_vr, 5, 'b', 'filled', 'MarkerFaceAlpha', 0.1);
hold on;

% Binned averages with error bars
errorbar(r_centers, vr_mean, vr_std, 'ro-', ...
         'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
         'CapSize', 8);

% Reference line
yline(0, 'k--', 'LineWidth', 1.5);

xlabel('Distance from center $r$ ($\mu$m)', 'FontSize', 16);
ylabel('Radial velocity $v_r$ ($\mu$m/s)', 'FontSize', 16);
title('Radial Velocity vs Distance', 'FontSize', 16);
legend('Individual measurements', 'Binned average', 'Location', 'best', 'FontSize', 12);
grid on; box on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'FontSize', 14);

%% Figure 2: Log-linear plot (like paper inset)
fig2 = figure('Position', [200, 200, 800, 600], 'Color', 'w');

% Only plot positive v_r for log scale
positive_idx = vr_mean > 0;
if sum(positive_idx) > 3
    semilogy(r_centers(positive_idx), vr_mean(positive_idx), 'ro-', ...
             'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    xlabel('Distance from center $r$ ($\mu$m)', 'FontSize', 16);
    ylabel('Radial velocity $v_r$ ($\mu$m/s)', 'FontSize', 16);
    title('Log-Linear Plot (for thermophoresis)', 'FontSize', 16);
    grid on; box on;
    set(gca, 'LineWidth', 1.5, 'FontSize', 14);
else
    text(0.5, 0.5, 'Not enough positive v_r values for log plot', ...
         'HorizontalAlignment', 'center', 'FontSize', 14);
    axis off;
end

%% Figure 3: Velocity magnitude vs distance
fig3 = figure('Position', [250, 250, 800, 600], 'Color', 'w');

scatter(all_r, all_v_mag, 5, 'b', 'filled', 'MarkerFaceAlpha', 0.1);
hold on;
plot(r_centers, vmag_mean, 'ro-', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r');

xlabel('Distance from center $r$ ($\mu$m)', 'FontSize', 16);
ylabel('Velocity magnitude $|v|$ ($\mu$m/s)', 'FontSize', 16);
title('Total Velocity vs Distance', 'FontSize', 16);
legend('Individual', 'Binned average', 'Location', 'best', 'FontSize', 12);
grid on; box on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'FontSize', 14);

%% ========================================================================
%  SUMMARY
%% ========================================================================
fprintf('\n=== SUMMARY ===\n');
fprintf('ROI radius: %.2f μm\n', roi_radius);
fprintf('Data points: %d\n', length(all_r));
fprintf('Mean radial velocity: %.3f μm/s\n', mean(all_vr));
fprintf('Mean velocity magnitude: %.3f μm/s\n', mean(all_v_mag));
fprintf('\nInterpretation:\n');
fprintf('  Negative v_r = particles moving TOWARD center\n');
fprintf('  Positive v_r = particles moving AWAY from center\n');
fprintf('  v_r increasing with r = thermophoretic signature\n');