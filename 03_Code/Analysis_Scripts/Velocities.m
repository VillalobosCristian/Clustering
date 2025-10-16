%% Particle Count in Circular ROI Over Time
clear all; clc; close all;

%% Parameters
MIN_TRACK_LENGTH = 50;
TIME_STEP = 1;  % Analyze every TIME_STEP frames

%% Load data
[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), return; end

fprintf('Loading data...\n');
data = readtable(fullfile(path_name, file_name));
data = data(4:end,:);

if iscell(data.TRACK_ID), data.TRACK_ID = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X), data.POSITION_X = cellfun(@str2double, data.POSITION_X); end
if iscell(data.POSITION_Y), data.POSITION_Y = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T), data.POSITION_T = cellfun(@str2double, data.POSITION_T); end

sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

%% Extract trajectories
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

%% Get all unique time points
all_times = [];
for i = 1:length(T)
    all_times = [all_times; T{i}];
end
unique_times = unique(all_times);
analysis_times = unique_times(1:TIME_STEP:end);

fprintf('Time points to analyze: %d\n', length(analysis_times));
fprintf('Time range: %.1f to %.1f s\n', min(analysis_times), max(analysis_times));

%% Show final positions and select ROI
final_time = max(analysis_times);
final_positions = [];

for i = 1:length(X)
    [~, timeIdx] = min(abs(T{i} - final_time));
    if ~isempty(timeIdx)
        final_positions = [final_positions; X{i}(timeIdx), Y{i}(timeIdx)];
    end
end

fprintf('\nParticles at final time: %d\n', size(final_positions, 1));

%% Interactive ROI selection
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

figure('Position', [100, 100, 900, 700], 'Color', 'w');
scatter(final_positions(:,1), final_positions(:,2), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.7);
axis equal;
xlabel('$x$ ($\mu$m)', 'FontSize', 14);
ylabel('$y$ ($\mu$m)', 'FontSize', 14);
title(sprintf('Final frame ($t = %.1f$ s) - Draw circular ROI', final_time), 'FontSize', 16);
grid on; box on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

fprintf('\n=== DRAW CIRCULAR ROI ===\n');
fprintf('Draw a circle around the region of interest\n');
fprintf('Double-click inside circle when done\n\n');

h = drawcircle('Color', 'r', 'LineWidth', 2);
wait(h);

roi_center = h.Center;
roi_radius = h.Radius;

fprintf('ROI Center: (%.2f, %.2f) μm\n', roi_center(1), roi_center(2));
fprintf('ROI Radius: %.2f μm\n', roi_radius);

%% Count particles in ROI over time
particle_count = zeros(length(analysis_times), 1);
time_values = analysis_times;

fprintf('\nCounting particles over time...\n');

for t_idx = 1:length(analysis_times)
    current_time = analysis_times(t_idx);
    positions = [];
    
    for i = 1:length(X)
        [~, timeIdx] = min(abs(T{i} - current_time));
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
        end
    end
    
    if ~isempty(positions)
        distances = sqrt((positions(:,1) - roi_center(1)).^2 + ...
                        (positions(:,2) - roi_center(2)).^2);
        particle_count(t_idx) = sum(distances <= roi_radius);
    end
    
    if mod(t_idx, 50) == 0
        fprintf('  Processed %d/%d time points\n', t_idx, length(analysis_times));
    end
end

fprintf('Done!\n');

%% Calculate density in ROI
roi_area = pi * roi_radius^2;
density_in_roi = particle_count / roi_area;

%% Plot results
figure('Position', [150, 150, 1400, 900], 'Color', 'w');

% Subplot 1: Particle count over time
subplot(2, 2, 1);
plot(time_values, particle_count, 'b-', 'LineWidth', 2.5);
xlabel('Time (s)', 'FontSize', 14);
ylabel('Number of Particles in ROI', 'FontSize', 14);
title('Particle Accumulation', 'FontSize', 15);
grid on; box on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Subplot 2: Density in ROI over time
subplot(2, 2, 2);
plot(time_values, density_in_roi, 'r-', 'LineWidth', 2.5);
xlabel('Time (s)', 'FontSize', 14);
ylabel('Particle Density ($\mu$m$^{-2}$)', 'FontSize', 14);
title('Density in ROI', 'FontSize', 15);
grid on; box on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Subplot 3: Rate of change
subplot(2, 2, 3);
dt = diff(time_values);
dN = diff(particle_count);
rate = dN ./ dt;
plot(time_values(1:end-1), rate, 'm-', 'LineWidth', 2);
hold on;
yline(0, 'k--', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 14);
ylabel('Rate (particles/s)', 'FontSize', 14);
title('Accumulation Rate', 'FontSize', 15);
grid on; box on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Subplot 4: Snapshots at different times
subplot(2, 2, 4);
n_snapshots = 4;
snapshot_times = linspace(min(time_values), max(time_values), n_snapshots);
colors_snap = parula(n_snapshots);

hold on;
for snap_idx = 1:n_snapshots
    [~, t_idx] = min(abs(analysis_times - snapshot_times(snap_idx)));
    current_time = analysis_times(t_idx);
    
    positions = [];
    for i = 1:length(X)
        [~, timeIdx] = min(abs(T{i} - current_time));
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
        end
    end
    
    if ~isempty(positions)
        distances = sqrt((positions(:,1) - roi_center(1)).^2 + ...
                        (positions(:,2) - roi_center(2)).^2);
        in_roi = distances <= roi_radius;
        
        scatter(positions(in_roi,1), positions(in_roi,2), 30, ...
                colors_snap(snap_idx,:), 'filled', 'MarkerFaceAlpha', 0.5);
    end
end

viscircles(roi_center, roi_radius, 'Color', 'k', 'LineWidth', 2);
axis equal;
xlabel('$x$ ($\mu$m)', 'FontSize', 14);
ylabel('$y$ ($\mu$m)', 'FontSize', 14);
title('Particle Positions Over Time', 'FontSize', 15);
grid on; box on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

%% Summary statistics
fprintf('\n=== SUMMARY ===\n');
fprintf('ROI area: %.2f μm²\n', roi_area);
fprintf('Initial count: %d particles\n', particle_count(1));
fprintf('Final count: %d particles\n', particle_count(end));
fprintf('Net change: %+d particles\n', particle_count(end) - particle_count(1));
fprintf('Initial density: %.4f particles/μm²\n', density_in_roi(1));
fprintf('Final density: %.4f particles/μm²\n', density_in_roi(end));
fprintf('Density increase: %.1fx\n', density_in_roi(end)/density_in_roi(1));

% Find time of fastest accumulation
[max_rate, max_idx] = max(rate);
fprintf('\nFastest accumulation: %.2f particles/s at t = %.1f s\n', ...
        max_rate, time_values(max_idx));

%% Optional: Save data
% save_file = fullfile(path_name, [file_name(1:end-4), '_roi_count.mat']);
% save(save_file, 'time_values', 'particle_count', 'density_in_roi', ...
%      'roi_center', 'roi_radius', 'roi_area');
% fprintf('Results saved to: %s\n', save_file);