clear all; clc; close all;

%% ========================================================================
%% VELOCITY FIELD VISUALIZATION FROM CSV
%% ========================================================================
%% Generates 2-panel figure: Quiver plot + Radial velocity profile
%% ========================================================================

%% GLOBAL FORMATTING PARAMETERS (CONSISTENT WITH MAIN ANALYSIS)
FIGURE_WIDTH = 12;      % inches
FIGURE_HEIGHT = 5;      % inches

FONT_SIZE_LABELS = 18;  % Axis labels
FONT_SIZE_AXES = 18;    % Tick labels
FONT_SIZE_LEGEND = 12;  % Legend

LINE_WIDTH = 2.5;       % Main plot lines
AXES_WIDTH = 1.5;       % Axes line width

SMOOTH_WINDOW = 5;      % For v_r smoothing

%% SETUP PLOTTING
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

%% LOAD CSV FILE
[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), error('No file selected'); end

[~, experiment_name, ~] = fileparts(file_name);
fprintf('\n=== Loading: %s ===\n', experiment_name);

%% GET EXPERIMENTAL PARAMETERS
prompt = {'Frame rate (fps):', 'Objective (e.g., 20x, 40x, 60x):', 'Power (mW):'};
dlgtitle = 'Experimental Parameters';
defaultans = {'10', '60x', '50'};
answer = inputdlg(prompt, dlgtitle, [1 50], defaultans);

if isempty(answer), error('Parameters not provided'); end

fps = str2double(answer{1});
objective = answer{2};
power_mW = str2double(answer{3});
time_per_frame = 1 / fps;

fprintf('Frame rate: %d fps\n', fps);
fprintf('Objective: %s\n', objective);
fprintf('Power: %.1f mW\n\n', power_mW);

%% DEFINE ILLUMINATION RADIUS BASED ON OBJECTIVE
if strcmp(objective, '60x')
    illumination_radius = 25;  % μm
elseif strcmp(objective, '40x')
    illumination_radius = 30;  % μm
elseif strcmp(objective, '20x')
    illumination_radius = 75;  % μm
else
    illumination_radius = 50;  % μm (default)
end
fprintf('Illumination zone radius: %.0f μm\n\n', illumination_radius);

%% LOAD AND PROCESS CSV DATA
fprintf('Loading CSV data...\n');
data = readtable(fullfile(path_name, file_name));
data = data(4:end,:);  % Skip TrackMate headers

% Convert to numeric if needed
if iscell(data.TRACK_ID)
    data.TRACK_ID = cellfun(@str2double, data.TRACK_ID);
    data.POSITION_X = cellfun(@str2double, data.POSITION_X);
    data.POSITION_Y = cellfun(@str2double, data.POSITION_Y);
    data.POSITION_T = cellfun(@str2double, data.POSITION_T);
end

sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

fprintf('Total tracks: %d\n', length(unique_tracks));

%% EXTRACT TRAJECTORIES
X = {}; Y = {}; T = {};
mean_speeds = [];

% Check if time is in frames or seconds
first_track = sorted_data(sorted_data.TRACK_ID == unique_tracks(1), :);
sample_times = first_track.POSITION_T(1:min(10, height(first_track)));
median_dt = median(diff(sample_times));

if median_dt > 0.5  % Likely frame numbers
    time_is_frames = true;
    fprintf('Time detected as FRAME NUMBERS → converting to seconds\n');
else
    time_is_frames = false;
    fprintf('Time detected as SECONDS → no conversion needed\n');
end

for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    
    if height(track_data) < 60, continue; end
    
    X{end+1} = track_data.POSITION_X;
    Y{end+1} = track_data.POSITION_Y;
    
    % Convert time if needed
    if time_is_frames
        T{end+1} = track_data.POSITION_T * time_per_frame;
    else
        T{end+1} = track_data.POSITION_T;
    end
    
    % Compute mean speed
    x_smooth = movmean(X{end}, 5);
    y_smooth = movmean(Y{end}, 5);
    dx = diff(x_smooth);
    dy = diff(y_smooth);
    dt = diff(T{end});
    speeds = sqrt(dx.^2 + dy.^2) ./ dt;
    valid_speeds = speeds(isfinite(speeds));
    
    if ~isempty(valid_speeds)
        mean_speeds(end+1) = mean(valid_speeds);
    end
end

n_tracks = length(X);
fprintf('Valid trajectories (≥60 points): %d\n\n', n_tracks);

%% COMPUTE VELOCITY GRID
n_grid = 40;
x_min = min(cellfun(@min, X));
x_max = max(cellfun(@max, X));
y_min = min(cellfun(@min, Y));
y_max = max(cellfun(@max, Y));

x_edges = linspace(x_min, x_max, n_grid+1);
y_edges = linspace(y_min, y_max, n_grid+1);
x_grid = (x_edges(1:end-1) + x_edges(2:end))/2;
y_grid = (y_edges(1:end-1) + y_edges(2:end))/2;

vx_grid = cell(n_grid, n_grid);
vy_grid = cell(n_grid, n_grid);

fprintf('Computing velocity field...\n');

for i = 1:length(X)
    x_traj = X{i};
    y_traj = Y{i};
    t_traj = T{i};
    
    % Smooth trajectories
    x_smooth = movmean(x_traj, 5);
    y_smooth = movmean(y_traj, 5);
    
    % Compute velocities
    vx = diff(x_smooth) ./ diff(t_traj);
    vy = diff(y_smooth) ./ diff(t_traj);
    x_mid = (x_smooth(1:end-1) + x_smooth(2:end))/2;
    y_mid = (y_smooth(1:end-1) + y_smooth(2:end))/2;
    
    % Assign to grid cells
    for j = 1:length(vx)
        if ~isfinite(vx(j)) || ~isfinite(vy(j)), continue; end
        
        ix = find(x_mid(j) >= x_edges(1:end-1) & x_mid(j) < x_edges(2:end), 1);
        iy = find(y_mid(j) >= y_edges(1:end-1) & y_mid(j) < y_edges(2:end), 1);
        
        if ~isempty(ix) && ~isempty(iy)
            vx_grid{iy, ix}(end+1) = vx(j);
            vy_grid{iy, ix}(end+1) = vy(j);
        end
    end
end

fprintf('Velocity grid computed.\n\n');

%% SELECT ZONES
fprintf('=== ZONE SELECTION ===\n');

% Use late frame for visualization
max_time = max(cellfun(@max, T));
snapshot_time = max_time * 0.9;

snapshot_x = [];
snapshot_y = [];

for i = 1:length(X)
    [~, idx] = min(abs(T{i} - snapshot_time));
    if idx > 0 && idx <= length(X{i})
        snapshot_x(end+1) = X{i}(idx);
        snapshot_y(end+1) = Y{i}(idx);
    end
end

% Create selection figure
figure('Position', [100, 100, 800, 700], 'Color', 'w');
scatter(snapshot_x, snapshot_y, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
hold on;

% Show illumination zone reference
x_center_guess = mean(snapshot_x);
y_center_guess = mean(snapshot_y);

viscircles([x_center_guess, y_center_guess], illumination_radius, ...
    'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

xlabel('$x$ ($\mu$m)', 'FontSize', 16);
ylabel('$y$ ($\mu$m)', 'FontSize', 16);
title('Draw circle around CRYSTAL zone', 'FontSize', 18);
axis equal tight;
% grid on;

fprintf('Red dashed circle = illumination zone (%.0f μm radius)\n', illumination_radius);
fprintf('Draw GREEN circle around the crystal\n');
fprintf('Then double-click inside it to continue...\n\n');

h = drawcircle('Color', 'g', 'LineWidth', 3, 'FaceAlpha', 0.1);
wait(h);

x_center = h.Center(1);
y_center = h.Center(2);
crystal_radius = h.Radius;

fprintf('Crystal zone - Center: (%.1f, %.1f) μm, Radius: %.1f μm\n', ...
    x_center, y_center, crystal_radius);
fprintf('Illumination zone - Radius: %.0f μm\n\n', illumination_radius);

close(gcf);

%% COMPUTE RADIAL VELOCITY PROFILE
fprintf('Computing radial velocity profile...\n');

% Average velocities in each grid cell
vx_avg_grid = nan(n_grid, n_grid);
vy_avg_grid = nan(n_grid, n_grid);

for ix = 1:n_grid
    for iy = 1:n_grid
        if length(vx_grid{iy, ix}) >= 10
            vx_avg_grid(iy, ix) = mean(vx_grid{iy, ix});
            vy_avg_grid(iy, ix) = mean(vy_grid{iy, ix});
        end
    end
end

[X_mesh, Y_mesh] = meshgrid(x_grid, y_grid);

% Compute radial velocity components
cell_distances = [];
cell_vr = [];

for ix = 1:n_grid
    for iy = 1:n_grid
        if length(vx_grid{iy, ix}) < 5, continue; end
        
        vx_avg = mean(vx_grid{iy, ix});
        vy_avg = mean(vy_grid{iy, ix});
        
        x_cell = x_grid(ix);
        y_cell = y_grid(iy);
        r_cell = sqrt((x_cell - x_center)^2 + (y_cell - y_center)^2);
        
        if r_cell > 1e-6
            rx = (x_cell - x_center) / r_cell;
            ry = (y_cell - y_center) / r_cell;
            vr = vx_avg * rx + vy_avg * ry;
            
            cell_distances(end+1) = r_cell;
            cell_vr(end+1) = vr;
        end
    end
end

% Bin radial velocities
n_bins = 20;
final_x = cellfun(@(x) x(end), X);
final_y = cellfun(@(y) y(end), Y);
distances = sqrt((final_x - x_center).^2 + (final_y - y_center).^2);
r_edges = linspace(0, max(distances), n_bins+1);
r_centers = (r_edges(1:end-1) + r_edges(2:end))/2;

vr_mean = zeros(n_bins, 1);
vr_sem = zeros(n_bins, 1);

for k = 1:n_bins
    idx = cell_distances >= r_edges(k) & cell_distances < r_edges(k+1);
    if sum(idx) > 2
        vr_mean(k) = mean(cell_vr(idx));
        vr_sem(k) = std(cell_vr(idx)) / sqrt(sum(idx));
    else
        vr_mean(k) = NaN;
        vr_sem(k) = NaN;
    end
end

% Apply smoothing
vr_mean_smooth = movmean(vr_mean, SMOOTH_WINDOW, 'omitnan');

% Find minimum
[vr_min, idx_min] = min(vr_mean_smooth);
r_min = r_centers(idx_min);

fprintf('Minimum v_r: %.4f μm/s at r = %.1f μm\n', vr_min, r_min);
fprintf('Ratio r_min/R_illum: %.3f\n\n', r_min/illumination_radius);

%% CREATE OUTPUT DIRECTORY
output_dir = fullfile(path_name, 'velocity_analysis');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ========================================================================
%% GENERATE FIGURE: VELOCITY FIELD + RADIAL PROFILE (2-PANEL)
%% ========================================================================

fprintf('Generating velocity field figure...\n');

figure('Position', [100, 100, 1600, 700], 'Color', 'w');

%% PANEL A: VELOCITY FIELD (QUIVER PLOT)
subplot(1,2,1);
hold on;

% Plot velocity field
quiver(X_mesh, Y_mesh, vx_avg_grid, vy_avg_grid, 2, ...
       'LineWidth', 1.5, 'Color', [0.2, 0.4, 0.8],'ShowArrowHead','off');

% Overlay zones
viscircles([x_center, y_center], crystal_radius, ...
           'Color', [0, 0.7, 0], 'LineWidth', 3);

viscircles([x_center, y_center], illumination_radius, ...
           'Color', [0.8, 0, 0], 'LineWidth', 3, 'LineStyle', '--');

viscircles([x_center, y_center], r_min, ...
           'Color', [0.8, 0, 0.8], 'LineWidth', 2.5, 'LineStyle', ':');

% Mark center
scatter(x_center, y_center, 300, 'k', 'x', 'LineWidth', 4);
scatter(x_center, y_center, 150, 'k', 'o', 'LineWidth', 2.5);

% Add text annotations
text(x_center + crystal_radius*0.7, y_center + crystal_radius*0.7, ...
     'Crystal', 'FontSize', 14, 'Color', [0, 0.7, 0], ...
     'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', [0, 0.7, 0]);

text(x_center + illumination_radius*0.7, y_center - illumination_radius*0.7, ...
     'Illumination', 'FontSize', 14, 'Color', [0.8, 0, 0], ...
     'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', [0.8, 0, 0]);

text(x_center - r_min*0.7, y_center + r_min*0.7, ...
     '$v_r$ min', 'FontSize', 13, 'Color', [0.8, 0, 0.8], ...
     'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', [0.8, 0, 0.8], ...
     'Interpreter', 'latex');

xlabel('$x\ (\mu\mathrm{m})$', 'FontSize', FONT_SIZE_LABELS);
ylabel('$y\ (\mu\mathrm{m})$', 'FontSize', FONT_SIZE_LABELS);

axis equal tight;
% grid on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);
box on;

%% PANEL B: RADIAL VELOCITY PROFILE v_r(r)
subplot(1,2,2);
hold on;

% Main curve with error bars
errorbar(r_centers, vr_mean_smooth, vr_sem, 'o-', ...
         'LineWidth', LINE_WIDTH, 'MarkerSize', 10, ...
         'Color', [0.2, 0.4, 0.8], 'MarkerFaceColor', [0.2, 0.4, 0.8], ...
         'CapSize', 8, 'DisplayName', '$\langle v_r \rangle \pm$ SEM');

% Mark minimum
plot(r_min, vr_min, 'p', 'MarkerSize', 22, ...
     'MarkerFaceColor', [0.8, 0, 0.8], 'MarkerEdgeColor', 'k', ...
     'LineWidth', 2.5, 'DisplayName', sprintf('Min at $r=%.1f$ $\\mu$m', r_min));

% Zone markers
xline(crystal_radius, '-', 'Crystal', 'LineWidth', 3, ...
      'Color', [0, 0.7, 0], 'Alpha', 0.8, ...
      'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center', ...
      'FontSize', 13, 'Interpreter', 'latex', 'HandleVisibility', 'off');

xline(illumination_radius, '--', 'Illumination', 'LineWidth', 3, ...
      'Color', [0.8, 0, 0], 'Alpha', 0.8, ...
      'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', ...
      'FontSize', 13, 'Interpreter', 'latex', 'HandleVisibility', 'off');

yline(0, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', FONT_SIZE_LABELS);
ylabel('$\langle v_r \rangle\ (\mu\mathrm{m/s})$', 'FontSize', FONT_SIZE_LABELS);

% legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Interpreter', 'latex');
% grid on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);
box on;
xlim([0, max(r_centers)]);

%% SAVE FIGURE
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

filename = sprintf('velocity_field_%s.pdf', experiment_name);
export_fig(fullfile(output_dir, filename), '-pdf', '-painters', '-r300', '-transparent');

fprintf('\n========================================\n');
fprintf('FIGURE GENERATED SUCCESSFULLY\n');
fprintf('========================================\n');
fprintf('Output: %s\n', fullfile(output_dir, filename));
fprintf('Size: %.0f x %.0f inches\n', FIGURE_WIDTH, FIGURE_HEIGHT);
fprintf('Font size: %d pt (labels), %d pt (axes)\n', FONT_SIZE_LABELS, FONT_SIZE_AXES);
fprintf('\nDone!\n');