% Clear workspace and command window
clear all;
clc;
close all

% Open file dialog to zselect CSV file
[file_name, path_name] = uigetfile('*.csv', 'Select a CSV file');

% Check if user didn't cancel the file dialog
% if isequal(file_name, 0)
    disp('User selected Cancel');
% else
    disp(['User selected ', fullfile(path_name, file_name)]);
    
    % Read the selected CSV file
    data = readtable(fullfile(path_name, file_name));
    
    % Remove non-data rows (assuming first 3 rows are headers or metadata)
    data = data(4:end,:);
    
    % Convert columns to numeric values if they are cell arrays
    if iscell(data.TRACK_ID)
        data.TRACK_ID = cellfun(@str2double, data.TRACK_ID);
    end
    
    if iscell(data.POSITION_X)
        data.POSITION_X = cellfun(@str2double, data.POSITION_X);
    end
    
    if iscell(data.POSITION_Y)
        data.POSITION_Y = cellfun(@str2double, data.POSITION_Y);
    end
    
    if iscell(data.POSITION_T)
        data.POSITION_T = cellfun(@str2double, data.POSITION_T);
    end
    
    % Sort data by TRACK_ID and POSITION_T
    sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
    
    % Get unique track IDs
    unique_tracks = unique(sorted_data.TRACK_ID);
    num_tracks = length(unique_tracks);
    
    % Initialize storage for trajectory data
    X = cell(num_tracks, 1);
    Y = cell(num_tracks, 1);
    T = cell(num_tracks, 1);
    mean_speeds = zeros(num_tracks, 1);
    
    % Process each trajectory
    valid_tracks = 0;
    for i = 1:length(unique_tracks)
        % Get data for this specific track
        track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
        
        if height(track_data) < 100
            continue;
        end
        
        valid_tracks = valid_tracks + 1;
        
        % Store trajectory data
        X{i} = track_data.POSITION_X;
        Y{i} = track_data.POSITION_Y;
        T{i} = track_data.POSITION_T;
        
        % Calculate speeds for this trajectory
        dx = diff(X{i});
        dy = diff(Y{i});
        dt = diff(T{i});
        
        % Calculate speeds for each segment
        speeds = sqrt((dx./dt).^2 + (dy./dt).^2);
        
        % Calculate mean speed (removing any NaN or Inf values)
        valid_speeds = speeds(isfinite(speeds));
        if ~isempty(valid_speeds)
            mean_speeds(i) = mean(valid_speeds);
        else
            % If no valid speeds, set to NaN
            mean_speeds(i) = NaN;
        end
    end
    
    fprintf('Processing %d valid trajectories.\n', valid_tracks);
    
    % Create a mask for valid trajectories
    valid_mask = ~cellfun('isempty', X);
    
    % Remove empty cells (trajectories that were filtered out)
    X = X(valid_mask);
    Y = Y(valid_mask);
    T = T(valid_mask);
    mean_speeds = mean_speeds(valid_mask);
    
    % Remove any NaN speeds
    valid_speed_mask = isfinite(mean_speeds);
    X = X(valid_speed_mask);
    Y = Y(valid_speed_mask);
    T = T(valid_speed_mask);
    mean_speeds = mean_speeds(valid_speed_mask);
    
    % Check if we have any valid trajectories
    if isempty(X) || isempty(Y)
        error('No valid trajectories found after filtering.');
    end
    
    % Get min and max speeds for color scaling
    min_speed = min(mean_speeds);
    max_speed = max(mean_speeds);
    
%%
figure('Position', [100, 100, 900, 700], 'Color', 'w');
    set(gcf, 'PaperPositionMode', 'auto');
    set(0, 'defaultTextInterpreter', 'latex');
    set(0, 'defaultAxesTickLabelInterpreter', 'latex');
    set(0, 'defaultLegendInterpreter', 'latex');
    set(0, 'defaultColorbarTickLabelInterpreter', 'latex');
    
    hold on;
    % title('$\mathbf{All\ Detected\ Particle\ Trajectories}$', 'FontSize', 16, 'Interpreter', 'latex');
    
    cmap = thermal(256);
    
    % Plot all trajectories
    for i = 1:length(X)
        % Normalize speed for color mapping
        norm_speed = (mean_speeds(i) - min_speed) / (max_speed - min_speed);
        
        % Map to colormap index
        color_idx = max(1, min(256, round(norm_speed * 255) + 1));
        trajectory_color = cmap(color_idx, :);
        
        plot(X{i}, Y{i}, '-', 'Color', trajectory_color, 'LineWidth', 1);
    end
    

    xlabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');

    grid on;
    set(gca, 'GridLineStyle', ':','xtick',[],'ytick',[]);
    set(gca, 'GridAlpha', 0.3);
    
    % Calculate and set axis limits
    x_min = min(cellfun(@min, X));
    x_max = max(cellfun(@max, X));
    y_min = min(cellfun(@min, Y));
    y_max = max(cellfun(@max, Y));
    
    % Add margin
    margin = 0.05 * max(x_max-x_min, y_max-y_min);
    axis([x_min-margin, x_max+margin, y_min-margin, y_max+margin]);
    
    % Improve axes appearance
    set(gca, 'LineWidth', 1.5);
    set(gca, 'TickDir', 'out');
    set(gca, 'TickLength', [0.02 0.02]);
    set(gca, 'FontSize', 12);
    
    % colormap(cmap);
    % cb = colorbar;
    % % cb.Label.String = '$Mean\ Speed\ (\mu m/s)$';
    % cb.Label.FontSize = 14;
    % cb.Label.Interpreter = 'latex';
    % cb.LineWidth = 1.5;
    % cb.TickDirection = 'out';
    % cb.Box = 'off';
    % cb.ticks=[]
    % cb.Ticks = [0, 0.25, 0.5, 0.75, 1];
    % cb.TickLabels = {['$' sprintf('%.2f', min_speed) '\ \mathrm{(Slow)}$'], '', ...
    %                  ['$' sprintf('%.2f', (min_speed+max_speed)/2) '$'], '', ...
    %                  ['$' sprintf('%.2f', max_speed) '\ \mathrm{(Fast)}$']};
    
    box on;
    set(gca, 'LineWidth', 1.5);
    
    set(gcf, 'Units', 'inches');
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1), pos(2), 8, 6]);
 
    baseDir = pwd;
fileName = 'Trajectories_fin';
% savefigures(baseDir, fileName);
% end

%%
figure('Position', [100, 100, 900, 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultColorbarTickLabelInterpreter', 'latex');

hold on;

for i = 1:length(X)
    % Normalize speed for color mapping
    norm_speed = (mean_speeds(i) - min_speed) / (max_speed - min_speed);
        color_idx = max(1, min(256, round(norm_speed * 255) + 1));
    marker_color = cmap(color_idx, :);
    
    scatter(X{i}(end), Y{i}(end), 40, marker_color, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
end

xlabel('$X (\mu m)$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$Y (\mu m)$', 'FontSize', 14, 'Interpreter', 'latex');
grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);

axis([x_min-margin, x_max+margin, y_min-margin, y_max+margin]);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02 0.02]);
set(gca, 'FontSize', 12);

colormap(cmap);
cb = colorbar;
cb.Label.String = '$\mathbf{Mean\ Speed\ (\mu m/s)}$';
cb.Label.FontSize = 14;
cb.Label.Interpreter = 'latex';
cb.LineWidth = 1.5;
cb.TickDirection = 'out';
cb.Box = 'off';
cb.Ticks = [0, 0.25, 0.5, 0.75, 1];
cb.TickLabels = {['$' sprintf('%.2f', min_speed) '\ \mathrm{(Slow)}$'], '', ...
                 ['$' sprintf('%.2f', (min_speed+max_speed)/2) '$'], '', ...
                 ['$' sprintf('%.2f', max_speed) '\ \mathrm{(Fast)}$']};

box on;
set(gca, 'LineWidth', 1.5);
    xlabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');

    grid on;
    set(gca, 'GridLineStyle', ':','xtick',[],'ytick',[]);
    set(gca, 'GridAlpha', 0.3);
set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 8, 6]);
    baseDir = pwd;
fileName = 'Trajectories_fin';
savefigures(baseDir, fileName);

%% 
final_positions_x = zeros(length(X), 1);
final_positions_y = zeros(length(X), 1);

for i = 1:length(X)
    final_positions_x(i) = X{i}(end);
    final_positions_y(i) = Y{i}(end);
end

%
figure('Position', [100, 100, 900, 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultColorbarTickLabelInterpreter', 'latex');

nbins = 500; % 
h = histogram2(final_positions_x, final_positions_y, nbins, ...
               'DisplayStyle', 'tile', ...
               'ShowEmptyBins', 'on', ...
               'XBinLimits', [x_min-margin, x_max+margin], ...
               'YBinLimits', [y_min-margin, y_max+margin], ...
               'Normalization', 'count');
x_edges = h.XBinEdges;
y_edges = h.YBinEdges;
counts = h.Values;

counts_smooth = imgaussfilt(counts, 1.5); %
clf;
[X_mesh, Y_mesh] = meshgrid(x_edges, y_edges);

pcolor(X_mesh, Y_mesh, [counts_smooth, zeros(size(counts_smooth,1),1); zeros(1,size(counts_smooth,2)+1)]);

colormap(magma(256));

shading interp;

cb = colorbar;
cb.Label.String = '$\mathbf{Particle\ Density}$';
cb.Label.FontSize = 14;
cb.Label.Interpreter = 'latex';
cb.LineWidth = 1.5;
cb.TickDirection = 'out';
cb.Box = 'off';

    xlabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');

axis([x_min-margin, x_max+margin, y_min-margin, y_max+margin]);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02 0.02]);
set(gca, 'FontSize', 12);

grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);

box on;


set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 8, 6]);
%     baseDir = pwd;
% fileName = 'Trajectories_fin';
% savefigures(baseDir, fileName);
%%
%% SPATIAL VELOCITY
% Create spatial bins
nx = 20; ny = 20;
x_bins = linspace(x_min, x_max, nx+1);
y_bins = linspace(y_min, y_max, ny+1);

% Initialize velocity arrays
vx_map = zeros(ny, nx);
vy_map = zeros(ny, nx);

count_map = zeros(ny, nx);

% Accumulate velocities in each spatial bin
for i = 1:length(X)
    % For each point in trajectory except the last
    for j = 1:length(X{i})-1
        % Calculate velocity
        dt = T{i}(j+1) - T{i}(j);
        vx = (X{i}(j+1) - X{i}(j))/dt;
        vy = (Y{i}(j+1) - Y{i}(j))/dt;
        
        % Find spatial bin
        x_idx = find(x_bins <= X{i}(j), 1, 'last');
        y_idx = find(y_bins <= Y{i}(j), 1, 'last');
        
        if ~isempty(x_idx) && ~isempty(y_idx) && x_idx <= nx && y_idx <= ny
            vx_map(y_idx, x_idx) = vx_map(y_idx, x_idx) + vx;
            vy_map(y_idx, x_idx) = vy_map(y_idx, x_idx) + vy;
            count_map(y_idx, x_idx) = count_map(y_idx, x_idx) + 1;
        end
    end
end

vx_map = vx_map ./ max(count_map, 1);
vy_map = vy_map ./ max(count_map, 1);

[X_grid, Y_grid] = meshgrid((x_bins(1:end-1) + x_bins(2:end))/2, (y_bins(1:end-1) + y_bins(2:end))/2);

% Calculate velocity magnitude map
v_magnitude = sqrt(vx_map.^2 + vy_map.^2);

final_positions_x = zeros(length(X), 1);
final_positions_y = zeros(length(X), 1);
for i = 1:length(X)
    final_positions_x(i) = X{i}(end);
    final_positions_y(i) = Y{i}(end);
end

% Create publication-quality figure
figure('Position', [100, 100, 900, 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultColorbarTickLabelInterpreter', 'latex');

% Mask bins with no data
mask = count_map > 0;
X_masked = X_grid(mask);
Y_masked = Y_grid(mask);
vx_masked = vx_map(mask);
vy_masked = vy_map(mask);

contourf(X_grid, Y_grid, v_magnitude, 20, 'LineStyle', 'none', 'FaceAlpha', 0.7);

colormap(hot(256));

cb = colorbar;
cb.Label.String = '$\mathbf{Velocity\ Magnitude\ (\mu m/s)}$';
cb.Label.FontSize = 14;
cb.Label.Interpreter = 'latex';
cb.LineWidth = 1.5;
cb.TickDirection = 'out';
cb.Box = 'off';

hold on;

q = quiver(X_masked, Y_masked, vx_masked, vy_masked, 1.5, 'LineWidth', 1.5, 'Color', 'k', 'MaxHeadSize', 0.5);
q.ShowArrowHead = 'on';

s = scatter(final_positions_x, final_positions_y, 10, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
s.MarkerFaceAlpha = 0.3;  % 50% transparency for marker fill
s.MarkerEdgeAlpha = 0.3;  % 70% transparency for marker edge


    xlabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');

axis([x_min-margin, x_max+margin, y_min-margin, y_max+margin]);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02 0.02]);
set(gca, 'FontSize', 16);

% Add grid but keep it subtle
grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);

% Add box with proper line width
box on;
set(gca, 'LineWidth', 1.5);

% % Add legend to explain the visual elements
% scatter_handle = scatter(NaN, NaN, 30, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
% quiver_handle = quiver(NaN, NaN, 1, 0, 'LineWidth', 1.75, 'Color', 'k', 'MaxHeadSize', 0.5);
% legend([quiver_handle, scatter_handle], {'$\mathbf{Velocity\ Vectors}$', '$\mathbf{Final\ Positions}$'}, ...
%        'Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');

set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 8, 6]);
baseDir = pwd;
fileName = 'final_pos_quiver';
savefigures(baseDir, fileName);
%%
%% Simple Velocity Distribution Histogram
% Remove invalid speeds
valid_speeds = mean_speeds(isfinite(mean_speeds) & mean_speeds > 0);

% Create figure
figure('Position', [100, 100, 900, 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

% Create histogram with proper PDF normalization
num_bins = round(sqrt(length(valid_speeds)));
[counts, bin_edges] = histcounts(valid_speeds, num_bins);
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
bin_width = bin_edges(2) - bin_edges(1);
pdf_values = counts / (length(valid_speeds) * bin_width);

% Plot histogram
bar(bin_centers, pdf_values, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'k');

% Labels
xlabel('$\mathrm{Velocity\ (\mu m/s)}$', 'FontSize', 14);
ylabel('$\mathrm{Probability\ Density}$', 'FontSize', 14);

grid on;
box on;
grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02 0.02]);
set(gca, 'FontSize', 16);
box on;

% Set figure size
set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 8, 6]);
baseDir = pwd;
fileName = 'vel_histo';
savefigures(baseDir, fileName);
%%
%% INTERACTIVE CENTER SELECTION AND RADIAL ANALYSIS

% Step 1: Show particle positions and let user click center
figure('Position', [100, 100, 900, 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

% Plot all final positions
scatter(final_positions_x, final_positions_y, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.7);

xlabel('$\mathrm{X\ (\mu m)}$', 'FontSize', 14);
ylabel('$\mathrm{Y\ (\mu m)}$', 'FontSize', 14);
title('$\mathrm{Click\ on\ the\ Center\ of\ the\ Cluster}$', 'FontSize', 16);

grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02 0.02]);
set(gca, 'FontSize', 16);
box on;
fileName = 'final_pos';
axis([x_min-margin, x_max+margin, y_min-margin, y_max+margin]);
set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 8, 6]);

fprintf('Click on the center of the cluster...\n');
[x_center, y_center] = ginput(1);

% Mark the selected center
hold on;
scatter(x_center, y_center, 200, 'r', 'x', 'LineWidth', 4);
scatter(x_center, y_center, 100, 'r', 'o', 'LineWidth', 3);

fprintf('Center selected at: (%.1f, %.1f) μm\n', x_center, y_center);

% Step 2: Calculate distances from selected center
distances_from_center = sqrt((final_positions_x - x_center).^2 + (final_positions_y - y_center).^2);

% Clean data (remove invalid speeds)
valid_idx = isfinite(mean_speeds) & mean_speeds > 0;
clean_distances = distances_from_center(valid_idx);
clean_speeds = mean_speeds(valid_idx);
baseDir = pwd;
savefigures(baseDir, fileName);
% Get user click
% Step 3: Plot velocity vs distance
figure('Position', [100, 100, 900, 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

scatter(clean_distances, clean_speeds, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);

xlabel('$\mathrm{Distance\ from\ Selected\ Center\ (\mu m)}$', 'FontSize', 14);
ylabel('$\mathrm{Mean\ Velocity\ (\mu m/s)}$', 'FontSize', 14);

grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02 0.02]);
set(gca, 'FontSize', 16);
box on;

% Set figure size
set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 8, 6]);
baseDir = pwd;
fileName = 'vel_radial';
savefigures(baseDir, fileName);