clear all; clc; close all;

%% Load CSV data 
[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), return; end

data = readtable(fullfile(path_name, file_name));
data = data(4:end,:); % Remove TrackMate headers

% Convert to numeric if needed
if iscell(data.TRACK_ID), data.TRACK_ID = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X), data.POSITION_X = cellfun(@str2double, data.POSITION_X); end
if iscell(data.POSITION_Y), data.POSITION_Y = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T), data.POSITION_T = cellfun(@str2double, data.POSITION_T); end

% Sort data for easier processing
sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

%% Extract trajectories
min_track_length = 100;
X = {}; Y = {}; T = {};
count = 0;

for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    if height(track_data) >= min_track_length
        count = count + 1;
        X{count} = track_data.POSITION_X;
        Y{count} = track_data.POSITION_Y;
        T{count} = track_data.POSITION_T;
    end
end

%% Set up time points for animation
all_times = [];
for i = 1:length(T)
    all_times = [all_times; T{i}];
end
unique_times = unique(all_times);

time_step = 5; 
analysis_times = unique_times(1:time_step:end);

%% Parameters
calib = 1.0;
distance_threshold = 5;

all_x = []; all_y = [];
for i = 1:length(X)
    all_x = [all_x; X{i}];
    all_y = [all_y; Y{i}];
end
all_x = all_x * calib;
all_y = all_y * calib;

x_padding = 0.05 * (max(all_x) - min(all_x));
y_padding = 0.05 * (max(all_y) - min(all_y));
x_limits = [min(all_x) - x_padding, max(all_x) + x_padding];
y_limits = [min(all_y) - y_padding, max(all_y) + y_padding];

%% g(r) calculation parameters
r_max = 30; % Maximum distance to analyze (μm)
dr = 0.1;   % Distance bin width (μm)
r_bins = 0:dr:r_max;
r_centers = r_bins(1:end-1) + dr/2;

% Storage for g(r) evolution
g_r_evolution = [];
clustering_metrics = [];

%% Create video writer
video_name = 'particle_clustering_with_gr.mp4';
v = VideoWriter(video_name, 'MPEG-4');
v.FrameRate = 8; % slower to see g(r) evolution
open(v);

%%
figure('Position', [100, 100, 1800, 600]); % Wider for 3 panels

% Animation loop
for t_idx = 1:length(analysis_times)
    current_time = analysis_times(t_idx);
    
    % Extract positions at this time point
    positions = [];
    particle_ids = [];
    
    for i = 1:length(X)
        timeIdx = find(T{i} == current_time, 1);
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
            particle_ids = [particle_ids; i];
        end
    end
    
    positions = positions * calib;
    N_particles = size(positions, 1);
    
    if N_particles < 3
        continue;
    end
    
    % Create Delaunay triangulation
    DT = delaunayTriangulation(positions(:,1), positions(:,2));
    triangles = DT.ConnectivityList;
    n_triangles = size(triangles, 1);
    
    valid_triangles = [];
    for i = 1:n_triangles
        p1 = triangles(i, 1);
        p2 = triangles(i, 2);
        p3 = triangles(i, 3);
        
        % Calculate the three side lengths
        side1 = norm(positions(p1,:) - positions(p2,:));
        side2 = norm(positions(p2,:) - positions(p3,:));
        side3 = norm(positions(p3,:) - positions(p1,:));
        
        max_side = max([side1, side2, side3]);
        
        % Keep triangle only if sides are below threshold
        if max_side <= distance_threshold
            valid_triangles = [valid_triangles; triangles(i,:)];
        end
    end
    
    n_valid_triangles = size(valid_triangles, 1);
    
    %% Calculate g(r) for current time point
    [g_r, r_centers_actual] = calculate_gr(positions, r_bins, x_limits, y_limits);
    g_r_evolution = [g_r_evolution; g_r'];
    
    % Store clustering metrics
    clustering_ratio = n_valid_triangles / N_particles;
    clustering_metrics = [clustering_metrics; clustering_ratio];
    
    clf;
    
    % Panel 1: Full Delaunay triangulation
    subplot(1, 3, 1);
    triplot(DT, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5);
    hold on;
    scatter(positions(:,1), positions(:,2), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.8);
    axis equal; grid on;
    xlim(x_limits); ylim(y_limits);
    title('Full Delaunay Triangulation');
    xlabel('x (μm)'); ylabel('y (μm)');
    
    % Panel 2: Distance-filtered triangles
    subplot(1, 3, 2);
    scatter(positions(:,1), positions(:,2), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.8);
    hold on;
    
    % Draw only the valid triangles
    for i = 1:n_valid_triangles
        p1 = valid_triangles(i, 1);
        p2 = valid_triangles(i, 2);
        p3 = valid_triangles(i, 3);
        
        % Draw the triangle edges
        plot([positions(p1,1), positions(p2,1)], [positions(p1,2), positions(p2,2)], 'r-', 'LineWidth', 1);
        plot([positions(p2,1), positions(p3,1)], [positions(p2,2), positions(p3,2)], 'r-', 'LineWidth', 1);
        plot([positions(p3,1), positions(p1,1)], [positions(p3,2), positions(p1,2)], 'r-', 'LineWidth', 1);
    end
    
    axis equal; grid on;
    xlim(x_limits); ylim(y_limits);
    title('Distance-Filtered Triangles');
    xlabel('x (μm)'); ylabel('y (μm)');
    
    % Panel 3: g(r) function
    subplot(1, 3, 3);
    plot(r_centers_actual, g_r, 'k-', 'LineWidth', 2);
    hold on;
    plot([0, r_max], [1, 1], 'r--', 'LineWidth', 1); % Random distribution line
    
    % Highlight clustering region
    clustering_region = r_centers_actual <= distance_threshold;
    if any(clustering_region)
        area(r_centers_actual(clustering_region), g_r(clustering_region), 'FaceColor', 'red', 'FaceAlpha', 0.2);
    end
    
    xlim([0, r_max]); ylim([0, max(3, max(g_r)*1.1)]);
    xlabel('Distance r (μm)'); ylabel('g(r)');
    title('Radial Distribution Function');
    grid on;
    
    % Add text annotations
    text(0.7*r_max, 0.9*max(3, max(g_r)*1.1), sprintf('Time = %d', current_time), 'FontSize', 10);
    text(0.7*r_max, 0.8*max(3, max(g_r)*1.1), sprintf('Particles = %d', N_particles), 'FontSize', 10);
    text(0.7*r_max, 0.7*max(3, max(g_r)*1.1), sprintf('Clustering = %.2f', clustering_ratio), 'FontSize', 10);
    
    % Overall title
    sgtitle(sprintf('Particle Clustering Analysis - Time = %d | Particles = %d | Valid/Total = %d/%d', ...
        current_time, N_particles, n_valid_triangles, n_triangles));
    
    drawnow;
    
    % Capture frame and write to video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close video file
close(v);

%% Plot g(r) evolution over time
figure('Position', [100, 100, 1200, 400]);

subplot(1, 2, 1);
imagesc(analysis_times(1:size(g_r_evolution,1)), r_centers_actual, g_r_evolution');
colorbar;
xlabel('Time'); ylabel('Distance r (μm)');
title('g(r) Evolution Over Time');
caxis([0, 3]); % Adjust color scale as needed

subplot(1, 2, 2);
plot(analysis_times(1:length(clustering_metrics)), clustering_metrics, 'b-', 'LineWidth', 2);
xlabel('Time'); ylabel('Clustering Ratio');
title('Clustering Metric Evolution');
grid on;

fprintf('Video saved as: %s\n', video_name);

%% Function to calculate g(r)
function [g_r, r_centers] = calculate_gr(positions, r_bins, x_limits, y_limits)
    N = size(positions, 1);
    
    % Calculate all pairwise distances
    distances = [];
    for i = 1:N
        for j = i+1:N
            dist = norm(positions(i,:) - positions(j,:));
            distances = [distances; dist];
        end
    end
    
    % Create histogram of distances
    [counts, edges] = histcounts(distances, r_bins);
    r_centers = edges(1:end-1) + diff(edges)/2;
    
    % Calculate system density
    area = (x_limits(2) - x_limits(1)) * (y_limits(2) - y_limits(1));
    rho = N / area;
    
    % Normalize by expected counts in random distribution
    g_r = zeros(size(r_centers));
    for i = 1:length(r_centers)
        r = r_centers(i);
        dr = edges(i+1) - edges(i);
        
        % Expected number of particles in annular ring
        ring_area = pi * ((r + dr/2)^2 - (r - dr/2)^2);
        expected_pairs = N * (N-1) / 2 * ring_area / area;
        
        % Avoid division by zero
        if expected_pairs > 0
            g_r(i) = counts(i) / expected_pairs;
        else
            g_r(i) = 0;
        end
    end
end