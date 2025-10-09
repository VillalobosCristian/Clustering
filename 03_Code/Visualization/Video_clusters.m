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
min_track_length = 10;
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

%% Set up video writer
video_filename = [file_name(1:end-4), '_delaunay_animation.mp4'];
video_path = fullfile(path_name, video_filename);
v = VideoWriter(video_path, 'MPEG-4');
v.FrameRate = 10; % Adjust frame rate as needed
v.Quality = 95;
open(v);

%% Create figure with improved aesthetics
fig = figure('Position', [100, 100, 1400, 600], 'Color', 'white');
set(fig, 'Renderer', 'painters'); % Better quality rendering

% Define aesthetic colors
particle_color = [0.2, 0.4, 0.8]; % Nice blue
delaunay_color = [0.8, 0.8, 0.8]; % Light gray
valid_triangle_color = [0.9, 0.2, 0.2]; % Red
background_color = 'white';

fprintf('Creating video: %s\n', video_filename);
fprintf('Processing %d time points...\n', length(analysis_times));

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
    
    % Skip if too few particles
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
    clf;
    
    % Left subplot: Full Delaunay Triangulation
    subplot(1, 2, 1);
    triplot(DT, 'Color', delaunay_color, 'LineWidth', 0.8);
    hold on;
    scatter(positions(:,1), positions(:,2), 20, particle_color, 'filled', ...
        'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.8, 'LineWidth', 1);
    
    axis equal; 
    grid on;
    set(gca, 'GridColor', [0.9, 0.9, 0.9], 'GridAlpha', 0.6);
    set(gca, 'Color', background_color);
    xlim(x_limits); ylim(y_limits);
    
    title('Full Delaunay Triangulation', 'FontSize', 14, 'FontWeight', 'normal', 'Color', [0.2, 0.2, 0.2]);
   xlabel('$x (\mu m)$', 'FontSize', 12, 'FontWeight', 'normal','Interpreter','latex');
    ylabel('$y  (\mu m)$', 'FontSize', 12, 'FontWeight', 'normal','Interpreter','latex');
    
    % Add triangle count text
    text(0.02, 0.98, sprintf('%d triangles', n_triangles), 'Units', 'normalized', ...
        'FontSize', 11, 'FontWeight', 'normal', 'Color', [0.3, 0.3, 0.3], ...
        'BackgroundColor', 'white', 'EdgeColor', [0.8, 0.8, 0.8], 'Margin', 5, ...
        'VerticalAlignment', 'top');
    
    % Right subplot: Distance-Filtered Triangles
    subplot(1, 2, 2);
    scatter(positions(:,1), positions(:,2), 60, particle_color, 'filled', ...
        'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.8, 'LineWidth', 1);
    hold on;
    
    % Draw only the valid triangles
    for i = 1:n_valid_triangles
        p1 = valid_triangles(i, 1);
        p2 = valid_triangles(i, 2);
        p3 = valid_triangles(i, 3);
        
        % Draw the triangle edges with improved styling
        plot([positions(p1,1), positions(p2,1)], [positions(p1,2), positions(p2,2)], ...
            'Color', valid_triangle_color, 'LineWidth', 1.5);
        plot([positions(p2,1), positions(p3,1)], [positions(p2,2), positions(p3,2)], ...
            'Color', valid_triangle_color, 'LineWidth', 1.5);
        plot([positions(p3,1), positions(p1,1)], [positions(p3,2), positions(p1,2)], ...
            'Color', valid_triangle_color, 'LineWidth', 1.5);
    end
    
    axis equal; 
    grid on;
    set(gca, 'GridColor', [0.9, 0.9, 0.9], 'GridAlpha', 0.6);
    set(gca, 'Color', background_color);
    xlim(x_limits); ylim(y_limits);
    
    title('Distance-Filtered Triangulation', 'FontSize', 14, 'FontWeight', 'normal', 'Color', [0.2, 0.2, 0.2]);

    % Add filtered triangle count and threshold text
    text(0.02, 0.98, sprintf('%d triangles', n_valid_triangles), 'Units', 'normalized', ...
        'FontSize', 11, 'FontWeight', 'normal', 'Color', [0.3, 0.3, 0.3], ...
        'BackgroundColor', 'white', 'EdgeColor', [0.8, 0.8, 0.8], 'Margin', 5, ...
        'VerticalAlignment', 'top');
    
    text(0.02, 0.88, sprintf('Threshold: %.1f μm', distance_threshold), 'Units', 'normalized', ...
        'FontSize', 10, 'FontWeight', 'normal', 'Color', [0.5, 0.5, 0.5], ...
        'BackgroundColor', 'white', 'EdgeColor', [0.8, 0.8, 0.8], 'Margin', 5, ...
        'VerticalAlignment', 'top');
    
    % Main title with enhanced styling
    % sgtitle(sprintf('Time: %d | Particles: %d | Efficiency: %.1f%%', ...
    %     current_time, N_particles, 100*n_valid_triangles/n_triangles), ...
    %     'FontSize', 16, 'FontWeight', 'normal', 'Color', [0.1, 0.1, 0.1]);
    % 
    % Improve overall figure appearance
    set(fig, 'Color', 'white');
    
    % Add progress indicator
    if mod(t_idx, 10) == 0 || t_idx == length(analysis_times)
        fprintf('Progress: %d/%d frames (%.1f%%)\n', t_idx, length(analysis_times), 100*t_idx/length(analysis_times));
    end
    
    drawnow;
    
    % Capture frame for video
    frame = getframe(fig);
    writeVideo(v, frame);
end

% Close video file
close(v);
close(fig);

fprintf('\nVideo saved successfully: %s\n', video_path);
fprintf('Video duration: %.2f seconds at %d fps\n', length(analysis_times)/v.FrameRate, v.FrameRate);