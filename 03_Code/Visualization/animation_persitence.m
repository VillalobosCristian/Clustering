clear all; clc; close all;

%% Load saved persistence data
[results_file, results_path] = uigetfile('*_persistence_results.mat', 'Select persistence results file');
if isequal(results_file, 0), return; end

fprintf('Loading saved results: %s\n', results_file);
load(fullfile(results_path, results_file));

%% particle positions
csv_file = strrep(results_file, '_persistence_results.mat', '.csv');
csv_path = fullfile(results_path);

if ~exist(csv_path, 'file')
    [csv_file, csv_path] = uigetfile('*.csv', 'Select original CSV file');
    if isequal(csv_file, 0), return; end
end

%% Quick data reload
data = readtable(fullfile(csv_path, csv_file));
data = data(4:end,:);

if iscell(data.TRACK_ID), data.TRACK_ID = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X), data.POSITION_X = cellfun(@str2double, data.POSITION_X); end  
if iscell(data.POSITION_Y), data.POSITION_Y = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T), data.POSITION_T = cellfun(@str2double, data.POSITION_T); end

sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

% Rebuild X, Y, T arrays (minimal version)
X = {}; Y = {}; T = {};
for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    if height(track_data) >= 10
        X{end+1} = track_data.POSITION_X;
        Y{end+1} = track_data.POSITION_Y; 
        T{end+1} = track_data.POSITION_T;
    end
end

%% Animation setup
analysis_times = unique(sorted_data.POSITION_T);
analysis_times = analysis_times(1:1:end); % Match your time_step=5

colors = [0.7,0.7,0.7; 0.2,0.4,0.8; 0.9,0.2,0.2]; % Gray, Blue, Red
state_names = {'Disconnected', 'Transient', 'Stable'};

% Get spatial limits
all_x = []; all_y = [];
for i = 1:length(X), all_x = [all_x; X{i}]; all_y = [all_y; Y{i}]; end
xlims = [min(all_x), max(all_x)] + [-5, 5];
ylims = [min(all_y), max(all_y)] + [-5, 5];

% Video setup
video_file = strrep(results_file, '_persistence_results.mat', '_animation.mp4');
v = VideoWriter(fullfile(results_path, video_file), 'MPEG-4');
v.FrameRate = 8;
open(v);

fig = figure('Position', [100, 100, 800, 600], 'Color', 'white');



%% Main animation loop
for t_idx = 1:2:length(analysis_times)  % Every 2nd frame for speed
    current_time = analysis_times(t_idx);
    
    % Get particles at this time
    positions = []; present_particles = [];
    for i = 1:length(X)
        timeIdx = find(T{i} == current_time, 1);
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
            present_particles = [present_particles; i];
        end
    end
    
    if size(positions,1) < 2, continue; end
    
    % Find time index in saved results  
    saved_t_idx = find(time_points == current_time, 1);
    if isempty(saved_t_idx), continue; end
    
    % Get particle states from saved data
    stable_count = stable_particles_over_time(saved_t_idx);
    transient_count = transient_particles_over_time(saved_t_idx);  
    connected_count = stable_count + transient_count;
    
    % Assign colors (simplified logic)
    particle_colors = repmat(colors(1,:), size(positions,1), 1); % Default: disconnected
    
    if connected_count > 0
        % First 'stable_count' connected particles = stable (red)
        % Next 'transient_count' connected particles = transient (blue)
        conn_idx = 1;
        for j = 1:size(positions,1)
            if conn_idx <= stable_count
                particle_colors(j,:) = colors(3,:); % Red = stable
                conn_idx = conn_idx + 1;
            elseif conn_idx <= connected_count
                particle_colors(j,:) = colors(2,:); % Blue = transient  
                conn_idx = conn_idx + 1;
            end
        end
    end
    
    % Find connections for lines
    distances = pdist(positions);
    distance_matrix = squareform(distances);
    
    clf;
    
    % Draw connection lines
    for i = 1:size(positions,1)
        for j = i+1:size(positions,1)
            if distance_matrix(i,j) <= distance_threshold
                plot([positions(i,1), positions(j,1)], [positions(i,2), positions(j,2)], ...
                     'k-', 'LineWidth', 0.8, 'Color', [0.6, 0.6, 0.6]);
                hold on;
            end
        end
    end
    
    % Draw particles
    for i = 1:size(positions,1)
        scatter(positions(i,1), positions(i,2), 60, particle_colors(i,:), 'filled', ...
               'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'white', 'LineWidth', 1);
        hold on;
    end
    
    % Formatting
    axis equal; xlim(xlims); ylim(ylims); grid on;
    xlabel('X Position'); ylabel('Y Position');
    title(sprintf('Time: %.0f | Stable: %d | Transient: %d | Disconnected: %d', ...
          current_time, stable_count, transient_count, size(positions,1)-connected_count));
    
    % Legend
    h1 = scatter(NaN,NaN,60,colors(3,:),'filled','MarkerEdgeColor','white');
    h2 = scatter(NaN,NaN,60,colors(2,:),'filled','MarkerEdgeColor','white');  
    h3 = scatter(NaN,NaN,60,colors(1,:),'filled','MarkerEdgeColor','white');
    legend([h1,h2,h3], {'Stable','Transient','Disconnected'}, 'Location','northeast');
    
    drawnow;
    writeVideo(v, getframe(fig));
    
    if mod(t_idx, 20) == 0
        fprintf('  Frame %d/%d\n', t_idx, length(analysis_times));
    end
end

close(v); close(fig);

