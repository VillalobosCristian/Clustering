clear all; clc; close all;

%% Load saved persistence data
[results_file, results_path] = uigetfile('*_persistence_results.mat', 'Select persistence results file');
if isequal(results_file, 0), return; end
fprintf('Loading saved results: %s\n', results_file);
load(fullfile(results_path, results_file));

if ~exist('particle_history', 'var')
    error('particle_history not found! Please re-run the persistence analysis.');
end

fprintf('  particle_history loaded: %d particles\n', length(particle_history));
fprintf('  Analysis frames: %d\n', length(analysis_times));
if exist('cluster_nucleation_time', 'var') && ~isnan(cluster_nucleation_time)
    fprintf('  Cluster nucleation at time: %.0f\n', cluster_nucleation_time);
end

%% Load particle positions
csv_file = strrep(results_file, '_persistence_results.mat', '.csv');
csv_full_path = fullfile(results_path, csv_file);

if ~exist(csv_full_path, 'file')
    [csv_file, csv_path] = uigetfile('*.csv', 'Select original CSV file');
    if isequal(csv_file, 0), return; end
    csv_full_path = fullfile(csv_path, csv_file);
end

%% Quick data reload
fprintf('Loading CSV data...\n');
data = readtable(csv_full_path);
data = data(4:end,:);
if iscell(data.TRACK_ID), data.TRACK_ID = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X), data.POSITION_X = cellfun(@str2double, data.POSITION_X); end
if iscell(data.POSITION_Y), data.POSITION_Y = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T), data.POSITION_T = cellfun(@str2double, data.POSITION_T); end
sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

% Rebuild X, Y, T arrays
X = {}; Y = {}; T = {};
for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    if height(track_data) >= 10
        X{end+1} = track_data.POSITION_X;
        Y{end+1} = track_data.POSITION_Y;
        T{end+1} = track_data.POSITION_T;
    end
end

fprintf('  Loaded %d trajectories\n', length(X));

%% Animation setup
colors = [0.7,0.7,0.7; 0.2,0.4,0.8; 0.9,0.2,0.2]; % Gray, Blue, Red
cluster_color = [1, 0.8, 0];  % Yellow for cluster highlight

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
fig = figure('Position', [100, 100, 1000, 700], 'Color', 'white');

fprintf('\nGenerating animation...\n');

%% Main animation loop
frame_count = 0;
for t_idx = 1:2:length(analysis_times)
    current_time = analysis_times(t_idx);
    
    % Get particles at this time
    positions = [];
    present_particles = [];
    for i = 1:length(X)
        timeIdx = find(T{i} == current_time, 1);
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
            present_particles = [present_particles; i];
        end
    end
    
    if size(positions,1) < 2, continue; end
    
    %% Determine particle states
    particle_colors = repmat(colors(1,:), size(positions,1), 1);
    stable_count = 0;
    transient_count = 0;
    
    for j = 1:length(present_particles)
        particle_idx = present_particles(j);
        
        if particle_history{particle_idx}(t_idx) == 1
            
            % Check persistence
            connected_count = 0;
            current_gap = 0;
            
            for k = t_idx:-1:1
                if particle_history{particle_idx}(k) == 1
                    connected_count = connected_count + 1;
                    current_gap = 0;
                elseif particle_history{particle_idx}(k) == 0
                    current_gap = current_gap + 1;
                    if current_gap > tolerance_frames
                        break;
                    end
                elseif particle_history{particle_idx}(k) == -1
                    current_gap = current_gap + 1;
                    if current_gap > tolerance_frames
                        break;
                    end
                end
            end
            
            % Assign color
            if connected_count >= persistence_frames
                particle_colors(j,:) = colors(3,:); % Red = stable
                stable_count = stable_count + 1;
            else
                particle_colors(j,:) = colors(2,:); % Blue = transient
                transient_count = transient_count + 1;
            end
        end
    end
    
    % Check if particles belong to cluster
    in_cluster = false(size(positions,1), 1);
    cluster_hull_x = [];
    cluster_hull_y = [];
    
    if exist('cluster_particle_ids', 'var') && ~isempty(cluster_particle_ids{t_idx})
        cluster_ids = cluster_particle_ids{t_idx};
        
        for j = 1:length(present_particles)
            if ismember(present_particles(j), cluster_ids)
                in_cluster(j) = true;
            end
        end
        
        % Get cluster hull for visualization
        cluster_pos_idx = find(in_cluster);
        if length(cluster_pos_idx) >= 3
            try
                cluster_positions = positions(cluster_pos_idx, :);
                [k, ~] = convhull(cluster_positions(:,1), cluster_positions(:,2));
                cluster_hull_x = cluster_positions(k, 1);
                cluster_hull_y = cluster_positions(k, 2);
            catch
            end
        end
    end
    
    % Calculate distances for lines
    distances = pdist(positions);
    distance_matrix = squareform(distances);
    
    clf;
    
    % Draw cluster boundary if it exists
    if ~isempty(cluster_hull_x) && cluster_fixed(t_idx) == 1
        fill(cluster_hull_x, cluster_hull_y, cluster_color, ...
            'FaceAlpha', 0.2, 'EdgeColor', cluster_color, 'LineWidth', 2);
        hold on;
    end
    
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
        if in_cluster(i) && cluster_fixed(t_idx) == 1
            % Cluster particles: thicker border
            scatter(positions(i,1), positions(i,2), 80, particle_colors(i,:), 'filled', ...
                'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', cluster_color, 'LineWidth', 2.5);
        else
            % Normal particles
            scatter(positions(i,1), positions(i,2), 60, particle_colors(i,:), 'filled', ...
                'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'white', 'LineWidth', 1);
        end
        hold on;
    end
    
    % Formatting
    axis equal; xlim(xlims); ylim(ylims); grid on;
    xlabel('X Position (\mum)', 'FontSize', 11);
    ylabel('Y Position (\mum)', 'FontSize', 11);
    
    % Title with cluster info
    if cluster_exists(t_idx) == 1
        if cluster_fixed(t_idx) == 1
            cluster_status = sprintf('FIXED CLUSTER: %d particles, %.1f µm²', ...
                cluster_size(t_idx), cluster_area(t_idx));
        else
            cluster_status = sprintf('Cluster forming: %d particles', cluster_size(t_idx));
        end
    else
        cluster_status = 'No cluster detected';
    end
    
    title({sprintf('Time: %.0f | Stable: %d | Transient: %d | Disconnected: %d', ...
        current_time, stable_count, transient_count, ...
        size(positions,1)-(stable_count+transient_count)), cluster_status}, 'FontSize', 11);
    
    % Legend
    h1 = scatter(NaN,NaN,60,colors(3,:),'filled','MarkerEdgeColor','white','LineWidth',1);
    h2 = scatter(NaN,NaN,60,colors(2,:),'filled','MarkerEdgeColor','white','LineWidth',1);
    h3 = scatter(NaN,NaN,60,colors(1,:),'filled','MarkerEdgeColor','white','LineWidth',1);
    if any(cluster_fixed)
        h4 = scatter(NaN,NaN,80,colors(3,:),'filled','MarkerEdgeColor',cluster_color,'LineWidth',2.5);
        legend([h1,h2,h3,h4], {'Stable','Transient','Disconnected','In Cluster'}, ...
            'Location','northeast');
    else
        legend([h1,h2,h3], {'Stable','Transient','Disconnected'}, 'Location','northeast');
    end
    
    drawnow;
    writeVideo(v, getframe(fig));
    
    frame_count = frame_count + 1;
    if mod(frame_count, 20) == 0
        fprintf('  Frame %d/%d (%.1f%%)\n', frame_count, floor(length(analysis_times)/2), ...
            100*frame_count/floor(length(analysis_times)/2));
    end
end

close(v); close(fig);
fprintf('\n=== Animation Complete! ===\n');
fprintf('Saved to: %s\n', fullfile(results_path, video_file));
fprintf('Total frames: %d\n', frame_count);