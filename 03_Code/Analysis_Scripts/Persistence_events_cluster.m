clear all; clc; close all;

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

%% Parameters 
calib = 1.0;                    % µm/pixel
distance_threshold = 5.0;       % µm (connection distance)
time_step = 1;                  % frames (temporal sampling)
persistence_frames = 5;         % frames (minimum for "stable")
tolerance_frames = 1;           % frames (allowed gaps)

% Cluster parameters
min_cluster_size = 8;           % Minimum stable particles to form a cluster
min_cluster_persistence = 10;   % Frames cluster must exist to be "fixed"

%% Get all unique time points
all_times = [];
for i = 1:length(T)
    all_times = [all_times; T{i}];
end
unique_times = unique(all_times);
analysis_times = unique_times(1:time_step:end);

%% Create particle history storage
particle_history = cell(length(X), 1);
for i = 1:length(X)
    particle_history{i} = -ones(length(analysis_times), 1);
end

connected_particles = [];
total_particles = [];
connection_ratio = [];
time_points = [];

%% MAIN LOOP: Calculate connections
fprintf('Calculating particle connections...\n');
for t_idx = 1:length(analysis_times)
    current_time = analysis_times(t_idx);
    
    % Extract positions at current time 
    positions = [];
    present_particles = [];
    
    for i = 1:length(X)
        timeIdx = find(T{i} == current_time, 1);
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
            present_particles = [present_particles; i];
        end
    end
    
    positions = positions * calib;
    N_particles = size(positions, 1);
    
    if N_particles < 2
        continue;
    end
    
    % Calculate pairwise distances
    distances = pdist(positions);
    distance_matrix = squareform(distances);  
    
    % Check connections
    connected_count = 0;
    for j = 1:N_particles
        particle_idx = present_particles(j);
        
        neighbors = distance_matrix(j, :);
        neighbors(j) = inf;
        
        if any(neighbors <= distance_threshold)
            connected_count = connected_count + 1;
            particle_history{particle_idx}(t_idx) = 1;
        else
            particle_history{particle_idx}(t_idx) = 0;
        end
    end
    
    % Store results
    time_points = [time_points; current_time];
    connected_particles = [connected_particles; connected_count];
    total_particles = [total_particles; N_particles];
    connection_ratio = [connection_ratio; connected_count / N_particles];
    
    if mod(t_idx, 50) == 0
        fprintf('  Progress: %d/%d frames\n', t_idx, length(analysis_times));
    end
end

%% Calculate persistence (stable vs transient)
fprintf('Calculating persistence...\n');
stable_particles_over_time = [];
transient_particles_over_time = [];
stable_particle_indices = cell(length(analysis_times), 1);

for t_idx = 1:length(analysis_times)
    current_time = analysis_times(t_idx);
    
    present_particles = [];
    for i = 1:length(X)
        timeIdx = find(T{i} == current_time, 1);
        if ~isempty(timeIdx)
            present_particles = [present_particles; i];
        end
    end
    
    if length(present_particles) < 2
        stable_particles_over_time = [stable_particles_over_time; 0];
        transient_particles_over_time = [transient_particles_over_time; 0];
        stable_particle_indices{t_idx} = [];
        continue;
    end
    
    stable_count = 0;
    transient_count = 0;
    stable_list = [];
    
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
            
            % Classify
            if connected_count >= persistence_frames
                stable_count = stable_count + 1;
                stable_list = [stable_list; particle_idx];
            else
                transient_count = transient_count + 1;
            end
        end
    end
    
    stable_particles_over_time = [stable_particles_over_time; stable_count];
    transient_particles_over_time = [transient_particles_over_time; transient_count];
    stable_particle_indices{t_idx} = stable_list;
end

%% CLUSTER DETECTION
fprintf('\n=== Detecting stable clusters ===\n');

cluster_exists = zeros(length(analysis_times), 1);
cluster_size = zeros(length(analysis_times), 1);
cluster_area = zeros(length(analysis_times), 1);
cluster_particle_ids = cell(length(analysis_times), 1);

for t_idx = 1:length(analysis_times)
    
    stable_list = stable_particle_indices{t_idx};
    
    if length(stable_list) < min_cluster_size
        continue;
    end
    
    % Get positions of stable particles
    current_time = analysis_times(t_idx);
    stable_positions = [];
    stable_pos_to_id = [];
    
    for i = 1:length(stable_list)
        particle_idx = stable_list(i);
        timeIdx = find(T{particle_idx} == current_time, 1);
        if ~isempty(timeIdx)
            stable_positions = [stable_positions; X{particle_idx}(timeIdx), Y{particle_idx}(timeIdx)];
            stable_pos_to_id = [stable_pos_to_id; particle_idx];
        end
    end
    
    if size(stable_positions, 1) < min_cluster_size
        continue;
    end
    
    stable_positions = stable_positions * calib;
    
    % Find connected components among stable particles
    distances = pdist(stable_positions);
    distance_matrix = squareform(distances);
    adjacency = distance_matrix <= distance_threshold;
    adjacency(logical(eye(size(adjacency)))) = 0;
    
    % Use graph to find connected components
    G = graph(adjacency);
    components = conncomp(G);
    
    % Find largest component
    [component_sizes, component_ids] = groupcounts(components');
    [max_size, max_idx] = max(component_sizes);
    
    if max_size >= min_cluster_size
        % We have a cluster!
        largest_component_id = component_ids(max_idx);
        cluster_particle_idx = find(components == largest_component_id);
        cluster_pos = stable_positions(cluster_particle_idx, :);
        cluster_ids = stable_pos_to_id(cluster_particle_idx);
        
        % Calculate area using convex hull
        if size(cluster_pos, 1) >= 3
            try
                [k, hull_area] = convhull(cluster_pos(:,1), cluster_pos(:,2));
                cluster_area(t_idx) = hull_area;
            catch
                cluster_area(t_idx) = 0;
            end
        end
        
        cluster_exists(t_idx) = 1;
        cluster_size(t_idx) = max_size;
        cluster_particle_ids{t_idx} = cluster_ids;
    end
end

% Identify when cluster becomes "fixed"
cluster_fixed = zeros(length(analysis_times), 1);
cluster_nucleation_time = NaN;

for t_idx = min_cluster_persistence:length(analysis_times)
    if all(cluster_exists(t_idx-min_cluster_persistence+1:t_idx) == 1)
        cluster_fixed(t_idx) = 1;
        
        if isnan(cluster_nucleation_time)
            cluster_nucleation_time = time_points(t_idx);
            fprintf('  Cluster nucleation detected at time %.0f (frame %d)\n', ...
                cluster_nucleation_time, t_idx);
        end
    end
end

%% Plots
stable_ratio = stable_particles_over_time ./ max(total_particles, 1);
transient_ratio = transient_particles_over_time ./ max(total_particles, 1);

figure('Position', [100, 100, 1600, 1000], 'Color', 'w');

% Plot 1: Particle counts
subplot(3, 3, 1);
plot(time_points, total_particles, 'k--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
hold on;
plot(time_points, connected_particles, 'b-', 'LineWidth', 2);
plot(time_points, stable_particles_over_time, 'r-', 'LineWidth', 2);
plot(time_points, transient_particles_over_time, 'g--', 'LineWidth', 1.5);
if ~isnan(cluster_nucleation_time)
    xline(cluster_nucleation_time, 'k--', 'LineWidth', 2, 'Label', 'Nucleation');
end
xlabel('Time', 'FontSize', 12);
ylabel('Number of Particles', 'FontSize', 12);
legend({'Total', 'Connected', 'Stable', 'Transient'}, 'Location', 'best');
title('Particle Counts with Persistence');
grid on; box on;

% Plot 2: Connection percentages
subplot(3, 3, 2);
plot(time_points, connection_ratio * 100, 'b-', 'LineWidth', 2);
hold on;
plot(time_points, stable_ratio * 100, 'r-', 'LineWidth', 2);
plot(time_points, transient_ratio * 100, 'g--', 'LineWidth', 1.5);
if ~isnan(cluster_nucleation_time)
    xline(cluster_nucleation_time, 'k--', 'LineWidth', 2);
end
xlabel('Time', 'FontSize', 12);
ylabel('Percentage (%)', 'FontSize', 12);
legend({'Total Connected', 'Stable', 'Transient'}, 'Location', 'best');
title('Connection Percentages');
ylim([0 100]);
grid on; box on;

% Plot 3: Stability ratio
subplot(3, 3, 3);
stability_index = stable_particles_over_time ./ max(connected_particles, 1);
plot(time_points, stability_index * 100, 'm', 'LineWidth', 2);
if ~isnan(cluster_nucleation_time)
    xline(cluster_nucleation_time, 'k--', 'LineWidth', 2);
end
xlabel('Time', 'FontSize', 12);
ylabel('Stability Index (%)', 'FontSize', 12);
title('Stability Index (Stable/Connected)');
ylim([0 100]);
grid on; box on;

% Plot 4: Cluster size over time
subplot(3, 3, 4);
plot(time_points, cluster_size, 'r-', 'LineWidth', 2);
hold on;
plot(time_points(cluster_fixed==1), cluster_size(cluster_fixed==1), 'ko', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 4);
yline(min_cluster_size, 'k--', 'LineWidth', 1);
if ~isnan(cluster_nucleation_time)
    xline(cluster_nucleation_time, 'k--', 'LineWidth', 2);
end
xlabel('Time', 'FontSize', 12);
ylabel('Cluster Size (particles)', 'FontSize', 12);
title(sprintf('Stable Cluster Size (threshold=%d)', min_cluster_size));
grid on; box on;

% Plot 5: Cluster area over time
subplot(3, 3, 5);
plot(time_points, cluster_area, 'b-', 'LineWidth', 2);
hold on;
plot(time_points(cluster_fixed==1), cluster_area(cluster_fixed==1), 'ko', ...
    'MarkerFaceColor', 'b', 'MarkerSize', 4);
if ~isnan(cluster_nucleation_time)
    xline(cluster_nucleation_time, 'k--', 'LineWidth', 2);
end
xlabel('Time', 'FontSize', 12);
ylabel('Cluster Area (\mum^2)', 'FontSize', 12);
title('Stable Cluster Area');
grid on; box on;

% Plot 6: Cluster existence
subplot(3, 3, 6);
area_plot = area(time_points, [cluster_exists, cluster_fixed-cluster_exists]);
area_plot(1).FaceColor = [0.8 0.8 0.8];
area_plot(2).FaceColor = [0.2 0.8 0.2];
area_plot(1).EdgeColor = 'none';
area_plot(2).EdgeColor = 'none';
hold on;
if ~isnan(cluster_nucleation_time)
    xline(cluster_nucleation_time, 'k--', 'LineWidth', 2);
end
xlabel('Time', 'FontSize', 12);
ylabel('Cluster State', 'FontSize', 12);
title('Cluster Detection');
legend({'Cluster Exists', 'Cluster Fixed'}, 'Location', 'best');
ylim([0 1.2]);
grid on; box on;

% Plot 7: Cluster size vs time (after nucleation)
if ~isnan(cluster_nucleation_time)
    subplot(3, 3, 7);
    nucleation_idx = find(time_points == cluster_nucleation_time, 1);
    post_nucleation_time = time_points(nucleation_idx:end) - cluster_nucleation_time;
    post_nucleation_size = cluster_size(nucleation_idx:end);
    
    plot(post_nucleation_time, post_nucleation_size, 'r-', 'LineWidth', 2);
    xlabel('Time after nucleation', 'FontSize', 12);
    ylabel('Cluster Size (particles)', 'FontSize', 12);
    title('Cluster Growth After Nucleation');
    grid on; box on;
    
    % Plot 8: Cluster area vs time (after nucleation)
    subplot(3, 3, 8);
    post_nucleation_area = cluster_area(nucleation_idx:end);
    plot(post_nucleation_time, post_nucleation_area, 'b-', 'LineWidth', 2);
    xlabel('Time after nucleation', 'FontSize', 12);
    ylabel('Cluster Area (\mum^2)', 'FontSize', 12);
    title('Cluster Area Growth After Nucleation');
    grid on; box on;
    
    % Plot 9: Cluster density (size/area)
    subplot(3, 3, 9);
    cluster_density = cluster_size(nucleation_idx:end) ./ max(cluster_area(nucleation_idx:end), 1);
    plot(post_nucleation_time, cluster_density, 'g-', 'LineWidth', 2);
    xlabel('Time after nucleation', 'FontSize', 12);
    ylabel('Cluster Density (particles/\mum^2)', 'FontSize', 12);
    title('Cluster Packing Density');
    grid on; box on;
end

%% Summary statistics
fprintf('\n=== Summary Statistics ===\n');
fprintf('Total frames analyzed: %d\n', length(analysis_times));
fprintf('Cluster detected in %d frames (%.1f%%)\n', sum(cluster_exists), ...
    100*sum(cluster_exists)/length(analysis_times));
fprintf('Cluster fixed in %d frames (%.1f%%)\n', sum(cluster_fixed), ...
    100*sum(cluster_fixed)/length(analysis_times));
if ~isnan(cluster_nucleation_time)
    fprintf('Nucleation time: %.0f\n', cluster_nucleation_time);
    max_size = max(cluster_size);
    max_area = max(cluster_area);
    fprintf('Maximum cluster size: %d particles\n', max_size);
    fprintf('Maximum cluster area: %.1f µm²\n', max_area);
    if max_area > 0
        fprintf('Average cluster density: %.2f particles/µm²\n', max_size/max_area);
    end
end

%% Save results
results_file = [file_name(1:end-4), '_persistence_results.mat'];
save(fullfile(path_name, results_file), 'time_points', 'stable_particles_over_time', ...
     'transient_particles_over_time', 'connection_ratio', 'stable_ratio', 'transient_ratio', ...
     'persistence_frames', 'tolerance_frames', 'distance_threshold', ...
     'particle_history', 'analysis_times', 'stable_particle_indices', ...
     'cluster_exists', 'cluster_size', 'cluster_area', 'cluster_particle_ids', ...
     'cluster_fixed', 'cluster_nucleation_time', 'min_cluster_size', 'min_cluster_persistence');

fprintf('\nResults saved to: %s\n', results_file);
fprintf('Ready for animation!\n');