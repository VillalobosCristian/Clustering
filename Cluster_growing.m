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

%% Parameters for triangle and cluster analysis
calib = 1.0;
distance_threshold = 3;  % Maximum edge length for valid triangles
min_cluster_size = 3;    % Minimum particles needed to form a cluster
overlap_threshold = 0.4; % Minimum overlap ratio for cluster tracking (40% of particles must overlap)

% Initialize global tracking structures
% These will accumulate information across all time points
triangle_database = []; % Format: [time, particle_id_1, particle_id_2, particle_id_3, area, perimeter]
all_cluster_history = []; % Format: [time, cluster_id, size, centroid_x, centroid_y, area, birth_time]
merging_events_log = []; % Format: [time, parent_cluster_1, parent_cluster_2, child_cluster, merged_size]
next_global_cluster_id = 1; % Global counter for assigning unique cluster IDs

% Calculate plot boundaries for consistent visualization across all frames
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

%% Main animation and analysis loop
figure('Position', [100, 100, 1400, 700]);

fprintf('=== CLUSTER MERGING ANALYSIS SYSTEM ===\n');
fprintf('Analyzing %d time points with distance threshold %.1f μm\n', ...
    length(analysis_times), distance_threshold);
fprintf('Minimum cluster size: %d particles\n', min_cluster_size);
fprintf('Cluster overlap threshold: %.1f%%\n', overlap_threshold * 100);
fprintf('\nStarting temporal analysis...\n\n');

% Main temporal loop - this processes each time point sequentially
% building up our understanding of triangle and cluster dynamics
for t_idx = 1:length(analysis_times)
    current_time = analysis_times(t_idx);
    
    %% STEP 1: Extract particle positions at current time (your original approach)
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
    
    % Skip time points with insufficient particles for meaningful analysis
    if N_particles < 3
        fprintf('Time %d: Skipping (only %d particles)\n', current_time, N_particles);
        continue;
    end
    
    %% STEP 2: Create Delaunay triangulation and filter by distance (your original method)
    % This step identifies which particles are close enough to potentially belong to the same cluster
    DT = delaunayTriangulation(positions(:,1), positions(:,2));
    triangles = DT.ConnectivityList;
    n_triangles = size(triangles, 1);
    
    % Filter triangles to keep only those where all edges are below our distance threshold
    % This is the fundamental step that defines what "close enough to cluster" means
    valid_triangles = [];
    for i = 1:n_triangles
        p1 = triangles(i, 1);
        p2 = triangles(i, 2);
        p3 = triangles(i, 3);
        
        % Calculate all three edge lengths for this triangle
        side1 = norm(positions(p1,:) - positions(p2,:));
        side2 = norm(positions(p2,:) - positions(p3,:));
        side3 = norm(positions(p3,:) - positions(p1,:));
        
        max_side = max([side1, side2, side3]);
        
        % Keep triangle only if its longest edge is below our threshold
        % This ensures all particles in the triangle are "close" to each other
        if max_side <= distance_threshold
            valid_triangles = [valid_triangles; triangles(i,:)];
        end
    end
    
    n_valid_triangles = size(valid_triangles, 1);
    
    %% STEP 3: Store triangle information for persistence analysis
    % This builds on your triangle tracking work to create a historical record
    if ~isempty(valid_triangles)
        for i = 1:size(valid_triangles, 1)
            % Extract the three particles that form this triangle
            p1_idx = valid_triangles(i, 1);
            p2_idx = valid_triangles(i, 2);
            p3_idx = valid_triangles(i, 3);
            
            % Convert to particle IDs (which remain constant as particles move)
            p1_id = particle_ids(p1_idx);
            p2_id = particle_ids(p2_idx);
            p3_id = particle_ids(p3_idx);
            
            % Get actual spatial coordinates for geometric calculations
            p1_pos = positions(p1_idx, :);
            p2_pos = positions(p2_idx, :);
            p3_pos = positions(p3_idx, :);
            
            % Calculate triangle area using the cross product formula
            triangle_area = 0.5 * abs((p2_pos(1) - p1_pos(1)) * (p3_pos(2) - p1_pos(2)) - ...
                                      (p3_pos(1) - p1_pos(1)) * (p2_pos(2) - p1_pos(2)));
            
            % Calculate triangle perimeter (total edge length)
            side1 = norm(p1_pos - p2_pos);
            side2 = norm(p2_pos - p3_pos);
            side3 = norm(p3_pos - p1_pos);
            triangle_perimeter = side1 + side2 + side3;
            
            % Create unique triangle fingerprint by sorting particle IDs
            % This allows us to recognize the same triangle relationship across time
            fingerprint = sort([p1_id, p2_id, p3_id]);
            
            % Store in our master triangle database
            triangle_database = [triangle_database; current_time, fingerprint, triangle_area, triangle_perimeter];
        end
    end
    
    %% STEP 4: Detect clusters from triangle network
    % This is where we move from individual triangles to connected cluster structures
    % A cluster is defined as a connected component in the triangle network
    current_clusters = [];
    
    if ~isempty(valid_triangles)
        % Build adjacency matrix showing which particles are connected through triangles
        % Think of this as creating a social network where particles are people
        % and triangle edges represent relationships
        adj_matrix = false(N_particles, N_particles);
        
        % Each triangle creates mutual connections between all three of its vertices
        for i = 1:size(valid_triangles, 1)
            p1 = valid_triangles(i, 1);
            p2 = valid_triangles(i, 2);
            p3 = valid_triangles(i, 3);
            
            % Create bidirectional connections for all pairs in this triangle
            adj_matrix(p1, p2) = true; adj_matrix(p2, p1) = true;
            adj_matrix(p2, p3) = true; adj_matrix(p3, p2) = true;
            adj_matrix(p1, p3) = true; adj_matrix(p3, p1) = true;
        end
        
        % Find connected components using depth-first search (flood fill algorithm)
        % Each connected component represents one distinct cluster
        visited = false(N_particles, 1);
        cluster_count = 0;
        
        for particle_idx = 1:N_particles
            if ~visited(particle_idx)
                % Start exploring a new cluster using depth-first search
                cluster_particles = [];
                stack = particle_idx; % Stack for iterative depth-first search
                
                % Explore all particles reachable from this starting particle
                while ~isempty(stack)
                    current = stack(end);
                    stack(end) = []; % Pop from stack
                    
                    if ~visited(current)
                        visited(current) = true;
                        cluster_particles = [cluster_particles; current];
                        
                        % Add all unvisited neighbors to the exploration stack
                        neighbors = find(adj_matrix(current, :));
                        for neighbor = neighbors
                            if ~visited(neighbor)
                                stack = [stack; neighbor];
                            end
                        end
                    end
                end
                
                % Only accept clusters that meet our minimum size requirement
                % This filters out noise and ensures we focus on meaningful structures
                if length(cluster_particles) >= min_cluster_size
                    cluster_count = cluster_count + 1;
                    
                    % Calculate comprehensive cluster properties
                    cluster_positions = positions(cluster_particles, :);
                    centroid = mean(cluster_positions, 1);
                    
                    % Calculate cluster area using convex hull
                    % This gives us the total spatial footprint of the cluster
                    try
                        hull_idx = convhull(cluster_positions(:,1), cluster_positions(:,2));
                        cluster_area = polyarea(cluster_positions(hull_idx,1), cluster_positions(hull_idx,2));
                    catch
                        cluster_area = 0; % Fallback for degenerate cases
                    end
                    
                    % Create a unique fingerprint for this cluster based on its particles
                    cluster_particle_ids = sort(particle_ids(cluster_particles));
                    
                    % Store all cluster information for tracking and analysis
                    current_clusters(cluster_count).particles = cluster_particles;
                    current_clusters(cluster_count).particle_ids = cluster_particle_ids;
                    current_clusters(cluster_count).positions = cluster_positions;
                    current_clusters(cluster_count).centroid = centroid;
                    current_clusters(cluster_count).size = length(cluster_particles);
                    current_clusters(cluster_count).area = cluster_area;
                    current_clusters(cluster_count).time = current_time;
                    current_clusters(cluster_count).global_id = []; % Will be assigned during tracking
                    current_clusters(cluster_count).birth_time = []; % Will be assigned during tracking
                end
            end
        end
        
        if cluster_count > 0
            cluster_sizes = [current_clusters.size];
            fprintf('Time %d: Found %d clusters (sizes: %s, total particles: %d)\n', ...
                current_time, cluster_count, mat2str(cluster_sizes), sum(cluster_sizes));
        else
            fprintf('Time %d: No clusters detected\n', current_time);
        end
    else
        fprintf('Time %d: No valid triangles found\n', current_time);
    end
    
    %% STEP 5: Track clusters and detect merging events
    % This is the sophisticated part where we maintain cluster identity across time
    % and detect when separate clusters combine into larger structures
    
    if t_idx == 1
        % First time point: simply assign initial IDs to all clusters
        for i = 1:length(current_clusters)
            current_clusters(i).global_id = next_global_cluster_id;
            current_clusters(i).birth_time = current_time;
            next_global_cluster_id = next_global_cluster_id + 1;
        end
        
    else
        % Subsequent time points: track clusters and detect merging
        previous_time = analysis_times(t_idx - 1);
        
        % Extract information about clusters that existed in the previous frame
        % SAFETY CHECK: only proceed if we have cluster history data to examine
        if ~isempty(all_cluster_history)
            previous_cluster_data = all_cluster_history(all_cluster_history(:,1) == previous_time, :);
        else
            previous_cluster_data = [];
        end
        
        if ~isempty(current_clusters) && ~isempty(previous_cluster_data)
            % Detect merging events and assign cluster IDs
            [current_clusters, new_merging_events] = detectMergingEvents(current_clusters, ...
                previous_cluster_data, current_time, next_global_cluster_id, overlap_threshold);
            
            % Update the global cluster ID counter
            if ~isempty(current_clusters)
                next_global_cluster_id = max([current_clusters.global_id]) + 1;
            end
            
            % Log any merging events that occurred
            if ~isempty(new_merging_events)
                merging_events_log = [merging_events_log; new_merging_events];
                
                % Report merging events as they happen
                for i = 1:size(new_merging_events, 1)
                    fprintf('>>> MERGING EVENT at time %d: Clusters %d + %d → Cluster %d (size %d)\n', ...
                        current_time, new_merging_events(i,2), new_merging_events(i,3), ...
                        new_merging_events(i,4), new_merging_events(i,5));
                end
            end
        elseif ~isempty(current_clusters)
            % No previous clusters exist, treat all current clusters as newly born
            for i = 1:length(current_clusters)
                current_clusters(i).global_id = next_global_cluster_id;
                current_clusters(i).birth_time = current_time;
                next_global_cluster_id = next_global_cluster_id + 1;
            end
        end
    end
    
    %% STEP 6: Store cluster information in global history
    % This creates a comprehensive record of all cluster evolution
    for i = 1:length(current_clusters)
        if ~isempty(current_clusters(i).global_id)
            cluster_info = [current_time, current_clusters(i).global_id, current_clusters(i).size, ...
                           current_clusters(i).centroid, current_clusters(i).area, current_clusters(i).birth_time];
            all_cluster_history = [all_cluster_history; cluster_info];
        end
    end
    
    %% STEP 7: Enhanced visualization showing triangles, clusters, and merging
    % This creates a comprehensive view of the multi-scale dynamics
    clf; % Clear figure for next frame
    
    % Panel 1: Triangle age visualization (based on your persistence work)
    subplot(2, 2, 1);
    scatter(positions(:,1), positions(:,2), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    
    % Calculate and visualize triangle ages
    if ~isempty(valid_triangles)
        triangle_ages = calculateTriangleAges(valid_triangles, particle_ids, triangle_database, current_time, time_step);
        
        if ~isempty(triangle_ages)
            max_age = max(triangle_ages);
            min_age = min(triangle_ages);
            age_colormap = hot(max(max_age - min_age + 1, 1));
            
            for i = 1:size(valid_triangles, 1)
                p1 = valid_triangles(i, 1); p2 = valid_triangles(i, 2); p3 = valid_triangles(i, 3);
                triangle_age = triangle_ages(i);
                color_idx = triangle_age - min_age + 1;
                triangle_color = age_colormap(color_idx, :);
                line_thickness = 1 + (triangle_age - min_age) * 0.3;
                
                plot([positions(p1,1), positions(p2,1)], [positions(p1,2), positions(p2,2)], ...
                    'Color', triangle_color, 'LineWidth', line_thickness);
                plot([positions(p2,1), positions(p3,1)], [positions(p2,2), positions(p3,2)], ...
                    'Color', triangle_color, 'LineWidth', line_thickness);
                plot([positions(p3,1), positions(p1,1)], [positions(p3,2), positions(p1,2)], ...
                    'Color', triangle_color, 'LineWidth', line_thickness);
            end
        end
    end
    
    axis equal; grid on; xlim(x_limits); ylim(y_limits);
    title('Triangle Persistence Network');
    xlabel('x (μm)'); ylabel('y (μm)');
    
    % Panel 2: Cluster identification with unique colors
    subplot(2, 2, 2);
    scatter(positions(:,1), positions(:,2), 20, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    
    if ~isempty(current_clusters)
        cluster_colors = lines(length(current_clusters));
        
        for i = 1:length(current_clusters)
            cluster_pos = current_clusters(i).positions;
            cluster_color = cluster_colors(i, :);
            
            % Draw cluster particles with unique color
            scatter(cluster_pos(:,1), cluster_pos(:,2), 80, cluster_color, 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            
            % Draw cluster boundary (convex hull)
            try
                hull_idx = convhull(cluster_pos(:,1), cluster_pos(:,2));
                plot(cluster_pos(hull_idx,1), cluster_pos(hull_idx,2), ...
                    'Color', cluster_color, 'LineWidth', 2);
            catch
                % Skip boundary for degenerate clusters
            end
            
            % Label cluster with ID and size
            centroid = current_clusters(i).centroid;
            text(centroid(1), centroid(2), sprintf('C%d\n(%d)', current_clusters(i).global_id, current_clusters(i).size), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'white', ...
                'BackgroundColor', cluster_color, 'EdgeColor', 'k');
        end
    end
    
    axis equal; grid on; xlim(x_limits); ylim(y_limits);
    title(sprintf('Cluster Detection (%d clusters)', length(current_clusters)));
    xlabel('x (μm)'); ylabel('y (μm)');
    
    % Panel 3: Cluster size evolution over time
    subplot(2, 2, 3);
    if ~isempty(all_cluster_history)
        % Plot size evolution for each cluster that has existed
        unique_cluster_ids = unique(all_cluster_history(:,2));
        cluster_colors_timeline = lines(length(unique_cluster_ids));
        
        for i = 1:length(unique_cluster_ids)
            cluster_id = unique_cluster_ids(i);
            cluster_data = all_cluster_history(all_cluster_history(:,2) == cluster_id, :);
            
            times = cluster_data(:, 1);
            sizes = cluster_data(:, 3);
            
            plot(times, sizes, 'o-', 'Color', cluster_colors_timeline(i, :), ...
                'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', cluster_colors_timeline(i, :));
            hold on;
            
            % Highlight current time point
            current_idx = find(times == current_time);
            if ~isempty(current_idx)
                plot(current_time, sizes(current_idx), 'o', 'Color', cluster_colors_timeline(i, :), ...
                    'MarkerSize', 12, 'LineWidth', 3);
            end
        end
        
        % Mark merging events
        if ~isempty(merging_events_log)
            merge_times = merging_events_log(:, 1);
            for i = 1:length(merge_times)
                xline(merge_times(i), 'r--', 'LineWidth', 2, 'Alpha', 0.7);
            end
        end
    end
    
    xlabel('Time'); ylabel('Cluster Size (particles)');
    title('Cluster Growth & Merging Timeline');
    grid on;
    
    % Panel 4: Merging events summary
    subplot(2, 2, 4);
    if ~isempty(merging_events_log)
        % Create a bar chart of merging events over time
        merge_times = merging_events_log(:, 1);
        merge_sizes = merging_events_log(:, 5);
        
        [unique_times, ~, idx] = unique(merge_times);
        merge_counts = accumarray(idx, 1);
        
        bar(unique_times, merge_counts, 'FaceColor', [0.8 0.4 0.4], 'EdgeColor', 'k');
        
        xlabel('Time'); ylabel('Number of Merging Events');
        title('Merging Event Timeline');
        grid on;
    else
        text(0.5, 0.5, 'No merging events detected yet', 'HorizontalAlignment', 'center', ...
            'Units', 'normalized', 'FontSize', 12);
        title('Merging Event Timeline');
    end
    
    % Overall figure title with comprehensive status
    sgtitle(sprintf('Time = %d | Particles = %d | Triangles = %d | Clusters = %d | Total Merges = %d', ...
        current_time, N_particles, n_valid_triangles, length(current_clusters), size(merging_events_log, 1)));
    
    drawnow;
    pause(0.15); % Control animation speed
end

%% COMPREHENSIVE POST-ANALYSIS REPORT
% After processing all time points, analyze the complete merging dynamics

fprintf('\n' + string(repmat('=', 1, 60)) + '\n');
fprintf('COMPREHENSIVE CLUSTER MERGING ANALYSIS REPORT\n');
fprintf(string(repmat('=', 1, 60)) + '\n');

fprintf('\nDataset Overview:\n');
fprintf('  Time range analyzed: %d - %d (%d time points)\n', ...
    min(analysis_times), max(analysis_times), length(analysis_times));
fprintf('  Total triangle observations: %d\n', size(triangle_database, 1));
fprintf('  Total cluster observations: %d\n', size(all_cluster_history, 1));

if ~isempty(all_cluster_history)
    % Cluster statistics
    unique_cluster_ids = unique(all_cluster_history(:,2));
    fprintf('\nCluster Statistics:\n');
    fprintf('  Total unique clusters detected: %d\n', length(unique_cluster_ids));
    
    cluster_lifespans = [];
    cluster_max_sizes = [];
    
    for i = 1:length(unique_cluster_ids)
        cluster_id = unique_cluster_ids(i);
        cluster_data = all_cluster_history(all_cluster_history(:,2) == cluster_id, :);
        
        lifespan = size(cluster_data, 1);
        max_size = max(cluster_data(:, 3));
        
        cluster_lifespans = [cluster_lifespans; lifespan];
        cluster_max_sizes = [cluster_max_sizes; max_size];
    end
    
    fprintf('  Average cluster lifespan: %.1f frames\n', mean(cluster_lifespans));
    fprintf('  Average maximum cluster size: %.1f particles\n', mean(cluster_max_sizes));
    fprintf('  Largest cluster observed: %d particles\n', max(cluster_max_sizes));
end

% Merging event analysis
if ~isempty(merging_events_log)
    fprintf('\nMerging Event Analysis:\n');
    fprintf('  Total merging events: %d\n', size(merging_events_log, 1));
    
    merge_sizes = merging_events_log(:, 5);
    fprintf('  Average size of merged clusters: %.1f particles\n', mean(merge_sizes));
    fprintf('  Largest merged cluster: %d particles\n', max(merge_sizes));
    
    % Temporal distribution of merging
    merge_times = merging_events_log(:, 1);
    fprintf('  Merging events span: time %d to %d\n', min(merge_times), max(merge_times));
    
    % Calculate merging rate
    time_span = max(analysis_times) - min(analysis_times);
    merge_rate = size(merging_events_log, 1) / (time_span / time_step);
    fprintf('  Average merging rate: %.3f events per frame\n', merge_rate);
    
    fprintf('\nDetailed Merging Events:\n');
    for i = 1:size(merging_events_log, 1)
        fprintf('  Time %d: Clusters %d + %d → Cluster %d (final size: %d particles)\n', ...
            merging_events_log(i,1), merging_events_log(i,2), merging_events_log(i,3), ...
            merging_events_log(i,4), merging_events_log(i,5));
    end
else
    fprintf('\nMerging Event Analysis:\n');
    fprintf('  No merging events detected during analysis period.\n');
    fprintf('  Consider adjusting distance_threshold or overlap_threshold parameters.\n');
end

fprintf('\n' + string(repmat('=', 1, 60)) + '\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf(string(repmat('=', 1, 60)) + '\n');

%% HELPER FUNCTIONS FOR MERGING DETECTION AND ANALYSIS

function triangle_ages = calculateTriangleAges(valid_triangles, particle_ids, triangle_database, current_time, time_step)
    % Calculate how long each current triangle has existed
    triangle_ages = [];
    
    for i = 1:size(valid_triangles, 1)
        p1_idx = valid_triangles(i, 1);
        p2_idx = valid_triangles(i, 2);
        p3_idx = valid_triangles(i, 3);
        
        current_fingerprint = sort([particle_ids(p1_idx), particle_ids(p2_idx), particle_ids(p3_idx)]);
        
        % Find earliest occurrence of this triangle
        triangle_birth_time = current_time;
        for db_row = 1:size(triangle_database, 1)
            db_time = triangle_database(db_row, 1);
            db_fingerprint = triangle_database(db_row, 2:4);
            
            if isequal(sort(db_fingerprint), current_fingerprint)
                triangle_birth_time = min(triangle_birth_time, db_time);
            end
        end
        
        triangle_age = (current_time - triangle_birth_time) / time_step + 1;
        triangle_ages = [triangle_ages; triangle_age];
    end
end

function [tracked_clusters, merging_events] = detectMergingEvents(current_clusters, previous_cluster_data, current_time, next_global_id, overlap_threshold)
    % Detect when multiple previous clusters merge into a single current cluster
    % Returns updated current_clusters with proper IDs and a list of merging events
    
    merging_events = [];
    
    % Extract previous cluster information
    prev_ids = unique(previous_cluster_data(:, 2));
    prev_centroids = [];
    prev_sizes = [];
    
    for i = 1:length(prev_ids)
        cluster_data = previous_cluster_data(previous_cluster_data(:, 2) == prev_ids(i), :);
        prev_centroids = [prev_centroids; cluster_data(1, 4:5)];
        prev_sizes = [prev_sizes; cluster_data(1, 3)];
    end
    
    % Track each current cluster
    for i = 1:length(current_clusters)
        current_centroid = current_clusters(i).centroid;
        current_size = current_clusters(i).size;
        
        % Calculate spatial distances to all previous clusters
        if ~isempty(prev_centroids)
            distances = sqrt(sum((prev_centroids - current_centroid).^2, 2));
            
            % Find previous clusters that are spatially close and size-compatible
            max_distance = 10; % Maximum distance for considering clusters as potentially related
            close_clusters = find(distances < max_distance);
            
            if length(close_clusters) >= 2
                % Potential merging event: multiple previous clusters are close to this current cluster
                % Check if the combined size makes sense
                combined_prev_size = sum(prev_sizes(close_clusters));
                size_ratio = current_size / combined_prev_size;
                
                % Accept as merging event if size ratio is reasonable (between 0.7 and 1.3)
                if size_ratio >= 0.7 && size_ratio <= 1.3
                    % Record merging event
                    if length(close_clusters) == 2
                        % Simple two-cluster merge
                        parent1_id = prev_ids(close_clusters(1));
                        parent2_id = prev_ids(close_clusters(2));
                        child_id = next_global_id;
                        
                        merging_events = [merging_events; current_time, parent1_id, parent2_id, child_id, current_size];
                        
                        current_clusters(i).global_id = child_id;
                        current_clusters(i).birth_time = current_time;
                        next_global_id = next_global_id + 1;
                    else
                        % Multi-cluster merge (more complex, assign new ID)
                        child_id = next_global_id;
                        
                        % For multi-mergers, record the two largest parent clusters
                        [~, sorted_idx] = sort(prev_sizes(close_clusters), 'descend');
                        parent1_id = prev_ids(close_clusters(sorted_idx(1)));
                        parent2_id = prev_ids(close_clusters(sorted_idx(2)));
                        
                        merging_events = [merging_events; current_time, parent1_id, parent2_id, child_id, current_size];
                        
                        current_clusters(i).global_id = child_id;
                        current_clusters(i).birth_time = current_time;
                        next_global_id = next_global_id + 1;
                    end
                else
                    % Size mismatch - treat as continuation of largest nearby cluster
                    [~, closest_idx] = min(distances);
                    current_clusters(i).global_id = prev_ids(closest_idx);
                    
                    % Find birth time from previous data
                    birth_data = previous_cluster_data(previous_cluster_data(:, 2) == prev_ids(closest_idx), :);
                    current_clusters(i).birth_time = birth_data(1, 7);
                end
            elseif length(close_clusters) == 1
                % Single cluster continuation
                current_clusters(i).global_id = prev_ids(close_clusters(1));
                
                % Find birth time from previous data
                birth_data = previous_cluster_data(previous_cluster_data(:, 2) == prev_ids(close_clusters(1)), :);
                current_clusters(i).birth_time = birth_data(1, 7);
            else
                % New cluster (no close previous clusters)
                current_clusters(i).global_id = next_global_id;
                current_clusters(i).birth_time = current_time;
                next_global_id = next_global_id + 1;
            end
        else
            % No previous clusters exist - assign new ID
            current_clusters(i).global_id = next_global_id;
            current_clusters(i).birth_time = current_time;
            next_global_id = next_global_id + 1;
        end
    end
    
    tracked_clusters = current_clusters;
end