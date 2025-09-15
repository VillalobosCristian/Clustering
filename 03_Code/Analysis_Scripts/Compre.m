clear all; clc; close all;

%% COMPREHENSIVE TRACKING DATA EXTRACTOR
% This script extracts all key metrics from your tracking data
% Run this on different power conditions and share the results!

%% Load CSV data
[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), return; end

fprintf('=== COMPREHENSIVE TRACKING ANALYSIS ===\n');
fprintf('Analyzing: %s\n', file_name);
fprintf('========================================\n\n');

% Parse experimental conditions from filename
exp_info = parse_filename(file_name);
fprintf('Experimental Conditions:\n');
fprintf('  Magnification: %s\n', exp_info.magnification);
fprintf('  Power Setting: %s\n', exp_info.power);
fprintf('  Trial Number: %s\n', exp_info.trial);
fprintf('\n');

data = readtable(fullfile(path_name, file_name));
data = data(4:end,:); % Remove TrackMate headers

% Convert to numeric if needed
if iscell(data.TRACK_ID), data.TRACK_ID = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X), data.POSITION_X = cellfun(@str2double, data.POSITION_X); end
if iscell(data.POSITION_Y), data.POSITION_Y = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T), data.POSITION_T = cellfun(@str2double, data.POSITION_T); end

sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

%% Extract trajectories
min_track_length = 50; % Lower threshold to capture more data
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

%% Basic Dataset Characteristics
all_times = [];
for i = 1:length(T)
    all_times = [all_times; T{i}];
end
unique_times = unique(all_times);

% Dataset summary
results.experiment_info = exp_info;
results.basic_stats.filename = file_name;
results.basic_stats.total_tracks = count;
results.basic_stats.time_range = [min(unique_times), max(unique_times)];
results.basic_stats.total_timepoints = length(unique_times);
results.basic_stats.time_duration = max(unique_times) - min(unique_times);

fprintf('BASIC DATASET CHARACTERISTICS:\n');
fprintf('  Valid trajectories: %d\n', count);
fprintf('  Time range: %d to %d (%d points)\n', min(unique_times), max(unique_times), length(unique_times));
fprintf('  Duration: %d time units\n', max(unique_times) - min(unique_times));

%% Spatial characteristics
all_x = []; all_y = [];
for i = 1:length(X)
    all_x = [all_x; X{i}];
    all_y = [all_y; Y{i}];
end

results.spatial.field_of_view = [max(all_x) - min(all_x), max(all_y) - min(all_y)];
results.spatial.particle_density = length(all_x) / ((max(all_x)-min(all_x)) * (max(all_y)-min(all_y)));
results.spatial.center_of_mass = [mean(all_x), mean(all_y)];

fprintf('  Field of view: %.1f × %.1f pixels\n', results.spatial.field_of_view);
fprintf('  Overall particle density: %.4f particles/pixel²\n', results.spatial.particle_density);
fprintf('\n');

%% Time-resolved analysis
analysis_parameters.distance_threshold = 5.0; % μm
analysis_parameters.calib = 1.0;
analysis_parameters.time_step = 10; % Analyze every 10th frame for speed

analysis_times = unique_times(1:analysis_parameters.time_step:end);

% Initialize time series data
time_series = struct();
time_series.time_points = [];
time_series.n_particles = [];
time_series.connected_particles = [];
time_series.n_clusters = [];
time_series.largest_cluster_size = [];
time_series.total_clustered_particles = [];
time_series.cluster_sizes = {};
time_series.cluster_areas = [];

fprintf('TIME-RESOLVED ANALYSIS:\n');
fprintf('  Distance threshold: %.1f\n', analysis_parameters.distance_threshold);
fprintf('  Analyzing %d time points...\n', length(analysis_times));

for t_idx = 1:length(analysis_times)
    current_time = analysis_times(t_idx);
    
    % Extract positions at current time
    positions = [];
    for i = 1:length(X)
        timeIdx = find(T{i} == current_time, 1);
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
        end
    end
    
    positions = positions * analysis_parameters.calib;
    N_particles = size(positions, 1);
    
    if N_particles < 3
        continue;
    end
    
    % Connected particles analysis
    distances = pdist(positions);
    distance_matrix = squareform(distances);
    
    connected_count = 0;
    for i = 1:N_particles
        neighbors = distance_matrix(i, :);
        neighbors(i) = inf;
        if any(neighbors <= analysis_parameters.distance_threshold)
            connected_count = connected_count + 1;
        end
    end
    
    % Cluster detection using Delaunay + distance filtering
    if N_particles >= 3
        DT = delaunayTriangulation(positions(:,1), positions(:,2));
        triangles = DT.ConnectivityList;
        
        % Filter triangles by distance
        valid_triangles = [];
        for i = 1:size(triangles, 1)
            p1 = triangles(i, 1); p2 = triangles(i, 2); p3 = triangles(i, 3);
            side1 = norm(positions(p1,:) - positions(p2,:));
            side2 = norm(positions(p2,:) - positions(p3,:));
            side3 = norm(positions(p3,:) - positions(p1,:));
            
            if max([side1, side2, side3]) <= analysis_parameters.distance_threshold
                valid_triangles = [valid_triangles; triangles(i,:)];
            end
        end
        
        % Find connected components (clusters)
        if ~isempty(valid_triangles)
            adj_matrix = false(N_particles, N_particles);
            for i = 1:size(valid_triangles, 1)
                p1 = valid_triangles(i, 1); p2 = valid_triangles(i, 2); p3 = valid_triangles(i, 3);
                adj_matrix(p1, p2) = true; adj_matrix(p2, p1) = true;
                adj_matrix(p2, p3) = true; adj_matrix(p3, p2) = true;
                adj_matrix(p1, p3) = true; adj_matrix(p3, p1) = true;
            end
            
            % Find clusters using connected components
            visited = false(N_particles, 1);
            clusters = [];
            
            for particle_idx = 1:N_particles
                if ~visited(particle_idx)
                    cluster_particles = [];
                    stack = particle_idx;
                    
                    while ~isempty(stack)
                        current = stack(end);
                        stack(end) = [];
                        
                        if ~visited(current)
                            visited(current) = true;
                            cluster_particles = [cluster_particles; current];
                            neighbors = find(adj_matrix(current, :));
                            for neighbor = neighbors
                                if ~visited(neighbor)
                                    stack = [stack; neighbor];
                                end
                            end
                        end
                    end
                    
                    if length(cluster_particles) >= 3
                        cluster_pos = positions(cluster_particles, :);
                        try
                            hull_idx = convhull(cluster_pos(:,1), cluster_pos(:,2));
                            cluster_area = polyarea(cluster_pos(hull_idx,1), cluster_pos(hull_idx,2));
                        catch
                            cluster_area = 0;
                        end
                        
                        clusters(end+1).size = length(cluster_particles);
                        clusters(end).area = cluster_area;
                        clusters(end).positions = cluster_pos;
                    end
                end
            end
        else
            clusters = [];
        end
    else
        clusters = [];
    end
    
    % Store time series data
    time_series.time_points = [time_series.time_points; current_time];
    time_series.n_particles = [time_series.n_particles; N_particles];
    time_series.connected_particles = [time_series.connected_particles; connected_count];
    
    if isempty(clusters)
        time_series.n_clusters = [time_series.n_clusters; 0];
        time_series.largest_cluster_size = [time_series.largest_cluster_size; 0];
        time_series.total_clustered_particles = [time_series.total_clustered_particles; 0];
        time_series.cluster_sizes{end+1} = [];
        time_series.cluster_areas = [time_series.cluster_areas; 0];
    else
        cluster_sizes = [clusters.size];
        cluster_areas = [clusters.area];
        
        time_series.n_clusters = [time_series.n_clusters; length(clusters)];
        time_series.largest_cluster_size = [time_series.largest_cluster_size; max(cluster_sizes)];
        time_series.total_clustered_particles = [time_series.total_clustered_particles; sum(cluster_sizes)];
        time_series.cluster_sizes{end+1} = cluster_sizes;
        time_series.cluster_areas = [time_series.cluster_areas; sum(cluster_areas)];
    end
    
    if mod(t_idx, 20) == 0 || t_idx == length(analysis_times)
        fprintf('    Progress: %d/%d (Time %d: %d particles, %d clusters)\n', ...
            t_idx, length(analysis_times), current_time, N_particles, length(clusters));
    end
end

results.time_series = time_series;
results.analysis_parameters = analysis_parameters;

%% Summary statistics
fprintf('\nSUMMARY STATISTICS:\n');

% Assembly metrics
if ~isempty(time_series.time_points)
    initial_connectivity = time_series.connected_particles(1) / time_series.n_particles(1) * 100;
    final_connectivity = time_series.connected_particles(end) / time_series.n_particles(end) * 100;
    max_cluster_size = max(time_series.largest_cluster_size);
    final_cluster_size = time_series.largest_cluster_size(end);
    max_clusters = max(time_series.n_clusters);
    
    results.summary.initial_connectivity_percent = initial_connectivity;
    results.summary.final_connectivity_percent = final_connectivity;
    results.summary.connectivity_change = final_connectivity - initial_connectivity;
    results.summary.max_cluster_size = max_cluster_size;
    results.summary.final_cluster_size = final_cluster_size;
    results.summary.max_number_clusters = max_clusters;
    results.summary.final_number_clusters = time_series.n_clusters(end);
    
    % Assembly kinetics
    mid_point = round(length(time_series.time_points)/2);
    if length(time_series.connected_particles) > mid_point
        early_connectivity = mean(time_series.connected_particles(1:mid_point) ./ time_series.n_particles(1:mid_point));
        late_connectivity = mean(time_series.connected_particles(mid_point:end) ./ time_series.n_particles(mid_point:end));
        results.summary.assembly_rate = (late_connectivity - early_connectivity) / (time_series.time_points(end) - time_series.time_points(mid_point));
    else
        results.summary.assembly_rate = 0;
    end
    
    fprintf('  Connectivity: %.1f%% → %.1f%% (Δ = %.1f%%)\n', ...
        initial_connectivity, final_connectivity, final_connectivity - initial_connectivity);
    fprintf('  Max cluster size: %d particles\n', max_cluster_size);
    fprintf('  Final cluster size: %d particles\n', final_cluster_size);
    fprintf('  Max simultaneous clusters: %d\n', max_clusters);
    fprintf('  Final number of clusters: %d\n', time_series.n_clusters(end));
    fprintf('  Assembly rate: %.6f connectivity/time\n', results.summary.assembly_rate);
end

%% Assembly classification
if results.summary.final_connectivity_percent > 80 && results.summary.final_cluster_size > 20
    assembly_outcome = 'Successful_Crystal';
elseif results.summary.final_connectivity_percent > 50
    assembly_outcome = 'Partial_Assembly';
elseif results.summary.connectivity_change > 20
    assembly_outcome = 'Dynamic_Clustering';
else
    assembly_outcome = 'Swarm_Behavior';
end

results.summary.assembly_outcome = assembly_outcome;
fprintf('  Assembly outcome: %s\n', assembly_outcome);

%% Create summary plots
figure('Position', [100, 100, 1200, 800]);

subplot(2, 3, 1);
plot(time_series.time_points, time_series.n_particles, 'b-', 'LineWidth', 2);
xlabel('Time'); ylabel('Number of Particles'); title('Particle Count');
grid on;

subplot(2, 3, 2);
plot(time_series.time_points, time_series.connected_particles ./ time_series.n_particles * 100, 'r-', 'LineWidth', 2);
xlabel('Time'); ylabel('Connected Particles (%)'); title('Assembly Progress');
grid on; ylim([0 100]);

subplot(2, 3, 3);
plot(time_series.time_points, time_series.largest_cluster_size, 'g-', 'LineWidth', 2);
xlabel('Time'); ylabel('Largest Cluster Size'); title('Cluster Growth');
grid on;

subplot(2, 3, 4);
plot(time_series.time_points, time_series.n_clusters, 'c-', 'LineWidth', 2);
xlabel('Time'); ylabel('Number of Clusters'); title('Cluster Count');
grid on;

subplot(2, 3, 5);
plot(time_series.time_points, time_series.total_clustered_particles, 'm-', 'LineWidth', 2);
xlabel('Time'); ylabel('Total Clustered Particles'); title('Clustering Efficiency');
grid on;

subplot(2, 3, 6);
if ~isempty(time_series.cluster_areas) && any(time_series.cluster_areas > 0)
    plot(time_series.time_points, time_series.cluster_areas, 'k-', 'LineWidth', 2);
    xlabel('Time'); ylabel('Total Cluster Area'); title('Spatial Clustering');
    grid on;
else
    text(0.5, 0.5, 'No cluster area data', 'HorizontalAlignment', 'center');
    title('Spatial Clustering');
end

sgtitle(sprintf('Assembly Analysis: %s (%s)', exp_info.power, assembly_outcome));

%% Save results
save_name = sprintf('%s_analysis_results.mat', file_name(1:end-4));
save(save_name, 'results');

% Create a text summary for easy sharing
summary_file = sprintf('%s_summary.txt', file_name(1:end-4));
fid = fopen(summary_file, 'w');

fprintf(fid, '=== TRACKING ANALYSIS SUMMARY ===\n');
fprintf(fid, 'File: %s\n', file_name);
fprintf(fid, 'Magnification: %s\n', exp_info.magnification);
fprintf(fid, 'Power: %s\n', exp_info.power);
fprintf(fid, 'Trial: %s\n', exp_info.trial);
fprintf(fid, '\nBasic Stats:\n');
fprintf(fid, 'Trajectories: %d\n', results.basic_stats.total_tracks);
fprintf(fid, 'Time Range: %d - %d (%d points)\n', results.basic_stats.time_range, results.basic_stats.total_timepoints);
fprintf(fid, 'Field of View: %.1f × %.1f\n', results.spatial.field_of_view);
fprintf(fid, '\nAssembly Metrics:\n');
fprintf(fid, 'Initial Connectivity: %.1f%%\n', results.summary.initial_connectivity_percent);
fprintf(fid, 'Final Connectivity: %.1f%%\n', results.summary.final_connectivity_percent);
fprintf(fid, 'Max Cluster Size: %d\n', results.summary.max_cluster_size);
fprintf(fid, 'Final Cluster Size: %d\n', results.summary.final_cluster_size);
fprintf(fid, 'Assembly Outcome: %s\n', results.summary.assembly_outcome);
fprintf(fid, 'Assembly Rate: %.6f\n', results.summary.assembly_rate);

fclose(fid);

fprintf('\nRESULTS SAVED:\n');
fprintf('  Full results: %s\n', save_name);
fprintf('  Text summary: %s\n', summary_file);
fprintf('  Plots: Currently displayed\n');

%% Display shareable summary
fprintf('\n=== SHAREABLE SUMMARY ===\n');
fprintf('Experiment: %s | %s | %s\n', exp_info.magnification, exp_info.power, exp_info.trial);
fprintf('Particles: %d | Duration: %d | Field: %.0f×%.0f\n', ...
    results.basic_stats.total_tracks, results.basic_stats.time_duration, results.spatial.field_of_view);
fprintf('Assembly: %.1f%% → %.1f%% | Max cluster: %d | Outcome: %s\n', ...
    results.summary.initial_connectivity_percent, results.summary.final_connectivity_percent, ...
    results.summary.max_cluster_size, results.summary.assembly_outcome);
fprintf('Rate: %.6f | Final clusters: %d\n', results.summary.assembly_rate, results.summary.final_number_clusters);
fprintf('========================\n');

% end

%% Helper function to parse filename
function info = parse_filename(filename)
    info.filename = filename;
    
    % Extract magnification
    if contains(filename, '20x')
        info.magnification = '20x';
    elseif contains(filename, '60x')
        info.magnification = '60x';
    else
        info.magnification = 'unknown';
    end
    
    % Extract power
    pw_match = regexp(filename, 'PW(\d+)', 'tokens');
    if ~isempty(pw_match)
        info.power = ['PW' pw_match{1}{1}];
    else
        info.power = 'unknown';
    end
    
    % Extract trial/spot number
    spots_match = regexp(filename, 'spots[_]*(\d+)', 'tokens');
    if ~isempty(spots_match)
        info.trial = ['spots_' spots_match{1}{1}];
    else
        info.trial = 'unknown';
    end
end