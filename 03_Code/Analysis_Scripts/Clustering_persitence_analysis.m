clear all; clc; close all;

%% ------------------------- Load CSV data -------------------------------
[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), return; end

data = readtable(fullfile(path_name, file_name));
data = data(4:end,:); % Remove TrackMate headers

% Convert columns to numeric if needed
if iscell(data.TRACK_ID),    data.TRACK_ID    = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X),  data.POSITION_X  = cellfun(@str2double, data.POSITION_X); end
if iscell(data.POSITION_Y),  data.POSITION_Y  = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T),  data.POSITION_T  = cellfun(@str2double, data.POSITION_T); end

% Sort for easier processing
sorted_data   = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

%% ------------------------- Extract trajectories ------------------------
min_track_length = 10;
X = {}; Y = {}; T = {}; particle_ids = [];
count = 0;

for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    if height(track_data) >= min_track_length
        count = count + 1;
        X{count} = track_data.POSITION_X;
        Y{count} = track_data.POSITION_Y;
        T{count} = track_data.POSITION_T;
        particle_ids(count) = unique_tracks(i);
    end
end

nTracks = numel(X);
if nTracks == 0
    error('No tracks with at least %d points were found.', min_track_length);
end

fprintf('Loaded %d tracks with ≥%d points\n', nTracks, min_track_length);

%% ------------------------- Time grid for analysis ----------------------
% Collect all time points present in kept tracks
all_times = [];
for i = 1:nTracks, all_times = [all_times; T{i}(:)]; end
unique_times = unique(all_times);

time_step = 1;  % analyze every frame
analysis_times = unique_times(1:time_step:end);
nT = numel(analysis_times);

fprintf('Analyzing %d time points (time step: %d)\n', nT, time_step);

%% ------------------------- Parameters ----------------------------------
calib = 1.0;                 % px -> µm
distance_threshold = 5;      % µm
persistence_frames = 5;      % min frames to be "stable"
tolerance_frames   = 1;      % allowed consecutive gaps

fprintf('\n=== ANALYSIS PARAMETERS ===\n');
fprintf('Distance threshold: %.1f μm\n', distance_threshold);
fprintf('Persistence threshold: %d frames\n', persistence_frames);
fprintf('Gap tolerance: %d frames\n', tolerance_frames);

%% ------------------------- ORIGINAL: Build connection history (slow) ---
fprintf('Building particle connection history...\n');
tic;

% ORIGINAL METHOD: Cell array of individual particle histories
particle_history = cell(nTracks, 1);
for i = 1:nTracks
    particle_history{i} = zeros(nT, 1); % 0 = not connected, 1 = connected
end

% ORIGINAL METHOD: Per-frame find() operations (SLOW!)
for t_idx = 1:nT
    current_time = analysis_times(t_idx);
    
    % ORIGINAL: Extract positions using repeated find() calls
    positions = [];
    present_particles = [];
    for i = 1:nTracks  % This loops through ALL tracks every frame
        timeIdx = find(T{i} == current_time, 1); % EXPENSIVE find() operation
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
            present_particles = [present_particles; i];
        end
    end
    
    positions = positions * calib;
    N_particles = size(positions, 1);
    
    if N_particles < 2
        % Mark all present particles as disconnected
        for j = 1:length(present_particles)
            particle_idx = present_particles(j);
            particle_history{particle_idx}(t_idx) = 0;
        end
        continue;
    end
    
    % ORIGINAL: Always use pdist/squareform (not optimized)
    distances = pdist(positions);
    distance_matrix = squareform(distances); % Creates full NxN matrix
    
    % Check connections for each present particle
    for j = 1:N_particles
        particle_idx = present_particles(j);
        
        % Check if particle j has any neighbor within threshold
        neighbors = distance_matrix(j, :);
        neighbors(j) = inf; % Exclude self
        
        if any(neighbors <= distance_threshold)
            particle_history{particle_idx}(t_idx) = 1; % Connected
        else
            particle_history{particle_idx}(t_idx) = 0; % Not connected
        end
    end
    
    % Progress indicator for slow operation
    if mod(t_idx, 100) == 0
        fprintf('  Processed %d/%d time points (%.1f%%)\n', t_idx, nT, 100*t_idx/nT);
    end
end

history_time = toc;
fprintf('Connection history built in %.2f seconds\n', history_time);

%% ------------------------- ORIGINAL: Persistence analysis (slow) -------
fprintf('Computing persistence classifications...\n');
tic;

% Storage for results
time_points = analysis_times;
total_particles = zeros(nT, 1);
connected_particles = zeros(nT, 1);
stable_clusters = zeros(nT, 1);
transient_particles = zeros(nT, 1);
disconnected_particles = zeros(nT, 1);

% ORIGINAL METHOD: Per-frame backtracking (VERY SLOW!)
for t_idx = 1:nT
    current_time = analysis_times(t_idx);
    
    % ORIGINAL: Count particles present using find() again
    present_particles = [];
    for i = 1:nTracks
        timeIdx = find(T{i} == current_time, 1); % Another expensive find()
        if ~isempty(timeIdx)
            present_particles = [present_particles; i];
        end
    end
    
    N_particles = length(present_particles);
    if N_particles < 2
        total_particles(t_idx) = N_particles;
        connected_particles(t_idx) = 0;
        stable_clusters(t_idx) = 0;
        transient_particles(t_idx) = 0;
        disconnected_particles(t_idx) = N_particles;
        continue;
    end
    
    % ORIGINAL: Per-particle persistence checking with backtracking
    stable_count = 0;
    transient_count = 0;
    total_connected = 0;
    disconnected_count = 0;
    
    for j = 1:length(present_particles)
        particle_idx = present_particles(j);
        
        % Check if currently connected
        if particle_history{particle_idx}(t_idx) == 1
            total_connected = total_connected + 1;
            
            % ORIGINAL: Backtrack through time (SLOW!)
            persistent_frames = 0;
            gap_tolerance = 0;
            
            % Look backwards from current time
            for k = t_idx:-1:1  % This is O(time) per particle per frame!
                if particle_history{particle_idx}(k) == 1
                    persistent_frames = persistent_frames + 1;
                    gap_tolerance = 0; % Reset gap counter
                elseif gap_tolerance < tolerance_frames
                    gap_tolerance = gap_tolerance + 1;
                    % Don't increment persistent_frames but don't break yet
                else
                    break; % Too many consecutive gaps
                end
            end
            
            % Classify as stable or transient
            if persistent_frames >= persistence_frames
                stable_count = stable_count + 1;
            else
                transient_count = transient_count + 1;
            end
        else
            disconnected_count = disconnected_count + 1;
        end
    end
    
    % Store results
    total_particles(t_idx) = N_particles;
    connected_particles(t_idx) = total_connected;
    stable_clusters(t_idx) = stable_count;
    transient_particles(t_idx) = transient_count;
    disconnected_particles(t_idx) = disconnected_count;
    
    % Progress indicator for slow operation
    if mod(t_idx, 100) == 0
        fprintf('  Analyzed %d/%d time points (%.1f%%)\n', t_idx, nT, 100*t_idx/nT);
    end
end

persistence_time = toc;
fprintf('Persistence analysis completed in %.2f seconds\n', persistence_time);

%% ------------------------- Calculate derived metrics -------------------
connection_ratio = connected_particles ./ max(total_particles, 1);
stable_ratio = stable_clusters ./ max(total_particles, 1);
stability_index = stable_clusters ./ max(connected_particles, 1);

%% ------------------------- Statistical summary -------------------------
fprintf('\n=== ANALYSIS RESULTS ===\n');

% Overall statistics
total_timepoints = length(time_points);
mean_particles = mean(total_particles);
mean_connected = mean(connected_particles);
mean_stable = mean(stable_clusters);

fprintf('Dataset overview:\n');
fprintf('  Total timepoints: %d\n', total_timepoints);
fprintf('  Mean particles per frame: %.1f\n', mean_particles);
fprintf('  Mean connected particles: %.1f (%.1f%%)\n', mean_connected, 100*mean_connected/mean_particles);
fprintf('  Mean stable clusters: %.1f (%.1f%%)\n', mean_stable, 100*mean_stable/mean_particles);

% Assembly dynamics
if length(connection_ratio) > 10
    initial_conn = mean(connection_ratio(1:min(3, end)));
    final_conn = mean(connection_ratio(max(1, end-2):end));
    initial_stable = mean(stable_ratio(1:min(3, end)));
    final_stable = mean(stable_ratio(max(1, end-2):end));
    
    fprintf('\nAssembly dynamics:\n');
    fprintf('  Initial connectivity: %.1f%% → Final: %.1f%% (Δ = %+.1f%%)\n', ...
        initial_conn*100, final_conn*100, (final_conn-initial_conn)*100);
    fprintf('  Initial stable: %.1f%% → Final: %.1f%% (Δ = %+.1f%%)\n', ...
        initial_stable*100, final_stable*100, (final_stable-initial_stable)*100);
end

% Peak values
[max_connected, max_conn_idx] = max(connected_particles);
[max_stable, max_stable_idx] = max(stable_clusters);

fprintf('\nPeak values:\n');
fprintf('  Max connected particles: %d at time %.0f\n', max_connected, time_points(max_conn_idx));
fprintf('  Max stable clusters: %d at time %.0f\n', max_stable, time_points(max_stable_idx));

%% ------------------------- Create publication plots --------------------
fprintf('\nGenerating plots...\n');

figure('Position', [100, 100, 1600, 900], 'Color', 'white');

% Colors
blue_color = [0.2, 0.4, 0.8];
red_color = [0.9, 0.2, 0.2]; 
green_color = [0.1, 0.6, 0.1];
gray_color = [0.5, 0.5, 0.5];

% Plot 1: Absolute numbers over time
subplot(2, 3, 1);
plot(time_points, total_particles, 'k--', 'LineWidth', 1.5, 'Color', gray_color);
hold on;
plot(time_points, connected_particles, 'b-', 'LineWidth', 2, 'Color', blue_color);
plot(time_points, stable_clusters, 'r-', 'LineWidth', 2, 'Color', red_color);
plot(time_points, transient_particles, 'g--', 'LineWidth', 1.5, 'Color', green_color);
xlabel('Time', 'FontSize', 12);
ylabel('Number of Particles', 'FontSize', 12);
legend({'Total', 'Connected', 'Stable', 'Transient'}, 'Location', 'best');
title('Particle Counts Over Time', 'FontSize', 14);
grid on; box on;

% Plot 2: Percentages over time
subplot(2, 3, 2);
plot(time_points, connection_ratio * 100, 'b-', 'LineWidth', 2, 'Color', blue_color);
hold on;
plot(time_points, stable_ratio * 100, 'r-', 'LineWidth', 2, 'Color', red_color);
xlabel('Time', 'FontSize', 12);
ylabel('Percentage (%)', 'FontSize', 12);
legend({'Connected', 'Stable'}, 'Location', 'best');
title('Connection Percentages', 'FontSize', 14);
ylim([0, 100]);
grid on; box on;

% Plot 3: Assembly rates (derivatives)
subplot(2, 3, 3);
if length(connection_ratio) > 5
    smooth_conn = smoothdata(connection_ratio, 'gaussian', 5);
    conn_rate = gradient(smooth_conn, time_points);
    plot(time_points, conn_rate * 100, 'b-', 'LineWidth', 2, 'Color', blue_color);
    hold on;
    
    smooth_stable = smoothdata(stable_ratio, 'gaussian', 5);
    stable_rate = gradient(smooth_stable, time_points);
    plot(time_points, stable_rate * 100, 'r-', 'LineWidth', 2, 'Color', red_color);
    
    xlabel('Time', 'FontSize', 12);
    ylabel('Assembly Rate (%/time)', 'FontSize', 12);
    legend({'Connection Rate', 'Stability Rate'}, 'Location', 'best');
    title('Assembly Rates', 'FontSize', 14);
    grid on; box on;
else
    text(0.5, 0.5, 'Insufficient data for rate calculation', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    title('Assembly Rates', 'FontSize', 14);
end

% Plot 4: Stability index over time
subplot(2, 3, 4);
plot(time_points, stability_index * 100, 'g-', 'LineWidth', 2, 'Color', green_color);
xlabel('Time', 'FontSize', 12);
ylabel('Stability Index (%)', 'FontSize', 12);
title('Fraction of Connected that are Stable', 'FontSize', 14);
ylim([0, 100]);
grid on; box on;

% Plot 5: Histogram of particle states
subplot(2, 3, 5);
total_stable = sum(stable_clusters);
total_transient = sum(transient_particles);
total_disconnected = sum(disconnected_particles);
total_observations = total_stable + total_transient + total_disconnected;

bar_data = [total_stable, total_transient, total_disconnected] / total_observations * 100;
bar_colors = [red_color; green_color; gray_color];
b = bar(bar_data, 'FaceColor', 'flat');
b.CData = bar_colors;

set(gca, 'XTickLabel', {'Stable', 'Transient', 'Disconnected'});
ylabel('Percentage of Total Observations (%)', 'FontSize', 12);
title('Overall State Distribution', 'FontSize', 14);
grid on; box on;

% Plot 6: Time series with moving averages
subplot(2, 3, 6);
window = min(20, ceil(nT/10)); % Adaptive window size
ma_connected = movmean(connection_ratio, window);
ma_stable = movmean(stable_ratio, window);

plot(time_points, connection_ratio * 100, 'b-', 'LineWidth', 1, 'Color', [blue_color, 0.3]);
hold on;
plot(time_points, stable_ratio * 100, 'r-', 'LineWidth', 1, 'Color', [red_color, 0.3]);
plot(time_points, ma_connected * 100, 'b-', 'LineWidth', 3, 'Color', blue_color);
plot(time_points, ma_stable * 100, 'r-', 'LineWidth', 3, 'Color', red_color);

xlabel('Time', 'FontSize', 12);
ylabel('Percentage (%)', 'FontSize', 12);
legend({'Connected (raw)', 'Stable (raw)', sprintf('Connected (MA-%d)', window), sprintf('Stable (MA-%d)', window)}, ...
       'Location', 'best', 'FontSize', 10);
title('Trends with Moving Averages', 'FontSize', 14);
ylim([0, 100]);
grid on; box on;

%% ------------------------- Save results --------------------------------
% Save analysis results
results_filename = [file_name(1:end-4), '_persistence_analysis_original.mat'];
results_path = fullfile(path_name, results_filename);

save(results_path, 'time_points', 'total_particles', 'connected_particles', ...
     'stable_clusters', 'transient_particles', 'connection_ratio', 'stable_ratio', ...
     'stability_index', 'distance_threshold', 'persistence_frames', 'tolerance_frames', ...
     'calib', 'particle_history', 'analysis_times');

fprintf('\nResults saved to: %s\n', results_path);

% Save figure
fig_filename = [file_name(1:end-4), '_persistence_plots_original.png'];
fig_path = fullfile(path_name, fig_filename);
saveas(gcf, fig_path, 'png');
fprintf('Plots saved to: %s\n', fig_path);

%% ------------------------- Performance summary -------------------------
total_analysis_time = history_time + persistence_time;

fprintf('\n=== PERFORMANCE SUMMARY ===\n');
fprintf('Connection history: %.2f s\n', history_time);
fprintf('Persistence analysis: %.2f s\n', persistence_time);
fprintf('Total analysis time: %.2f s\n', total_analysis_time);
fprintf('Analysis rate: %.0f timepoints/second\n', nT / total_analysis_time);
fprintf('Used basic pdist/squareform for distance calculations (not optimized)\n');

fprintf('\nOriginal (slow) analysis complete!\n');