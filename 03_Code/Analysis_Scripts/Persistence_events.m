clear all; clc; close all;

[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), return; end

data = readtable(fullfile(path_name, file_name));
data = data(4:end,:); % Remove TrackMate headers
% data = data(4:end,:); % Remove TrackMate headers
% data = data(4:end,:);
% Convert frame numbers to seconds
fps = 5;
% fprintf('Converting POSITION_T from frames to seconds (fps = %d)...\n', fps);
data.POSITION_T = data.POSITION_T / fps;
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

%% Parameters 
calib = 1.0;              % scale
distance_threshold = 5.0;  % distance en micro
time_step = 1;            % 

%% Get all unique time points
all_times = [];
for i = 1:length(T)
    all_times = [all_times; T{i}];
end
unique_times = unique(all_times);
analysis_times = unique_times(1:time_step:end);

%%  Create particle history storage
particle_history = cell(length(X), 1);  % One cell per particle
for i = 1:length(X)
    % Each particle gets an array: 1 = connected, 0 = not connected, -1 = not present
    particle_history{i} = -ones(length(analysis_times), 1);  % Start as "not present"
end

connected_particles = [];
total_particles = [];
connection_ratio = [];
time_points = [];

%  Arrays to store persistence results (we'll fill these later)
stable_particles = [];
transient_particles = [];

%%
for t_idx = 1:length(analysis_times)
    current_time = analysis_times(t_idx);
    
    % Extract positions at current time 
    positions = [];
    present_particles = []; % track which particles are present
    
    for i = 1:length(X)
        timeIdx = find(T{i} == current_time, 1);
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
            present_particles = [present_particles; i];  % Remember which particle this is
        end
    end
    
    positions = positions * calib;
    N_particles = size(positions, 1);
    
    if N_particles < 2
        continue; % Skip if not enough particles
    end
    
    % Calculate pairwise distances (pdist already computes Euclidean distances)
    distances = pdist(positions);
    distance_matrix = squareform(distances);  
    
    %  Check connections and record for each individual particle
    connected_count = 0;
    for j = 1:N_particles
        particle_idx = present_particles(j);  % Which particle is this?
        
        % Check if this particle has any neighbor within threshold
        neighbors = distance_matrix(j, :);
        neighbors(j) = inf; % Exclude self
        
        if any(neighbors <= distance_threshold)
            connected_count = connected_count + 1;
            % Record that this particle was connected at this time
            particle_history{particle_idx}(t_idx) = 1;
        else
            %Record that this particle was present but not connected
            particle_history{particle_idx}(t_idx) = 0;
        end
    end
    
    % Store results
    time_points = [time_points; current_time];
    connected_particles = [connected_particles; connected_count];
    total_particles = [total_particles; N_particles];
    connection_ratio = [connection_ratio; connected_count / N_particles];
    
end

%%

% Look at the first few particles
for i = 1:min(3, length(particle_history))
    history = particle_history{i};
    n_present = sum(history >= 0);  % Times when particle was present
    n_connected = sum(history == 1);  % Times when particle was connected
    n_disconnected = sum(history == 0);  % Times when particle was present but disconnected
    
    fprintf('  Particle %d: present %d times, connected %d times (%.1f%%)\n', ...
        i, n_present, n_connected, 100*n_connected/max(n_present,1));
end

%% 
figure('Position', [100, 100, 1200, 400], 'Color', 'w');

% Plot 1: Number of connected particles over time
subplot(1, 2, 1);
plot(time_points, connected_particles, 'b-', 'LineWidth', 2);
hold on;
plot(time_points, total_particles, 'k--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
xlabel('Time', 'FontSize', 12);
ylabel('Number of Particles', 'FontSize', 12);
legend({'Connected Particles', 'Total Particles'}, 'Location', 'best');
title('Current Analysis (No Persistence Yet)');
grid on; box on;

% Plot 2: Connection ratio over time
subplot(1, 2, 2);
plot(time_points, connection_ratio * 100, 'r-', 'LineWidth', 2);
xlabel('Time', 'FontSize', 12);
ylabel('Connected Particles (%)', 'FontSize', 12);
ylim([0 100]);
grid on; box on;
%%
% For each connected particle:
%1. Look backwards from current time
%2. Count consecutive connected frames 
%3. Allow some gap tolerance (1 frame disconnection is OK)
%4. If >= 5 consecutive frames = "stable"
%5. If < 5 consecutive frames = "transient"
persistence_frames = 5;     % Minimum frames to be considered "stable"
tolerance_frames = 1;       % Allowed gap frames in connection sequence
stable_particles_over_time = [];
transient_particles_over_time = [];

for t_idx = 1:length(analysis_times)
    current_time = analysis_times(t_idx);
    
    % Find which particles are present at this time
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
        continue;
    end
    
    % For each present particle, check if it's stably connected
    stable_count = 0;
    transient_count = 0;
    
    for j = 1:length(present_particles)
        particle_idx = present_particles(j);
        
        % Is this particle currently connected?
        if particle_history{particle_idx}(t_idx) == 1 %yes
            
            % Check persistence: count actual connected frames with tolerance
            connected_count = 0;
            current_gap = 0;
            
            % Look backwards from current time
            for k = t_idx:-1:1
                if particle_history{particle_idx}(k) == 1
                    % Particle was connected
                    connected_count = connected_count + 1;
                    current_gap = 0;  % Reset gap counter
                elseif particle_history{particle_idx}(k) == 0
                    % Particle was present but not connected - this is a gap
                    current_gap = current_gap + 1;
                    if current_gap > tolerance_frames
                        % Too many consecutive gaps - break the sequence
                        break;
                    end
                elseif particle_history{particle_idx}(k) == -1
                    % Particle was not present - this is also a gap
                    current_gap = current_gap + 1;
                    if current_gap > tolerance_frames
                        break;
                    end
                end
            end
            
            % Classify based on persistence
            if connected_count >= persistence_frames
                stable_count = stable_count + 1;
            else
                transient_count = transient_count + 1;
            end
        end
    end
    
    stable_particles_over_time = [stable_particles_over_time; stable_count];
    transient_particles_over_time = [transient_particles_over_time; transient_count];

end
%%
stable_ratio = stable_particles_over_time ./ max(total_particles, 1);
transient_ratio = transient_particles_over_time ./ max(total_particles, 1);

figure('Position', [100, 100, 1600, 800], 'Color', 'w');

% Plot 1: Particle counts over time
subplot(2, 3, 1);
plot(time_points, total_particles, 'k--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
hold on;
plot(time_points, connected_particles, 'b-', 'LineWidth', 2);
plot(time_points, stable_particles_over_time, 'r-', 'LineWidth', 2);
plot(time_points, transient_particles_over_time, 'g--', 'LineWidth', 1.5);
xlabel('Time', 'FontSize', 12);
ylabel('Number of Particles', 'FontSize', 12);
legend({'Total', 'Connected', 'Stable', 'Transient'}, 'Location', 'best');
title('Particle Counts with Persistence');
grid on; box on;

% Plot 2: Connection percentages
subplot(2, 3, 2);
plot(time_points, connection_ratio * 100, 'b-', 'LineWidth', 2);
hold on;
plot(time_points, stable_ratio * 100, 'r-', 'LineWidth', 2);
plot(time_points, transient_ratio * 100, 'g--', 'LineWidth', 1.5);
xlabel('Time', 'FontSize', 12);
ylabel('Percentage (%)', 'FontSize', 12);
legend({'Total Connected', 'Stable', 'Transient'}, 'Location', 'best');
title('Connection Percentages');
ylim([0 100]);
grid on; box on;

% Plot 3: Stability ratio
subplot(2, 3, 3);
stability_index = stable_particles_over_time ./ max(connected_particles, 1);
plot(time_points, stability_index * 100, 'm', 'LineWidth', 2);
xlabel('Time', 'FontSize', 12);
ylabel('Stability Index (%)', 'FontSize', 12);
title('Stability Index (Stable/Connected)');
ylim([0 100]);
grid on; box on;

%%
% Save results - NOW INCLUDING particle_history and analysis_times
results_file = [file_name(1:end-4), '_persistence_results.mat'];
save(fullfile(path_name, results_file), 'time_points', 'stable_particles_over_time', ...
     'transient_particles_over_time', 'connection_ratio', 'stable_ratio', 'transient_ratio', ...
     'persistence_frames', 'tolerance_frames', 'distance_threshold', ...
     'particle_history', 'analysis_times');

fprintf('\nResults saved to: %s\n', results_file);
fprintf('particle_history has been saved for animation!\n');