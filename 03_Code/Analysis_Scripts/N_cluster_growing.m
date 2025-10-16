clear all; clc; close all;

[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), return; end

data = readtable(fullfile(path_name, file_name));
data = data(4:end,:); % Remove TrackMate headers
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
time_step = 5;            % 

%% Get all unique time points
all_times = [];
for i = 1:length(T)
    all_times = [all_times; T{i}];
end
unique_times = unique(all_times);
analysis_times = unique_times(1:time_step:end);

%% Initialize storage arrays
connected_particles = [];
total_particles = [];
connection_ratio = [];
time_points = [];


%% Main analysis loop
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
    
    positions = positions * calib;
    N_particles = size(positions, 1);
    
    if N_particles < 2
        continue; % Skip if not enough particles
    end
    
    % Calculate pairwise distances
    distances = pdist(positions);  % Pairwise distances
    distance_matrix = squareform(distances);  % Convert to matrix form
    
    % Count particles with at least one neighbor within threshold
    connected_count = 0;
    for i = 1:N_particles
        % Check if particle i has any neighbor within threshold
        neighbors = distance_matrix(i, :);
        neighbors(i) = inf; % Exclude self
        
        if any(neighbors <= distance_threshold)
            connected_count = connected_count + 1;
        end
    end
    
    % Store results
    time_points = [time_points; current_time];
    connected_particles = [connected_particles; connected_count];
    total_particles = [total_particles; N_particles];
    connection_ratio = [connection_ratio; connected_count / N_particles];
    

end

%% Create publication-quality plot
figure('Position', [100, 100, 1200, 400], 'Color', 'w');


% figure('Position', [100, 100, 900, 700], 'Color', 'w');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

% Plot 1: Number of connected particles over time
subplot(1, 2, 1);
plot(time_points, connected_particles, 'b-', 'LineWidth', 2);
hold on;
plot(time_points, total_particles, 'k--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
xlabel('Time', 'FontSize', 16);
ylabel('Number of Particles', 'FontSize', 16);
legend({'Connected Particles', 'Total Particles'}, 'Location', 'best');
grid on; box on;

% Plot 2: Connection ratio over time
subplot(1, 2, 2);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');

plot(time_points, connection_ratio * 100, 'r-', 'LineWidth', 2);
xlabel('Time', 'FontSize', 16);
ylabel(' Percentage Connected Particles', 'FontSize', 16);
ylim([0 100]);
grid on; box on;

baseDir = pwd;
fileName = 'Cluster_grow';
savefigures(baseDir, fileName);
% subplot(1, 2, 3);
% if length(connection_ratio) > 5
%     % Calculate smoothed derivative
%     smooth_ratio = smoothdata(connection_ratio, 'gaussian', 5);
%     assembly_rate = gradient(smooth_ratio, time_points);
%     plot(time_points, assembly_rate, 'g-', 'LineWidth', 2);
%     xlabel('Time', 'FontSize', 12);
%     ylabel('Assembly Rate (% / time)', 'FontSize', 12);
%     grid on; box on;
% else
%     text(0.5, 0.5, 'Insufficient data\nfor rate calculation', ...
%         'HorizontalAlignment', 'center', 'FontSize', 12);
% end


% 
% % Find assembly phases
% if length(connection_ratio) > 10
%     initial_ratio = mean(connection_ratio(1:3));
%     final_ratio = mean(connection_ratio(end-2:end));
%     fprintf('Assembly efficiency: %.1f%% → %.1f%% (Δ = %.1f%%)\n', ...
%         initial_ratio*100, final_ratio*100, (final_ratio-initial_ratio)*100);
% end

