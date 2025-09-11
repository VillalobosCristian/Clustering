clear all; clc; close all;

%% Load data
[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
data = readtable(fullfile(path_name, file_name));
data = data(4:end,:);
if iscell(data.TRACK_ID), data.TRACK_ID = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X), data.POSITION_X = cellfun(@str2double, data.POSITION_X); end
if iscell(data.POSITION_Y), data.POSITION_Y = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T), data.POSITION_T = cellfun(@str2double, data.POSITION_T); end
sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

%% Extract trajectories
X = {}; Y = {}; T = {};
count = 0;
for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    if height(track_data) >= 100
        count = count + 1;
        X{count} = track_data.POSITION_X;
        Y{count} = track_data.POSITION_Y;
        T{count} = track_data.POSITION_T;
    end
end

%% Find time range and define three time points
all_times = [];
for i = 1:length(T)
    all_times = [all_times; T{i}];
end
t_min = min(all_times);
t_max = max(all_times);
t_mid = round((t_min + t_max) / 2);

target_times = [t_min, t_mid, t_max];
time_labels = {'Beginning', 'Middle', 'End'};

fprintf('Time analysis:\n');
fprintf('Start time: %d\n', t_min);
fprintf('Middle time: %d\n', t_mid);
fprintf('End time: %d\n', t_max);

%% Calculate g(r) for all three time points
dr = 0.3;
r_max = 200;
r_values = dr:dr:r_max;
g_r_all = zeros(length(target_times), length(r_values));
positions_all = cell(3, 1);
N_all = zeros(3, 1);

for t_idx = 1:3
    target_time = target_times(t_idx);
    positions = [];
    
    % Get positions at target time
    for i = 1:length(X)
        [~, timeIdx] = min(abs(T{i} - target_time)); % Find closest time point
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
        end
    end
    
    positions_all{t_idx} = positions;
    N = size(positions, 1);
    N_all(t_idx) = N;
    
    if N > 1
        area = (max(positions(:,1)) - min(positions(:,1))) * (max(positions(:,2)) - min(positions(:,2)));
        rho = N / area;
        
        % Calculate g(r)
        g_r = zeros(size(r_values));
        for r_idx = 1:length(r_values)
            r = r_values(r_idx);
            total_count = 0;
            for i = 1:N
                for j = 1:N
                    if i ~= j
                        distance = norm(positions(j, :) - positions(i, :));
                        if distance >= r && distance < r + dr
                            total_count = total_count + 1;
                        end
                    end
                end
            end
            average_count = total_count / N;
            shell_area = 2 * pi * r * dr;
            g_r(r_idx) = average_count / shell_area / rho;
        end
        g_r_all(t_idx, :) = g_r;
    end
    
    fprintf('%s (t=%d): %d particles\n', time_labels{t_idx}, target_time, N);
end

%% Plot results
figure('Position', [100, 100, 1400, 800], 'Color', 'white');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

% Colors for the three time points
colors = {'blue', 'green', 'red'};
line_styles = {'-', '--', ':'};

% Top row: Particle positions at three times
for t_idx = 1:3
    subplot(2, 3, t_idx);
    positions = positions_all{t_idx};
    scatter(positions(:,1), positions(:,2), 50, colors{t_idx}, 'filled', 'MarkerEdgeColor', 'k');
    axis equal;
    xlabel('$x$ position ($\mu$m)', 'FontSize', 12);
    ylabel('$y$ position ($\mu$m)', 'FontSize', 12);
    title(sprintf('%s ($t = %d$)\n$N = %d$ particles', time_labels{t_idx}, target_times(t_idx), N_all(t_idx)), 'FontSize', 12);
    set(gca, 'FontSize', 10, 'LineWidth', 1.2);
    box on;
end

% Bottom row: g(r) comparison
subplot(2, 3, [4, 5, 6]);
hold on;
for t_idx = 1:3
    plot(r_values, g_r_all(t_idx, :), 'Color', colors{t_idx}, 'LineStyle', line_styles{t_idx}, 'LineWidth', 2, ...
         'DisplayName', sprintf('%s (t=%d)', time_labels{t_idx}, target_times(t_idx)));
end
plot([0, r_max], [1, 1], 'k--', 'LineWidth', 1, 'DisplayName', 'Random ($g(r) = 1$)');
xlim([0, r_max]);
ylim([0, max(g_r_all(:))*1.1]);
xlabel('Distance $r$ ($\mu$m)', 'FontSize', 14);
ylabel('$g(r)$', 'FontSize', 14);
title('Radial Distribution Function Evolution', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
box on;
grid on;
set(gca, 'GridAlpha', 0.3);

% Set figure size
set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 12, 8]);

baseDir = pwd;
fileName = 'gDR_time_evolution';
% savefigures(baseDir, fileName);