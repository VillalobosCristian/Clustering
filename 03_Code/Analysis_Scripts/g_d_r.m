clear all; clc; close all;

%% PARAMETERS
PARTICLE_DIAMETER = 1.45;
DR = 0.2;
R_MAX = 20;
MIN_TRACK_LENGTH = 10;
NUM_TIME_SLICES = 3;

%% LOAD DATA
[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), return; end

data = readtable(fullfile(path_name, file_name));
data = data(4:end,:);

fps = 1;
data.POSITION_T = data.POSITION_T / fps;

if iscell(data.TRACK_ID), data.TRACK_ID = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X), data.POSITION_X = cellfun(@str2double, data.POSITION_X); end
if iscell(data.POSITION_Y), data.POSITION_Y = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T), data.POSITION_T = cellfun(@str2double, data.POSITION_T); end

sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

%% EXTRACT TRAJECTORIES
X = {}; Y = {}; T = {};
count = 0;

for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    if height(track_data) >= MIN_TRACK_LENGTH
        count = count + 1;
        X{count} = track_data.POSITION_X;
        Y{count} = track_data.POSITION_Y;
        T{count} = track_data.POSITION_T;
    end
end

fprintf('Loaded %d valid trajectories\n', count);

%% DETERMINE TIME POINTS
all_times = [];
for i = 1:length(T)
    all_times = [all_times; T{i}];
end

t_min = min(all_times);
t_max = max(all_times);
target_times = linspace(t_min, t_max, NUM_TIME_SLICES);

%% COMPUTE g(r) AT EACH TIME
r_values = DR:DR:R_MAX;
g_r_all = zeros(NUM_TIME_SLICES, length(r_values));
positions_all = cell(NUM_TIME_SLICES, 1);
N_all = zeros(NUM_TIME_SLICES, 1);

for t_idx = 1:NUM_TIME_SLICES
    target_time = target_times(t_idx);
    positions = [];
    
    for i = 1:length(X)
        [~, timeIdx] = min(abs(T{i} - target_time));
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
        end
    end
    
    positions_all{t_idx} = positions;
    N = size(positions, 1);
    N_all(t_idx) = N;
    
    if N > 1
        area = (max(positions(:,1)) - min(positions(:,1))) * ...
               (max(positions(:,2)) - min(positions(:,2)));
        rho = N / area;
        
        distances = pdist(positions);
        
        for r_idx = 1:length(r_values)
            r = r_values(r_idx);
            count_in_shell = sum(distances >= r-DR/2 & distances < r+DR/2);
            average_count = (2 * count_in_shell) / N;
            shell_area = pi * ((r+DR/2)^2 - (r-DR/2)^2);
            g_r_all(t_idx, r_idx) = average_count / (shell_area * rho);
        end
    end
end

%% PLOT - IMPROVED VISUALS
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');

% Use better colormap
cmap = viridis(NUM_TIME_SLICES);

figure('Position', [50, 50, 1400, 800], 'Color', 'white');

% Top row: positions at each time
for t_idx = 1:NUM_TIME_SLICES
    subplot(2, NUM_TIME_SLICES, t_idx);
    positions = positions_all{t_idx};
    
    if ~isempty(positions)
        scatter(positions(:,1), positions(:,2), 10, cmap(t_idx,:), ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, ...
                'MarkerFaceAlpha', 0.8);
        axis equal;
        xlabel('$x\ (\mu\mathrm{m})$', 'FontSize', 14);
        ylabel('$y\ (\mu\mathrm{m})$', 'FontSize', 14);
        title(sprintf('$t = %.1f$ s ($N = %d$)', target_times(t_idx), N_all(t_idx)), ...
              'FontSize', 14);
        set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'TickDir', 'out');
        box on;
    end
end

% Bottom: g(r) comparison
subplot(2, NUM_TIME_SLICES, NUM_TIME_SLICES+1:2*NUM_TIME_SLICES);
hold on;

for t_idx = 1:NUM_TIME_SLICES
    if N_all(t_idx) > 1
        plot(r_values, g_r_all(t_idx, :), ...
             'Color', cmap(t_idx,:), ...
             'LineWidth', 3, ...
             'DisplayName', sprintf('$t = %.1f$ s', target_times(t_idx)));
    end
end

yline(1, 'k--', 'LineWidth', 2, 'DisplayName', 'Random');
xline(PARTICLE_DIAMETER, 'Color', [0.8 0.2 0.2], 'LineStyle', '--', ...
      'LineWidth', 2, 'DisplayName', sprintf('$d = %.2f$ $\\mu$m', PARTICLE_DIAMETER));

xlim([0, R_MAX]);
ylim([0, max(g_r_all(:))*1.1]);
xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$g(r)$', 'FontSize', 18);

legend('Location', 'northeast', 'FontSize', 14, 'Box', 'off');
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickDir', 'out');
box on;
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 12, 7]);


baseDir = pwd;
fileName = 'GDR';
savefigures(baseDir, fileName);