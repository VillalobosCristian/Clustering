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

%% Get positions at t=300
target_time =  5500;
positions = [];
for i = 1:length(X)
    timeIdx = find(T{i} == target_time, 1);
    if ~isempty(timeIdx)
        positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
    end
end

N = size(positions, 1);
area = (max(positions(:,1)) - min(positions(:,1))) * (max(positions(:,2)) - min(positions(:,2)));
rho = N / area;

%% Calculate g(r)
dr = 0.3;
r_max = 20;
r_values = dr:dr:r_max;
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

%% Plot

figure('Position', [100, 100, 1200, 500], 'Color', 'white');

% Set default interpreter to LaTeX
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
subplot(1,2,1);
scatter(positions(:,1), positions(:,2), 100, 'b', 'filled','MarkerEdgeColor','k');
axis equal;
xlabel('$x$ position ($\mu$m)', 'FontSize', 14);
ylabel('$y$ position ($\mu$m)', 'FontSize', 14);
title(sprintf('Particle Positions at $t = %d$\n$N = %d$ particles', target_time, N), 'FontSize', 14);
set(gca, 'FontSize', 16, 'LineWidth', 1.2);
box on

subplot(1,2,2);
plot(r_values, g_r, 'k-', 'LineWidth', 2);
hold on;
plot([0, r_max], [1, 1], 'r--');
xlim([0, r_max]);
ylim([0, max(g_r)*1.1]);
xlabel('Distance $r$ ($\mu$m)', 'FontSize', 14);
ylabel('$g(r)$', 'FontSize', 14);set(gca, 'FontSize', 12, 'LineWidth', 1.2);
box on
set(gca, 'FontSize', 16, 'LineWidth', 1.2);
% Set figure size
set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), 8, 6]);
baseDir = pwd;
fileName = 'gDR';
savefigures(baseDir, fileName);