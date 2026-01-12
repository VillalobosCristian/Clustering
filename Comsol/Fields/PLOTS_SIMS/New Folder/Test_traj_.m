clear all;
clc;
close all

[file_name, path_name] = uigetfile('*.csv', 'Select a CSV file');

if isequal(file_name, 0)
    return;
end

data = readtable(fullfile(path_name, file_name));
data = data(4:end,:);

if iscell(data.TRACK_ID)
    data.TRACK_ID = cellfun(@str2double, data.TRACK_ID);
end
if iscell(data.POSITION_X)
    data.POSITION_X = cellfun(@str2double, data.POSITION_X);
end
if iscell(data.POSITION_Y)
    data.POSITION_Y = cellfun(@str2double, data.POSITION_Y);
end
if iscell(data.POSITION_T)
    data.POSITION_T = cellfun(@str2double, data.POSITION_T) / 24;
else
    data.POSITION_T = data.POSITION_T / 24;
end

sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);
num_tracks = length(unique_tracks);

X = cell(num_tracks, 1);
Y = cell(num_tracks, 1);
T = cell(num_tracks, 1);
mean_speeds = zeros(num_tracks, 1);

valid_tracks = 0;
for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    
    if height(track_data) < 1
        continue;
    end
    
    valid_tracks = valid_tracks + 1;
    
    X{i} = track_data.POSITION_X;
    Y{i} = track_data.POSITION_Y;
    T{i} = track_data.POSITION_T;
    
    dx = diff(X{i});
    dy = diff(Y{i});
    dt = diff(T{i});
    speeds = sqrt((dx./dt).^2 + (dy./dt).^2);
    
    valid_speeds = speeds(isfinite(speeds));
    if ~isempty(valid_speeds)
        mean_speeds(i) = mean(valid_speeds);
    else
        mean_speeds(i) = NaN;
    end
end

valid_mask = ~cellfun('isempty', X);
X = X(valid_mask);
Y = Y(valid_mask);
T = T(valid_mask);
mean_speeds = mean_speeds(valid_mask);

valid_speed_mask = isfinite(mean_speeds);
X = X(valid_speed_mask);
Y = Y(valid_speed_mask);
T = T(valid_speed_mask);
mean_speeds = mean_speeds(valid_speed_mask);

if isempty(X) || isempty(Y)
    error('No valid trajectories found after filtering.');
end

min_speed = min(mean_speeds);
max_speed = max(mean_speeds);

x_min = min(cellfun(@min, X));
x_max = max(cellfun(@max, X));
y_min = min(cellfun(@min, Y));
y_max = max(cellfun(@max, Y));
margin = 0.05 * max(x_max-x_min, y_max-y_min);

figure('Position', [100, 100, 900, 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultColorbarTickLabelInterpreter', 'latex');

hold on;

cmap = thermal(256);

for i = 1:length(X)
    norm_speed = (mean_speeds(i) - min_speed) / (max_speed - min_speed);
    color_idx = max(1, min(256, round(norm_speed * 255) + 1));
    trajectory_color = cmap(color_idx, :);
    plot(X{i}, Y{i}, '-', 'Color', trajectory_color, 'LineWidth', 1);
end

xlabel('$x\ (\mu m)$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$y\ (\mu m)$', 'FontSize', 14, 'Interpreter', 'latex');

grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);
axis([x_min-margin, x_max+margin, y_min-margin, y_max+margin]);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02 0.02]);
set(gca, 'FontSize', 12);

colormap(cmap);
cb = colorbar;
cb.Label.String = '$\mathrm{Mean\ Speed\ (\mu m/s)}$';
cb.Label.FontSize = 14;
cb.Label.Interpreter = 'latex';
cb.LineWidth = 1.5;
cb.TickDirection = 'out';
cb.Box = 'off';
cb.Ticks = [0, 0.5, 1];
cb.TickLabels = {sprintf('%.2f', min_speed), sprintf('%.2f', (min_speed+max_speed)/2), sprintf('%.2f', max_speed)};

box on;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

valid_speeds = mean_speeds(isfinite(mean_speeds) & mean_speeds > 0);

figure('Position', [100, 100, 900, 700], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

num_bins = round(sqrt(length(valid_speeds)));
[counts, bin_edges] = histcounts(valid_speeds, num_bins);
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
bin_width = bin_edges(2) - bin_edges(1);
pdf_values = counts / (length(valid_speeds) * bin_width);

bar(bin_centers, pdf_values, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'k');

xlabel('$\mathrm{Velocity\ (\mu m/s)}$', 'FontSize', 14);
ylabel('$\mathrm{Probability\ Density}$', 'FontSize', 14);

grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.02 0.02]);
set(gca, 'FontSize', 16);
box on;

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);