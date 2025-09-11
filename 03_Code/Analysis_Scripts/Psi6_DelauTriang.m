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

targetTime = 1;  
calib = 1.0;      

positions = [];
for i = 1:length(X)
    timeIdx = find(T{i} == targetTime, 1);
    if ~isempty(timeIdx)
        positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
    end
end
positions = positions * calib;

%% 
DT = delaunayTriangulation(positions(:,1), positions(:,2));
nPoints = size(DT.Points, 1);
Q6 = zeros(nPoints, 1);

for k = 1:nPoints
    orderIndex6 = 0;
    neighCount = 0;
    
    for j = 1:nPoints
        if k ~= j && isConnected(DT, k, j)
            alpha = atan2(DT.Points(j,2) - DT.Points(k,2), DT.Points(j,1) - DT.Points(k,1));
            orderIndex6 = orderIndex6 + exp(6*1i*alpha);
            neighCount = neighCount + 1;
        end
    end
    
    if neighCount > 0
        Q6(k) = abs(orderIndex6 / neighCount);
    end
end

%% 

    figure;
    scatter(positions(:,1), positions(:,2), 50, Q6, 'filled');
        xMin = min(positions(:,1)); xMax = max(positions(:,1));
    yMin = min(positions(:,2)); yMax = max(positions(:,2));
    padding = 0.1 * max(xMax - xMin, yMax - yMin);
    xLimits = [xMin - padding, xMax + padding];
    yLimits = [yMin - padding, yMax + padding];
    axis equal;
    xlim(xLimits);
    ylim(yLimits);
    hold on;
    [vx, vy] = voronoi(positions(:,1), positions(:,2));
    plot(vx, vy, 'k-', 'LineWidth', 0.5);
    hold off;
    
    colormap(internet); colorbar;

    %%
    figure('Position', [100, 100, 800, 600], 'Color', 'white');

xMin = min(positions(:,1)); xMax = max(positions(:,1));
yMin = min(positions(:,2)); yMax = max(positions(:,2));
padding = 0.08 * max(xMax - xMin, yMax - yMin);
xLimits = [xMin - padding, xMax + padding];
yLimits = [yMin - padding, yMax + padding];

scatter(positions(:,1), positions(:,2), 60, Q6, 'filled', ...
        'MarkerEdgeColor', [0.2, 0.2, 0.2], 'LineWidth', 0.8, ...
        'MarkerFaceAlpha', 1);

% Axes setup
axis equal;
xlim(xLimits);
ylim(yLimits);

% Voronoi cells
hold on;
[vx, vy] = voronoi(positions(:,1), positions(:,2));
plot(vx, vy, '-', 'Color', [0.4, 0.4, 0.4, 0.6], 'LineWidth', 1);
hold off;

colormap(rocket);
caxis([0, 1]);
cb = colorbar;
cb.Label.String = '$\mathrm{Order\ Parameter}$';  
cb.Label.Interpreter = 'latex';                   
cb.Label.FontWeight = 'normal';
cb.TickLabelInterpreter='latex'
% Professional formatting
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 16;
ax.LineWidth = 1;
ax.TickDir = 'out';
grid on;
ax.GridAlpha = 0.1;
ax.TickLabelInterpreter='latex';
box on

xlabel('$x \ (\mu m)$', 'FontSize', 18,'Interpreter','latex');
ylabel('$y \ (\mu m)$', 'FontSize', 18,'Interpreter','latex');
title('Hexagonal  Order Parameter', 'FontSize', 16,'Interpreter','latex');
% baseDir = pwd;
% fileName = 'Q6_test';
% savefigures(baseDir, fileName);
%%
% figure('Position', [100, 100, 800, 600], 'Color', 'white');
% imagen = imread('exp0000.tif');
% imagen = imadjust(imagen,[0.1 0.7],[]);
% if size(imagen, 3) == 3
%     imagen = rgb2gray(imagen);  % Convert to grayscale if needed
% end
% scale = 5.5667;  % px/um
% pixel_to_um = 1/scale;  % um/px
% [imgHeight, imgWidth] = size(imagen);
% imgWidth_um = imgWidth * pixel_to_um;
% imgHeight_um = imgHeight * pixel_to_um;
% 
% xMin = min(positions(:,1)); xMax = max(positions(:,1));
% yMin = min(positions(:,2)); yMax = max(positions(:,2));
% padding = 0.08 * max(xMax - xMin, yMax - yMin);
% xLimits = [xMin - padding, xMax + padding];
% yLimits = [yMin - padding, yMax + padding];
% 
% % Convert grayscale image to RGB so it won't be affected by colormap changes
% imagen_rgb = repmat(imagen, [1, 1, 3]);  %
% 
% imshow(imagen_rgb, 'XData', [0, imgWidth_um], 'YData', [0, imgHeight_um]);
% set(gca, 'YDir', 'normal');  % Correct Y-axis direction
% 
% hold on;
% 
% scatter(positions(:,1), positions(:,2), 25, Q6, 'filled', ...
%         'MarkerEdgeColor', [0.1, 0.1, 0.1], 'LineWidth', 0.8, ...
%         'MarkerFaceAlpha', 0.95);
% % 
% % scatter(positions(:,1), positions(:,2), 125, Q6, ...
% %         'MarkerEdgeColor', 'flat', 'LineWidth', 1.2, ...
% %         'MarkerFaceColor', 'none');
% [vx, vy] = voronoi(positions(:,1), positions(:,2));
% plot(vx, vy, '-', 'Color', 'w', 'LineWidth', 0.5);
% 
% colormap(turbo);
% caxis([0, 1]);
% 
% hold off;
% axis equal;
% xlim(xLimits);
% ylim(yLimits);
% 
% cb = colorbar;
% cb.Label.String = '$\mathrm{Order\ Parameter}$';
% cb.Label.Interpreter = 'latex';
% cb.Label.FontWeight = 'normal';
% cb.TickLabelInterpreter = 'latex';
% 
% ax = gca;
% ax.FontName = 'Arial';
% ax.FontSize = 16;
% ax.LineWidth = 1;
% ax.TickDir = 'out';
% grid on;
% ax.GridAlpha = 0.1;
% ax.GridColor = [0.8, 0.8, 0.8];  % Light grid to show on dark background
% ax.TickLabelInterpreter = 'latex';
% box on;
% 
% xlabel('$x \ (\mu m)$', 'FontSize', 18, 'Interpreter', 'latex');
% ylabel('$y \ (\mu m)$', 'FontSize', 18, 'Interpreter', 'latex');
% title('Hexagonal Order Parameter', 'FontSize', 16, 'Interpreter', 'latex');
% baseDir = pwd;
% fileName = 'Q6_test2';
% savefigures(baseDir, fileName);
%%
% Add this diagnostic after your Q6 calculation to see the issue:
figure;
scatter(positions(:,1), positions(:,2), 50, 'k', 'filled');
hold on;

% Plot ALL Delaunay connections
edges_list = edges(DT);
for i = 1:size(edges_list, 1)
    p1 = edges_list(i, 1);
    p2 = edges_list(i, 2);
    plot([DT.Points(p1,1), DT.Points(p2,1)], ...
         [DT.Points(p1,2), DT.Points(p2,2)], 'r-', 'LineWidth', 0.5);
end

N = size(positions,1);
k = 6;  % six nearest neighbors

% find the k+1 nearest neighbors (first neighbor is self, so we'll skip it)
[idxs, dists] = knnsearch(positions, positions, 'K', k+1);

psi6 = zeros(N,1);

for i = 1:N
    neigh = idxs(i,2:end);   % skip the first column = self
    dx = positions(neigh,1) - positions(i,1);
    dy = positions(neigh,2) - positions(i,2);
    thetas = atan2(dy, dx);   % bond angles from i to each neighbor
    
    psi6(i) = abs( mean( exp(6i * thetas) ) );
end

% --- now plot colored by |psi6| ---
figure; hold on;
colormap(jet);
scatter(positions(:,1), positions(:,2), 50, psi6, 'filled');
colorbar;
caxis([0 1]);
axis equal; box on;
xlabel('x (\mum)');
ylabel('y (\mum)');
title('Local Hexatic Order |\psi_6|');
%%
