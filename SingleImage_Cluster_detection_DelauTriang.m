clear all; clc; close all;

%%  Load CSV data 
[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), return; end

data = readtable(fullfile(path_name, file_name));
data = data(4:end,:); % Remove TrackMate headers

% Convert to numeric if needed (exactly as in your original code)
if iscell(data.TRACK_ID), data.TRACK_ID = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X), data.POSITION_X = cellfun(@str2double, data.POSITION_X); end
if iscell(data.POSITION_Y), data.POSITION_Y = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T), data.POSITION_T = cellfun(@str2double, data.POSITION_T); end

% Sort data for easier processing
sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

%% : Extract trajectories
min_track_length = 100;  % Keep tracks longer than 100 time points
X = {}; Y = {}; T = {};
count = 0;

fprintf('Processing %d tracks...\n', length(unique_tracks));

for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    if height(track_data) >= min_track_length
        count = count + 1;
        X{count} = track_data.POSITION_X;
        Y{count} = track_data.POSITION_Y;
        T{count} = track_data.POSITION_T;
    end
end

fprintf('Found %d tracks with length >= %d\n', count, min_track_length);

%%
% For now, let's work with a single time point to develop the method
targetTime = 1000;  % Starting 
calib = 1.0;        % 
positions = [];
particle_ids = [];  % Keep track of which particle is which

for i = 1:length(X)
    timeIdx = find(T{i} == targetTime, 1);
    if ~isempty(timeIdx)
        positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
        particle_ids = [particle_ids; i];  % Store which track this position belongs to
    end
end

positions = positions * calib;
N_particles = size(positions, 1);


%%  
figure('Position', [100, 100, 800, 600]);
scatter(positions(:,1), positions(:,2), 50, 'b', 'filled', 'MarkerFaceAlpha', 0.7);
axis equal;
grid on;
xlabel('x (\mum)');
ylabel('y (\mum)');
xRange = max(positions(:,1)) - min(positions(:,1));
yRange = max(positions(:,2)) - min(positions(:,2));
area = xRange * yRange;
density = N_particles / area;
 
fprintf('  Density: %.4f particles/μm²\n', density);
%%
DT = delaunayTriangulation(positions(:,1), positions(:,2));

% Get all triangles from the triangulation
triangles = DT.ConnectivityList;  % Each row is [p1, p2, p3] indices
n_triangles = size(triangles, 1);


%% distance filtering to triangles
% Following the paper's approach keep only triangles where ALL sides are shorter than threshold
distance_threshold = 4.5;  % 

% For each triangle, calculate all three side lengths
valid_triangles = [];
triangle_max_distances = [];

for i = 1:n_triangles
    p1 = triangles(i, 1);
    p2 = triangles(i, 2);
    p3 = triangles(i, 3);
    
    % Calculate the three side lengths
    side1 = norm(positions(p1,:) - positions(p2,:));  % Distance p1-p2
    side2 = norm(positions(p2,:) - positions(p3,:));  % Distance p2-p3
    side3 = norm(positions(p3,:) - positions(p1,:));  % Distance p3-p1
    
    max_side = max([side1, side2, side3]);
    triangle_max_distances = [triangle_max_distances; max_side];
    
    % Keep triangle only if ALL sides are below threshold
    if max_side <= distance_threshold
        valid_triangles = [valid_triangles; triangles(i,:)];
    end
end

n_valid_triangles = size(valid_triangles, 1);
fprintf('  Valid triangles (all sides <= %.1f μm): %d (%.1f%%)\n', ...
    distance_threshold, n_valid_triangles, 100*n_valid_triangles/n_triangles);

%% 
figure('Position', [100, 100, 1200, 500]);
subplot(1, 2, 1);
triplot(DT, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5);
hold on;
scatter(positions(:,1), positions(:,2), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.8);
axis equal; grid on;
title(sprintf('Full Delaunay Triangulation\n(%d triangles)', n_triangles));
xlabel('x (μm)'); ylabel('y (μm)');
subplot(1, 2, 2);
scatter(positions(:,1), positions(:,2), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.8);
hold on;

% Draw only the valid triangles
for i = 1:n_valid_triangles
    p1 = valid_triangles(i, 1);
    p2 = valid_triangles(i, 2);
    p3 = valid_triangles(i, 3);
    
    % Draw the triangle edges
    plot([positions(p1,1), positions(p2,1)], [positions(p1,2), positions(p2,2)], 'r-', 'LineWidth', 1);
    plot([positions(p2,1), positions(p3,1)], [positions(p2,2), positions(p3,2)], 'r-', 'LineWidth', 1);
    plot([positions(p3,1), positions(p1,1)], [positions(p3,2), positions(p1,2)], 'r-', 'LineWidth', 1);
end

axis equal; grid on;
title(sprintf('Distance-Filtered Triangles\n(%d triangles, threshold = %.1f μm)', ...
    n_valid_triangles, distance_threshold));
xlabel('x (μm)'); ylabel('y (μm)');

