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

%% Video parameters
calib = 1.0;
startTime = 1;
endTime = 1000;
timeStep = 10; % Process every 10th frame
frameRate = 10; % frames per second for the output video (adjusted for time step)

% Create video writer object
outputVideo = VideoWriter('hexatic_order_evolution.mp4', 'MPEG-4');
outputVideo.FrameRate = frameRate;
outputVideo.Quality = 100; % Maximum quality
open(outputVideo);

% Create figure for video frames
fig = figure('Position', [100, 100, 800, 600], 'Color', 'white');

% Find global position limits across all time points for consistent scaling
allPositions = [];
for t = startTime:timeStep:endTime
    positions_temp = [];
    for i = 1:length(X)
        timeIdx = find(T{i} == t, 1);
        if ~isempty(timeIdx)
            positions_temp = [positions_temp; X{i}(timeIdx), Y{i}(timeIdx)];
        end
    end
    if ~isempty(positions_temp)
        allPositions = [allPositions; positions_temp];
    end
end
allPositions = allPositions * calib;

% Calculate global limits for consistent view
xMin_global = min(allPositions(:,1)); 
xMax_global = max(allPositions(:,1));
yMin_global = min(allPositions(:,2)); 
yMax_global = max(allPositions(:,2));
padding = 0.08 * max(xMax_global - xMin_global, yMax_global - yMin_global);
xLimits = [xMin_global - padding, xMax_global + padding];
yLimits = [yMin_global - padding, yMax_global + padding];

%% Main loop to create video frames
fprintf('Creating video frames...\n');
frameCount = 0;
for targetTime = startTime:timeStep:endTime
    frameCount = frameCount + 1;
    % Clear the current figure
    clf(fig);
    
    % Extract positions at current time
    positions = [];
    for i = 1:length(X)
        timeIdx = find(T{i} == targetTime, 1);
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
        end
    end
    positions = positions * calib;
    
    % Skip if no particles at this time
    if isempty(positions)
        continue;
    end
    
    % Compute Delaunay triangulation
    DT = delaunayTriangulation(positions(:,1), positions(:,2));
    nPoints = size(DT.Points, 1);
    
    % Calculate hexatic order parameter Q6
    Q6 = zeros(nPoints, 1);
    for k = 1:nPoints
        orderIndex6 = 0;
        neighCount = 0;
        for j = 1:nPoints
            if k ~= j && isConnected(DT, k, j)
                alpha = atan2(DT.Points(j,2) - DT.Points(k,2), ...
                             DT.Points(j,1) - DT.Points(k,1));
                orderIndex6 = orderIndex6 + exp(6*1i*alpha);
                neighCount = neighCount + 1;
            end
        end
        if neighCount > 0
            Q6(k) = abs(orderIndex6 / neighCount);
        end
    end
    
    % Create the plot
    scatter(positions(:,1), positions(:,2), 60, Q6, 'filled', ...
        'MarkerEdgeColor', [0.2, 0.2, 0.2], 'LineWidth', 0.8, ...
        'MarkerFaceAlpha', 1);
    
    % Set consistent axes limits
    axis equal;
    xlim(xLimits);
    ylim(yLimits);
    
    % Add Voronoi cells
    hold on;
    [vx, vy] = voronoi(positions(:,1), positions(:,2));
    plot(vx, vy, '-', 'Color', [0.4, 0.4, 0.4, 0.6], 'LineWidth', 1);
    hold off;
    
    % Colormap and color scale
    colormap(rocket);
    caxis([0, 1]); % Fixed color scale for consistency
    cb = colorbar;
    cb.Label.String = '$\mathrm{Order\ Parameter}$';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontWeight = 'normal';
    cb.TickLabelInterpreter = 'latex';
    
    % Format the axes
    ax = gca;
    ax.FontName = 'Arial';
    ax.FontSize = 16;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    grid on;
    ax.GridAlpha = 0.1;
    ax.TickLabelInterpreter = 'latex';
    box on
    
    % Labels and title with time information
    xlabel('$x \ (\mu m)$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$y \ (\mu m)$', 'FontSize', 18, 'Interpreter', 'latex');
    title(sprintf('Hexagonal Order Parameter (t = %d)', targetTime), ...
          'FontSize', 16, 'Interpreter', 'latex');
    
    % Capture the frame
    drawnow; % Ensure the plot is fully rendered
    frame = getframe(fig);
    writeVideo(outputVideo, frame);
    
    % Display progress
    if mod(targetTime, 50) == 0
        fprintf('Processed frame %d of %d\n', targetTime, endTime);
    end
end

% Close the video file
close(outputVideo);
fprintf('Video saved as hexatic_order_evolution.mp4\n');

% Close the figure
close(fig);

%% Helper function to check if two points are connected in Delaunay triangulation
function connected = isConnected(DT, i, j)
    % Find all triangles containing point i
    triangles = vertexAttachments(DT, i);
    triangles = triangles{1};
    
    % Check if point j is in any of these triangles
    connected = false;
    for t = 1:length(triangles)
        if any(DT.ConnectivityList(triangles(t), :) == j)
            connected = true;
            break;
        end
    end
end