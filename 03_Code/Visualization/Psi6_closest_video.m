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
frameRate = 10; % frames per second for the output video

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
totalFrames = length(startTime:timeStep:endTime);
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
    if isempty(positions) || size(positions,1) < 7
        continue; % Need at least 7 particles for k=6 neighbors
    end
    
    %% Calculate hexatic order parameter using k-nearest neighbors
    N = size(positions,1);
    k = 6; % six nearest neighbors for hexatic order
    
    % Find the k+1 nearest neighbors (first neighbor is self, so we'll skip it)
    % knnsearch returns indices and distances to k+1 nearest points
    [idxs, dists] = knnsearch(positions, positions, 'K', min(k+1, N));
    
    % Initialize the order parameter array
    psi6 = zeros(N,1);
    
    % Calculate psi6 for each particle
    for i = 1:N
        % Get indices of the k nearest neighbors (excluding self)
        if size(idxs,2) > 1
            neigh = idxs(i, 2:min(end, k+1)); % skip the first column (self)
            
            % Calculate bond vectors from particle i to each neighbor
            dx = positions(neigh,1) - positions(i,1);
            dy = positions(neigh,2) - positions(i,2);
            
            % Calculate bond angles using atan2 (gives angles in [-π, π])
            thetas = atan2(dy, dx);
            
            % Calculate the complex hexatic order parameter
            % psi6 = |⟨exp(6iθ)⟩| where ⟨⟩ denotes average over bonds
            psi6(i) = abs(mean(exp(6i * thetas)));
        else
            psi6(i) = 0; % If not enough neighbors, set to 0
        end
    end
    
    %% Create the visualization
    % Scatter plot colored by |psi6|
    scatter(positions(:,1), positions(:,2), 60, psi6, 'filled', ...
        'MarkerEdgeColor', [0.2, 0.2, 0.2], 'LineWidth', 0.8, ...
        'MarkerFaceAlpha', 1);
    
    % Set consistent axes limits
    axis equal;
    xlim(xLimits);
    ylim(yLimits);
    
    % Add bounded Voronoi cells
    hold on;
    
    % Use voronoin for Voronoi diagram
    [V, C] = voronoin(positions);
    
    % Calculate a reasonable boundary for clipping
    hull_indices = convhull(positions(:,1), positions(:,2));
    hull_points = positions(hull_indices, :);
    
    % Expand the hull slightly for boundary
    center_x = mean(positions(:,1));
    center_y = mean(positions(:,2));
    expanded_hull = zeros(size(hull_points));
    expansion_factor = 1.1; % 10% expansion
    
    for i = 1:size(hull_points, 1)
        dx = hull_points(i, 1) - center_x;
        dy = hull_points(i, 2) - center_y;
        expanded_hull(i, 1) = center_x + dx * expansion_factor;
        expanded_hull(i, 2) = center_y + dy * expansion_factor;
    end
    
    % Plot each Voronoi cell with clipping
    for k = 1:length(C)
        vertices = C{k};
        
        % Skip cells with infinite vertices
        if any(vertices == 1)
            continue;
        end
        
        % Get the coordinates of this cell's vertices
        cell_x = V(vertices, 1);
        cell_y = V(vertices, 2);
        
        % Check if all vertices are within reasonable bounds
        max_distance = 2 * max(xMax_global - xMin_global, yMax_global - yMin_global);
        distances = sqrt((cell_x - center_x).^2 + (cell_y - center_y).^2);
        
        if all(distances < max_distance)
            % Check if the cell is inside or near the expanded hull
            in_region = false;
            for v = 1:length(vertices)
                if inpolygon(cell_x(v), cell_y(v), expanded_hull(:,1), expanded_hull(:,2))
                    in_region = true;
                    break;
                end
            end
            
            % Only plot cells that are in the region of interest
            if in_region
                plot([cell_x; cell_x(1)], [cell_y; cell_y(1)], '-', ...
                     'Color', [0.4, 0.4, 0.4, 0.6], 'LineWidth', 1);
            end
        end
    end
    
    % Optional: Draw k-nearest neighbor bonds to visualize the calculation
    % Uncomment to see which bonds are being used for psi6 calculation
    for i = 1:N
        neigh = idxs(i, 2:min(end, k+1));
        for j = 1:length(neigh)
            plot([positions(i,1), positions(neigh(j),1)], ...
                 [positions(i,2), positions(neigh(j),2)], ...
                 'b-', 'LineWidth', 0.5, 'Color', [0 0 1 0.2]);
        end
    end
    
    hold off;
    
    % Colormap and color scale
    colormap(rocket);
    caxis([0, 1]); % Fixed color scale for consistency
    cb = colorbar;
    cb.Label.String = '$|\psi_6|$';
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
    title(sprintf('Local Hexatic Order $|\\psi_6|$ (t = %d)', targetTime), ...
          'FontSize', 16, 'Interpreter', 'latex');
    
    % Capture the frame
    drawnow; % Ensure the plot is fully rendered
    frame = getframe(fig);
    writeVideo(outputVideo, frame);
    
    % Display progress
    if mod(frameCount, 10) == 0 || frameCount == totalFrames
        fprintf('Processed frame %d of %d (time = %d)\n', ...
                frameCount, totalFrames, targetTime);
    end
end

% Close the video file
close(outputVideo);
fprintf('\nVideo generation complete!\n');
fprintf('Total frames processed: %d\n', frameCount);
fprintf('Video saved as: hexatic_order_evolution.mp4\n');
fprintf('Video duration: %.1f seconds at %d fps\n', frameCount/frameRate, frameRate);

% Close the figure
close(fig);