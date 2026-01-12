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

% Sort data for easier processing
sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

%% Load original images
fprintf('Selecting image directory...\n');
img_path = uigetdir(path_name, 'Select directory containing exp10000.tif image sequence');
if isequal(img_path, 0), return; end

% Find the image sequence starting from exp10000.tif
start_frame = 10000;
img_files = {};
frame_numbers = [];

% Check for .tif and .tiff extensions
for i = 0:10000  % Check up to 10000 frames
    current_frame = start_frame + i;
    img_name_tif = sprintf('exp%05d.tif', current_frame);
    img_name_tiff = sprintf('exp%05d.tiff', current_frame);
    
    if exist(fullfile(img_path, img_name_tif), 'file')
        img_files{end+1} = img_name_tif;
        frame_numbers(end+1) = current_frame;
    elseif exist(fullfile(img_path, img_name_tiff), 'file')
        img_files{end+1} = img_name_tiff;
        frame_numbers(end+1) = current_frame;
    else
        % Stop when we can't find the next frame
        break;
    end
end

if isempty(img_files)
    error('No image sequence found starting with exp10000.tif in the selected directory');
end

fprintf('Found %d images in sequence: exp%05d to exp%05d\n', length(img_files), frame_numbers(1), frame_numbers(end));

% Load first image to get dimensions
first_img = imread(fullfile(img_path, img_files{1}));
[img_height, img_width] = size(first_img);
fprintf('Image dimensions: %d x %d\n', img_width, img_height);

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

%% Set up time points for animation
all_times = [];
for i = 1:length(T)
    all_times = [all_times; T{i}];
end
unique_times = unique(all_times);

time_step = 1; 
analysis_times = unique_times(1:time_step:end);

%% Parameters
calib = 15.5; % 16.7 pixels per micrometer - convert from μm to pixels
distance_threshold = 3; % Keep threshold in μm for analysis

all_x = []; all_y = [];
for i = 1:length(X)
    all_x = [all_x; X{i}];
    all_y = [all_y; Y{i}];
end

fprintf('Tracking data ranges in μm - X: %.1f to %.1f, Y: %.1f to %.1f\n', ...
    min(all_x), max(all_x), min(all_y), max(all_y));

% Convert to pixels for display
all_x_pixels = all_x * calib;
all_y_pixels = all_y * calib;

fprintf('Converted to pixels - X: %.1f to %.1f, Y: %.1f to %.1f\n', ...
    min(all_x_pixels), max(all_x_pixels), min(all_y_pixels), max(all_y_pixels));
fprintf('Image dimensions: %d x %d pixels\n', img_width, img_height);

%% Set up video writer
video_filename = [file_name(1:end-4), '_triangles_overlay.mp4'];
video_path = fullfile(path_name, video_filename);
v = VideoWriter(video_path, 'MPEG-4');
v.FrameRate = 8;
v.Quality = 95;
open(v);

%%
figure('Position', [100, 100, 800, 800]);

fprintf('Creating triangle overlay video: %s\n', video_filename);
fprintf('Processing %d time points...\n', length(analysis_times));

% Animation loop
for t_idx = 1:length(analysis_times)
    current_time = analysis_times(t_idx);
    
    % Find corresponding image frame
    target_frame = start_frame + current_time;
    frame_idx = find(frame_numbers == target_frame, 1);
    
    if ~isempty(frame_idx)
        current_img = imread(fullfile(img_path, img_files{frame_idx}));
    else
        % If exact frame not found, use the closest available frame
        [~, closest_idx] = min(abs(frame_numbers - target_frame));
        current_img = imread(fullfile(img_path, img_files{closest_idx}));
        if t_idx == 1
            fprintf('Note: Using closest available frames when exact match not found\n');
        end
    end
    
    % Extract positions at this time point (in micrometers)
    positions_um = [];
    particle_ids = [];
    
    for i = 1:length(X)
        timeIdx = find(T{i} == current_time, 1);
        if ~isempty(timeIdx)
            positions_um = [positions_um; X{i}(timeIdx), Y{i}(timeIdx)];
            particle_ids = [particle_ids; i];
        end
    end
    
    N_particles = size(positions_um, 1);
    
    % Skip if too few particles
    if N_particles < 3
        continue;
    end
    
    % Create Delaunay triangulation in micrometers (for correct distance calculations)
    DT = delaunayTriangulation(positions_um(:,1), positions_um(:,2));
    triangles = DT.ConnectivityList;
    n_triangles = size(triangles, 1);
    
    valid_triangles = [];
    for i = 1:n_triangles
        p1 = triangles(i, 1);
        p2 = triangles(i, 2);
        p3 = triangles(i, 3);
        
        % Calculate the three side lengths IN MICROMETERS
        side1 = norm(positions_um(p1,:) - positions_um(p2,:));
        side2 = norm(positions_um(p2,:) - positions_um(p3,:));
        side3 = norm(positions_um(p3,:) - positions_um(p1,:));
        
        max_side = max([side1, side2, side3]);
        
        % Keep triangle only if sides are below threshold (in μm)
        if max_side <= distance_threshold
            valid_triangles = [valid_triangles; triangles(i,:)];
        end
    end
    
    % Convert positions to pixels for display
    positions_pixels = positions_um * calib;
    
    n_valid_triangles = size(valid_triangles, 1);
    
    % Clear only the axes content, not the entire figure
    cla;
    
    % Display the microscopy image as background
    imshow(current_img, []);
    hold on;
    
    % Draw particles using pixel coordinates
    scatter(positions_pixels(:,1), positions_pixels(:,2), 120, 'w', 'filled', 'MarkerFaceAlpha', 0.8);
    
    % Draw only the valid triangles using pixel coordinates
    for i = 1:n_valid_triangles
        p1 = valid_triangles(i, 1);
        p2 = valid_triangles(i, 2);
        p3 = valid_triangles(i, 3);
        
        % Draw the triangle edges using pixel coordinates
        plot([positions_pixels(p1,1), positions_pixels(p2,1)], [positions_pixels(p1,2), positions_pixels(p2,2)], 'r-', 'LineWidth', 2);
        plot([positions_pixels(p2,1), positions_pixels(p3,1)], [positions_pixels(p2,2), positions_pixels(p3,2)], 'r-', 'LineWidth', 2);
        plot([positions_pixels(p3,1), positions_pixels(p1,1)], [positions_pixels(p3,2), positions_pixels(p1,2)], 'r-', 'LineWidth', 2);
    end
    
    % Set axis limits to image dimensions and hold them steady
    xlim([1, img_width]);
    ylim([1, img_height]);
    axis off; % Remove axis ticks for cleaner look
    
    title(sprintf('Time = %d | Particles = %d | Valid Triangles = %d', ...
        current_time, N_particles, n_valid_triangles), 'Color', 'white', 'FontSize', 12);
    
    % Add progress indicator with coordinate info
    if mod(t_idx, 10) == 0 || t_idx == length(analysis_times)
        fprintf('Progress: %d/%d frames (%.1f%%) - Particles: %d, Triangles: %d\n', ...
            t_idx, length(analysis_times), 100*t_idx/length(analysis_times), N_particles, n_valid_triangles);
        
        if t_idx == 1 && N_particles > 0
            fprintf('First particle: μm(%.1f,%.1f) -> pixels(%.1f,%.1f)\n', ...
                positions_um(1,1), positions_um(1,2), positions_pixels(1,1), positions_pixels(1,2));
        end
    end
    
    % Reduced screen updates - only update every few frames for smoother video
    if mod(t_idx, 2) == 0  % Update screen every 2nd frame
        drawnow;
    end
    
    % Capture frame for video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close video file
close(v);
close(gcf);

fprintf('\nTriangle overlay video saved successfully: %s\n', video_path);
fprintf('Video duration: %.2f seconds at %d fps\n', length(analysis_times)/v.FrameRate, v.FrameRate);