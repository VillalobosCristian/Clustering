clear all; clc; close all;

%% PARAMETERS
pixels_per_micron = 16;  % Calibration for 40x objective
fps = 10;  % Frame rate
time_per_frame = 1 / fps;

% PARTICLE DETECTION PARAMETERS
particle_radius_um = 1.0;  % Particle radius in microns (adjust for your particles)
particle_radius_px = particle_radius_um * pixels_per_micron;

%% COORDINATE ADJUSTMENT OPTIONS
% If positions are off, try adjusting these:
FLIP_Y = false;  % Set to true if Y-axis is inverted
X_OFFSET = 0;    % Offset in microns to add to X coordinates
Y_OFFSET = 0;    % Offset in microns to add to Y coordinates

%% IMAGE FILES AND CORRESPONDING TIMES (in frames)
image_files = {'t_0.jpg', 't_1266.jpg', 't_2531.jpg', 't_3796.jpg'};
frame_times = [0, 1266, 2531, 3796];  % Frame numbers
real_times = frame_times * time_per_frame;  % Convert to seconds

fprintf('Processing %d time points:\n', length(frame_times));
for i = 1:length(frame_times)
    fprintf('  %s: frame %d (t = %.1f s)\n', image_files{i}, frame_times(i), real_times(i));
end

%% LOAD CSV TRACKING DATA
csv_file = 'ps_particles_milliQ_60x_spots.csv';

if ~exist(csv_file, 'file')
    error('CSV file not found: %s', csv_file);
end

fprintf('\nLoading tracking data from: %s\n', csv_file);
data = readtable(csv_file);
data = data(4:end,:);  % Remove TrackMate headers

% Convert to numeric if needed
if iscell(data.TRACK_ID)
    data.TRACK_ID = cellfun(@str2double, data.TRACK_ID);
    data.POSITION_X = cellfun(@str2double, data.POSITION_X);
    data.POSITION_Y = cellfun(@str2double, data.POSITION_Y);
    data.POSITION_T = cellfun(@str2double, data.POSITION_T);
end

sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

fprintf('Loaded %d tracks\n', length(unique_tracks));

%% DETECT TIME FORMAT (frames vs seconds)
first_track = sorted_data(sorted_data.TRACK_ID == unique_tracks(1), :);
sample_times = first_track.POSITION_T(1:min(10, height(first_track)));
median_dt = median(diff(sample_times));

if median_dt > 0.5  % Likely frame numbers
    time_is_frames = true;
    fprintf('Time detected as FRAME NUMBERS → converting to seconds\n');
else  % Already in seconds
    time_is_frames = false;
    fprintf('Time detected as SECONDS → no conversion needed\n');
end

%% EXTRACT ALL TRAJECTORIES
X = {};
Y = {};
T = {};

for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    
    X{end+1} = track_data.POSITION_X + X_OFFSET;  % Apply offset
    Y{end+1} = track_data.POSITION_Y + Y_OFFSET;  % Apply offset
    
    if time_is_frames
        T{end+1} = track_data.POSITION_T * time_per_frame;
    else
        T{end+1} = track_data.POSITION_T;
    end
end

fprintf('Extracted %d trajectories\n', length(X));

%% DIAGNOSTIC: Print coordinate ranges
all_x = cell2mat(X');
all_y = cell2mat(Y');
fprintf('\n=== COORDINATE DIAGNOSTICS ===\n');
fprintf('TrackMate X range: %.1f to %.1f μm\n', min(all_x), max(all_x));
fprintf('TrackMate Y range: %.1f to %.1f μm\n', min(all_y), max(all_y));
fprintf('Image dimensions at %.1f pixels/μm:\n', pixels_per_micron);
fprintf('  X: %.1f to %.1f pixels\n', min(all_x)*pixels_per_micron, max(all_x)*pixels_per_micron);
fprintf('  Y: %.1f to %.1f pixels\n\n', min(all_y)*pixels_per_micron, max(all_y)*pixels_per_micron);

%% CREATE OUTPUT DIRECTORY
output_dir = 'trajectory_overlay_results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% SET UP PLOTTING DEFAULTS
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

%% PROCESS FIRST IMAGE FOR COORDINATE SYSTEM CHECK
fprintf('=== CHECKING FIRST IMAGE FOR COORDINATE ALIGNMENT ===\n');
first_img_file = image_files{1};
first_img = imread(first_img_file);
if size(first_img, 3) == 3
    first_img = rgb2gray(first_img);
end
[img_height, img_width] = size(first_img);
fprintf('Image size: %d x %d pixels\n', img_width, img_height);
fprintf('Expected size from data: %.1f x %.1f pixels\n\n', ...
    max(all_x)*pixels_per_micron, max(all_y)*pixels_per_micron);

%% PROCESS EACH TIME POINT
for t_idx = 1:length(frame_times)
    target_time = real_times(t_idx);
    img_file = image_files{t_idx};
    
    fprintf('Processing %s (t = %.1f s)...\n', img_file, target_time);
    
    % Load image
    if ~exist(img_file, 'file')
        fprintf('  WARNING: Image file not found: %s\n', img_file);
        continue;
    end
    
    img = imread(img_file);
    [img_height, img_width, ~] = size(img);
    
    % Convert to grayscale if color
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    
    % Extract positions at target time
    positions_um = [];
    for i = 1:length(X)
        % Find closest time point in trajectory
        [time_diff, timeIdx] = min(abs(T{i} - target_time));
        
        % Only include if within 0.5 seconds of target time
        if time_diff < 0.5 && ~isempty(timeIdx)
            x_pos = X{i}(timeIdx);
            y_pos = Y{i}(timeIdx);
            
            % Apply Y-flip if needed
            if FLIP_Y
                y_pos = max(all_y) - y_pos + min(all_y);
            end
            
            positions_um = [positions_um; x_pos, y_pos];
        end
    end
    
    % Convert positions from microns to pixels
    positions_px = positions_um * pixels_per_micron;
    
    fprintf('  Found %d particles at t = %.1f s\n', size(positions_um, 1), target_time);
    fprintf('  (This is from CSV - used only for reference/diagnostics)\n');
    
    % Create figure
    fig = figure('Position', [100, 100, 1000, 900], 'Color', 'white');
    
    % Show image as background
    imshow(img, []);
    hold on;
    
    % Plot trajectories UP TO this time point (in pixels)
    for i = 1:length(X)
        % Get all points up to target time
        valid_idx = T{i} <= target_time;
        
        if sum(valid_idx) >= 2  % Need at least 2 points for a trajectory
            x_traj_um = X{i}(valid_idx);
            y_traj_um = Y{i}(valid_idx);
            
            % Apply Y-flip if needed
            if FLIP_Y
                y_traj_um = max(all_y) - y_traj_um + min(all_y);
            end
            
            % Convert to pixels
            x_traj_px = x_traj_um * pixels_per_micron;
            y_traj_px = y_traj_um * pixels_per_micron;
            
            % Plot trajectory (thin, semi-transparent line)
            plot(x_traj_px, y_traj_px, '-', ...
                 'Color', [0.3, 0.6, 1, 0.4], 'LineWidth', 1);
        end
    end
    
    % Overlay current particle positions
    if t_idx == 1
        % First image: use CSV positions (ground truth)
        fprintf('  Using CSV positions for first image\n');
        if ~isempty(positions_px)
            scatter(positions_px(:,1), positions_px(:,2), 25, [1, 0.2, 0.2], 'filled', ...
                'MarkerEdgeColor', [1, 1, 1], 'LineWidth', 1.5, ...
                'MarkerFaceAlpha', 0.8);
        end
        n_particles = size(positions_px, 1);
    else
        % Other images: detect particles from image
        fprintf('  Detecting particles from image...\n');
        radius_range = round([particle_radius_px*0.7, particle_radius_px*1.3]);
        sensitivity = 0.92;
        edge_threshold = 0.1;
        
        [centers_px, radii_px] = imfindcircles(img, radius_range, ...
            'ObjectPolarity', 'bright', ...
            'Sensitivity', sensitivity, ...
            'EdgeThreshold', edge_threshold);
        
        fprintf('  Detected %d particles\n', length(radii_px));
        
        if isempty(centers_px)
            warning('No particles detected in %s!', img_file);
            n_particles = 0;
        else
            % Plot detected particles
            scatter(centers_px(:,1), centers_px(:,2), 25, [1, 0.2, 0.2], "filled", ...
                'MarkerEdgeColor', [1, 1, 1], 'LineWidth', 1);
            n_particles = size(centers_px, 1);
        end
    end
    
    % Add time label
    text(50, 50, sprintf('$t = %.1f$ s', target_time), ...
        'Color', 'white', 'FontSize', 24, 'Interpreter', 'latex', ...
        'BackgroundColor', [0, 0, 0, 0.5], 'EdgeColor', 'white', ...
        'FontWeight', 'bold', 'Margin', 5);
    
    % Add scale bar
    scale_bar_length_um = 20;  % 20 μm
    scale_bar_length_px = scale_bar_length_um * pixels_per_micron;
    x_bar = [50, 50 + scale_bar_length_px];
    y_bar = [img_height - 30, img_height - 30];
    plot(x_bar, y_bar, 'w-', 'LineWidth', 4);
    text(50 + scale_bar_length_px/2, img_height - 50, ...
        sprintf('%d $\\mu$m', scale_bar_length_um), ...
        'Color', 'white', 'FontSize', 24, 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    
    % % Add particle count
    % text(img_width - 50, 50, sprintf('$N = %d$', n_particles), ...
    %     'Color', 'white', 'FontSize', 16, 'Interpreter', 'latex', ...
    %     'HorizontalAlignment', 'right', ...
    %     'BackgroundColor', [0, 0, 0, 0.5], 'EdgeColor', 'white', ...
    %     'Margin', 4);
    % 
    % Format axes
    ax = gca;
    ax.FontSize = 14;
    ax.TickLabelInterpreter = 'latex';
    ax.XColor = 'white';
    ax.YColor = 'white';
    
    hold off;
    
    % Save figure
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [1, 1, 9, 9]);
    
    output_filename = sprintf('overlay_frame_%04d_t_%.1fs', frame_times(t_idx), target_time);
    exportgraphics(gcf, fullfile(output_dir, [output_filename, '.png']), 'Resolution', 300);
    export_fig(output_filename, '-pdf', '-painters', '-r300', '-transparent');

    fprintf('  Saved: %s.png\n\n', output_filename);
    
    % Don't close first figure - keep it open for inspection
    if t_idx > 1
        close(fig);
    end
end

% %% CREATE A 4-PANEL COMPARISON FIGURE
% fprintf('Creating 4-panel comparison figure...\n');
% 
% fig_compare = figure('Position', [100, 100, 150, 150], 'Color', 'white');
% 
% for t_idx = 1:length(frame_times)
%     target_time = real_times(t_idx);
%     img_file = image_files{t_idx};
% 
%     if ~exist(img_file, 'file')
%         continue;
%     end
% 
%     % Load image
%     img = imread(img_file);
%     if size(img, 3) == 3
%         img = rgb2gray(img);
%     end
%     [img_height, img_width] = size(img);
% 
%     % Extract positions
%     positions_um = [];
%     for i = 1:length(X)
%         [time_diff, timeIdx] = min(abs(T{i} - target_time));
%         if time_diff < 0.5 && ~isempty(timeIdx)
%             x_pos = X{i}(timeIdx);
%             y_pos = Y{i}(timeIdx);
% 
%             if FLIP_Y
%                 y_pos = max(all_y) - y_pos + min(all_y);
%             end
% 
%             positions_um = [positions_um; x_pos, y_pos];
%         end
%     end
%     positions_px = positions_um * pixels_per_micron;
% 
%     % Create subplot
%     subplot(2, 2, t_idx);
%     imshow(img, []);
%     hold on;
% 
%     % Plot trajectories
%     for i = 1:length(X)
%         valid_idx = T{i} <= target_time;
%         if sum(valid_idx) >= 2
%             x_traj_um = X{i}(valid_idx);
%             y_traj_um = Y{i}(valid_idx);
% 
%             if FLIP_Y
%                 y_traj_um = max(all_y) - y_traj_um + min(all_y);
%             end
% 
%             x_traj_px = x_traj_um * pixels_per_micron;
%             y_traj_px = y_traj_um * pixels_per_micron;
% 
%             plot(x_traj_px, y_traj_px, '-', ...
%                  'Color', [0.3, 0.6, 1, 0.3], 'LineWidth', 0.8);
%         end
%     end
% 
%     % Plot positions or detect particles
%     if t_idx == 1
%         % First panel: use CSV positions
%         if ~isempty(positions_px)
%             scatter(positions_px(:,1), positions_px(:,2), 5, [1, 0.2, 0.2], 'filled', ...
%                 'MarkerEdgeColor', [1, 1, 1], 'LineWidth', 1.2, ...
%                 'MarkerFaceAlpha', 0.8);
%         end
%         n_particles_panel = size(positions_px, 1);
%     else
%         % Other panels: detect particles
%         radius_range = round([particle_radius_px*0.7, particle_radius_px*1.3]);
%         sensitivity = 0.92;
%         edge_threshold = 0.1;
% 
%         [centers_px, radii_px] = imfindcircles(img, radius_range, ...
%             'ObjectPolarity', 'bright', ...
%             'Sensitivity', sensitivity, ...
%             'EdgeThreshold', edge_threshold);
% 
%         if ~isempty(centers_px)
%             scatter(centers_px(:,1), centers_px(:,2), 5, [1, 0.2, 0.2], 'filled', ...
%                 'MarkerEdgeColor', [1, 1, 1], 'LineWidth', 1.2, ...
%                 'MarkerFaceAlpha', 0.8);
%             n_particles_panel = size(centers_px, 1);
%         else
%             n_particles_panel = 0;
%         end
%     end
% 
%     % Time label
%     text(50, 50, sprintf('$t = %.1f$ s', target_time), ...
%         'Color', 'white', 'FontSize', 18, 'Interpreter', 'latex', ...
%         'BackgroundColor', [0, 0, 0, 0.5], 'EdgeColor', 'white', ...
%         'FontWeight', 'bold', 'Margin', 4);
% 
%     % Particle count
%     text(img_width - 50, 50, sprintf('$N = %d$', n_particles_panel), ...
%         'Color', 'white', 'FontSize', 14, 'Interpreter', 'latex', ...
%         'HorizontalAlignment', 'right', ...
%         'BackgroundColor', [0, 0, 0, 0.5], 'EdgeColor', 'white', ...
%         'Margin', 3);
% 
%     % Scale bar (only on first panel)
%     if t_idx == 1
%         scale_bar_length_um = 20;
%         scale_bar_length_px = scale_bar_length_um * pixels_per_micron;
%         x_bar = [50, 50 + scale_bar_length_px];
%         y_bar = [img_height - 30, img_height - 30];
%         plot(x_bar, y_bar, 'w-', 'LineWidth', 3);
%         text(50 + scale_bar_length_px/2, img_height - 50, ...
%             sprintf('%d $\\mu$m', scale_bar_length_um), ...
%             'Color', 'white', 'FontSize', 12, 'Interpreter', 'latex', ...
%             'HorizontalAlignment', 'center', 'FontWeight', 'bold');
%     end
% 
%     ax = gca;
%     ax.FontSize = 12;
%     ax.TickLabelInterpreter = 'latex';
%     ax.XColor = 'white';
%     ax.YColor = 'white';
% 
%     hold off;
% end
% 
% % Save 4-panel figure
% set(gcf, 'Units', 'inches');
% set(gcf, 'Position', [1, 1, 18, 18]);
% exportgraphics(gcf, fullfile(output_dir, 'comparison_4panel.png'), 'Resolution', 300);
% 
% fprintf('Saved: comparison_4panel.png\n');
% 
% fprintf('\n=== COMPLETE ===\n');
% fprintf('All figures saved to: %s\n', output_dir);
% fprintf('\nApproach used:\n');
% fprintf('  - Image 1 (t=0): Particles from CSV (ground truth)\n');
% fprintf('  - Images 2-4: Particles detected from images using imfindcircles\n');
% fprintf('  - All images: Trajectories from CSV tracking data\n');
% fprintf('\nIf positions are misaligned:\n');
% fprintf('  1. Check the first figure window\n');
% fprintf('  2. Adjust FLIP_Y, X_OFFSET, Y_OFFSET at top of script\n');
% fprintf('  3. Re-run the script\n');
% fprintf('\nIf particle detection needs adjustment:\n');
% fprintf('  1. Adjust particle_radius_um at top of script\n');
% fprintf('  2. Modify sensitivity (0.85-0.95) or edge_threshold (0.05-0.2)\n');
% fprintf('  3. Re-run the script\n');