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

%% Find actual time values in the data
all_times = [];
for i = 1:length(T)
    all_times = [all_times; T{i}];
end
unique_times = unique(all_times);
unique_times = sort(unique_times);

fprintf('Found %d unique time points\n', length(unique_times));
fprintf('Time range: %.1f to %.1f\n', min(unique_times), max(unique_times));

% Use actual times from data instead of assuming they start at 1
timeStep = 10; % Process every 10th time point
time_indices = 1:timeStep:length(unique_times);
selected_times = unique_times(time_indices);

fprintf('Will process %d time points\n', length(selected_times));

%% Video parameters
calib = 1.0;           % Calibration factor (μm/pixel)
particle_radius = 0.5; % Particle radius in μm (ADJUST THIS FOR YOUR PARTICLES!)
frameRate = 10;        % frames per second for the output video

% Create video writer object
outputVideo = VideoWriter('local_density_evolution.mp4', 'MPEG-4');
outputVideo.FrameRate = frameRate;
outputVideo.Quality = 100;
open(outputVideo);

% Create figure for video frames
fig = figure('Position', [100, 100, 800, 600], 'Color', 'white');

%% Find global position limits across all time points for consistent scaling
fprintf('Finding global position limits...\n');
allPositions = [];
for t = selected_times'
    positions_temp = [];
    for i = 1:length(X)
        % Use tolerance for floating point comparison
        timeIdx = find(abs(T{i} - t) < 0.01, 1);
        if ~isempty(timeIdx)
            positions_temp = [positions_temp; X{i}(timeIdx), Y{i}(timeIdx)];
        end
    end
    if ~isempty(positions_temp)
        allPositions = [allPositions; positions_temp];
    end
end

% Check if we found any positions
if isempty(allPositions)
    error('No particle positions found! Check that time values match between tracks.');
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

fprintf('Position limits: X=[%.1f, %.1f], Y=[%.1f, %.1f]\n', ...
        xLimits(1), xLimits(2), yLimits(1), yLimits(2));

% Calculate reference area for close-packed 2D hexagonal lattice
a = particle_radius;
A_ocp = 2 * sqrt(3) * a^2;
phi_ocp = pi / (2*sqrt(3));

fprintf('\nParticle radius: %.2f μm\n', a);
fprintf('Close-packed reference area A_ocp: %.3f μm²\n', A_ocp);
fprintf('Close-packing fraction φ_ocp: %.3f\n', phi_ocp);

%% Main loop to create video frames
fprintf('\nCreating video frames...\n');
totalFrames = length(selected_times);
frameCount = 0;

% Storage for time evolution data
time_array = [];
rho_avg_array = [];
rho_std_array = [];
N_particles_array = [];

for targetTime = selected_times'
    frameCount = frameCount + 1;
    
    % Clear the current figure
    clf(fig);
    
    % Extract positions at current time
    positions = [];
    for i = 1:length(X)
        % Use tolerance for floating point comparison
        timeIdx = find(abs(T{i} - targetTime) < 0.01, 1);
        if ~isempty(timeIdx)
            positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
        end
    end
    positions = positions * calib;
    
    % Skip if not enough particles at this time
    if isempty(positions) || size(positions,1) < 4
        fprintf('Skipping time %.1f: only %d particles\n', targetTime, size(positions,1));
        continue;
    end
    
    %% Calculate local density parameter using Voronoi tessellation
    N = size(positions,1);
    
    % Compute Voronoi tessellation
    try
        [V, C] = voronoin(positions);
    catch ME
        fprintf('Voronoi failed at time %.1f: %s\n', targetTime, ME.message);
        continue;
    end
    
    % Initialize the density parameter array
    rho = zeros(N,1);
    valid_particles = false(N,1);
    
    % Calculate a reasonable boundary for filtering
    if N >= 3
        try
            hull_indices = convhull(positions(:,1), positions(:,2));
            hull_points = positions(hull_indices, :);
        catch
            % If convhull fails, use all points
            hull_points = positions;
        end
    else
        hull_points = positions;
    end
    
    center_x = mean(positions(:,1));
    center_y = mean(positions(:,2));
    
    % Maximum reasonable distance for a Voronoi vertex
    max_distance = 2 * max(xMax_global - xMin_global, yMax_global - yMin_global);
    
    % Calculate ρ for each particle
    for i = 1:N
        vertices = C{i};
        
        % Skip cells with infinite vertices (indicated by vertex index 1)
        if any(vertices == 1)
            rho(i) = NaN;
            continue;
        end
        
        % Get the coordinates of this cell's vertices
        cell_x = V(vertices, 1);
        cell_y = V(vertices, 2);
        
        % Check if all vertices are within reasonable bounds
        distances = sqrt((cell_x - center_x).^2 + (cell_y - center_y).^2);
        
        if all(distances < max_distance) && length(vertices) >= 3
            % Calculate Voronoi cell area
            A_j = polyarea(cell_x, cell_y);
            
            % Only accept positive, finite areas
            if A_j > 0 && isfinite(A_j)
                % Local density parameter: ρⱼ = A_ocp / Aⱼ
                rho(i) = A_ocp / A_j;
                valid_particles(i) = true;
            else
                rho(i) = NaN;
            end
        else
            rho(i) = NaN;
        end
    end
    
    % Filter out invalid particles for statistics and plotting
    valid_rho = rho(valid_particles);
    valid_positions = positions(valid_particles, :);
    
    if isempty(valid_rho)
        fprintf('No valid particles at time %.1f\n', targetTime);
        continue;
    end
    
    % Store time evolution data
    time_array = [time_array; targetTime];
    rho_avg_array = [rho_avg_array; mean(valid_rho)];
    rho_std_array = [rho_std_array; std(valid_rho)];
    N_particles_array = [N_particles_array; length(valid_rho)];
    
    %% Create the visualization
    % Scatter plot colored by ρ
    scatter(valid_positions(:,1), valid_positions(:,2), 60, valid_rho, 'filled', ...
        'MarkerEdgeColor', [0.2, 0.2, 0.2], 'LineWidth', 0.8, ...
        'MarkerFaceAlpha', 1);
    
    % Set consistent axes limits
    axis equal;
    xlim(xLimits);
    ylim(yLimits);
    
    % Add bounded Voronoi cells
    hold on;
    
    % Expand the hull slightly for boundary
    if size(hull_points, 1) > 2
        expanded_hull = zeros(size(hull_points));
        expansion_factor = 1.1;
        
        for i = 1:size(hull_points, 1)
            dx = hull_points(i, 1) - center_x;
            dy = hull_points(i, 2) - center_y;
            expanded_hull(i, 1) = center_x + dx * expansion_factor;
            expanded_hull(i, 2) = center_y + dy * expansion_factor;
        end
    else
        expanded_hull = hull_points;
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
        distances = sqrt((cell_x - center_x).^2 + (cell_y - center_y).^2);
        
        if all(distances < max_distance) && length(vertices) >= 3
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
    
    % Add reference lines for phase boundaries
    rho_fluid_max = 0.70 / phi_ocp;    % ≈ 0.77
    rho_hexatic_max = 0.72 / phi_ocp;  % ≈ 0.79
    
    hold off;
    
    % Colormap and color scale
    colormap(hot);
    caxis([0, 1.0]);
    cb = colorbar;
    cb.Label.String = '$\rho = A_{\mathrm{ocp}}/A$';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontWeight = 'normal';
    cb.Label.FontSize = 14;
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
    box on;
    
    % Labels and title
    xlabel('$x \ (\mu\mathrm{m})$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$y \ (\mu\mathrm{m})$', 'FontSize', 18, 'Interpreter', 'latex');
    title_str = sprintf('Local Density $\\rho$ (t = %.1f) | $\\langle\\rho\\rangle = %.2f$', ...
                        targetTime, mean(valid_rho));
    title(title_str, 'FontSize', 16, 'Interpreter', 'latex');
    
    % Add phase info
    phi_avg = mean(valid_rho) * phi_ocp;
    if phi_avg < 0.70
        phase_str = 'FLUID';
        phase_color = [0, 0.7, 1];
    elseif phi_avg >= 0.70 && phi_avg < 0.72
        phase_str = 'HEXATIC';
        phase_color = [1, 0.7, 0];
    else
        phase_str = 'CRYSTAL';
        phase_color = [1, 0.2, 0.2];
    end
    
    text(0.02, 0.98, sprintf('Phase: %s ($\\phi \\approx %.2f$)', phase_str, phi_avg), ...
         'Units', 'normalized', 'FontSize', 14, 'Interpreter', 'latex', ...
         'VerticalAlignment', 'top', 'BackgroundColor', [1 1 1 0.8], ...
         'EdgeColor', phase_color, 'LineWidth', 2);
    
    % Capture the frame
    drawnow;
    frame = getframe(fig);
    writeVideo(outputVideo, frame);
    
    % Display progress
    if mod(frameCount, 10) == 0 || frameCount == totalFrames
        fprintf('Frame %d/%d (t=%.1f, <ρ>=%.2f, N=%d)\n', ...
                frameCount, totalFrames, targetTime, mean(valid_rho), length(valid_rho));
    end
end

% Close the video file
close(outputVideo);
fprintf('\nVideo generation complete!\n');
fprintf('Total frames processed: %d\n', frameCount);
fprintf('Video saved as: local_density_evolution.mp4\n');

%% Create summary plot of time evolution
if ~isempty(time_array)
    figure('Position', [100, 100, 1000, 400], 'Color', 'white');
    
    subplot(1,2,1)
    errorbar(time_array, rho_avg_array, rho_std_array, 'o-', ...
             'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'auto');
    hold on;
    yline(rho_fluid_max, '--', 'Color', 'cyan', 'LineWidth', 2, ...
          'Label', 'Fluid/Hexatic', 'LabelHorizontalAlignment', 'left');
    yline(rho_hexatic_max, '--', 'Color', [1 0.5 0], 'LineWidth', 2, ...
          'Label', 'Hexatic/Crystal', 'LabelHorizontalAlignment', 'left');
    hold off;
    xlabel('Time (frames)', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$\langle\rho\rangle$', 'FontSize', 16, 'Interpreter', 'latex');
    title('Average Local Density vs Time', 'FontSize', 14, 'Interpreter', 'latex');
    grid on;
    set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
    
    subplot(1,2,2)
    plot(time_array, N_particles_array, 'o-', 'LineWidth', 2, ...
         'MarkerSize', 6, 'MarkerFaceColor', 'auto');
    xlabel('Time (frames)', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Number of Particles', 'FontSize', 14, 'Interpreter', 'latex');
    title('Particle Count vs Time', 'FontSize', 14, 'Interpreter', 'latex');
    grid on;
    set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
    
    saveas(gcf, 'density_time_evolution.png');
    fprintf('Summary plot saved as: density_time_evolution.png\n');
    
    %% Save data
    output_data = table(time_array, rho_avg_array, rho_std_array, N_particles_array, ...
                        'VariableNames', {'Time', 'Rho_Mean', 'Rho_Std', 'N_Particles'});
    writetable(output_data, 'density_evolution_data.csv');
    fprintf('Data saved as: density_evolution_data.csv\n');
end

close(fig);