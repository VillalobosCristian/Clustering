%% SIMPLE ANIMATED TRIANGULATION (Fixed for duplicates)
fprintf('Starting triangulation animation...\n');

% Create figure
figure('Position', [100 100 800 600]);

% Loop through every 5th time point
for t = 1:50:length(unique_times)
    current_time = unique_times(t);
    
    % Get all particle positions at this time
    pos_x = [];
    pos_y = [];
    

    for i = 1:num_tracks
        [~, time_idx] = min(abs(T{i} - current_time));
        pos_x = [pos_x; X{i}(time_idx)];
        pos_y = [pos_y; Y{i}(time_idx)];
    end
    
    % Remove duplicate positions
    positions = [pos_x, pos_y];
    [unique_positions, ~] = unique(positions, 'rows');
    pos_x = unique_positions(:, 1);
    pos_y = unique_positions(:, 2);
    
    % Make triangulation
    DT = delaunayTriangulation(pos_x, pos_y);
    triangles = DT.ConnectivityList;
    
    % Clear and redraw
    clf;
    hold on;
    
    % Draw triangles
    for tri = 1:size(triangles, 1)
        triangle_x = pos_x(triangles(tri, [1,2,3,1]));
        triangle_y = pos_y(triangles(tri, [1,2,3,1]));
        plot(triangle_x, triangle_y, 'b-', 'LineWidth', 0.8);
    end
    
    % Draw particles
    scatter(pos_x, pos_y, 50, 'red', 'filled');
    
    title(sprintf('Time = %.0f (%d particles, %d triangles)', ...
          current_time, length(pos_x), size(triangles,1)));
    axis equal;
    grid on;
    
    pause(0.2);
    
    fprintf('t = %.0f (%d unique particles)\n', current_time, length(pos_x));
end

fprintf('Animation complete!\n');
%%
%% SIMPLE TRAJECTORY ANIMATION
fprintf('Creating trajectory animation...\n');

% Find time range
all_times = [];
for i = 1:num_tracks
    all_times = [all_times; T{i}];
end
time_min = min(all_times);
time_max = max(all_times);

fprintf('Time range: %.1f to %.1f\n', time_min, time_max);

% Create time vector for animation
num_frames = 10; % Number of animation frames
time_vector = linspace(time_min, time_max, num_frames);

% Create figure
figure('Position', [100 100 800 600]);

% Animation loop
for frame = 1:num_frames
    current_time = time_vector(frame);
    
    clf; % Clear figure
    hold on;
    
    % For each trajectory, plot up to current time
    for i = 1:min(50, num_tracks) % Show first 50 trajectories to avoid clutter
        
        % Find points up to current time
        time_mask = T{i} <= current_time;
        
        if sum(time_mask) > 0
            % Plot trajectory path up to current time
            plot(X{i}(time_mask), Y{i}(time_mask), 'b-', 'LineWidth', 1, 'Color', [0.5 0.5 0.8]);
            
            % Plot current particle position
            if sum(time_mask) > 0
                current_x = X{i}(time_mask);
                current_y = Y{i}(time_mask);
                scatter(current_x(end), current_y(end), 40, 'red', 'filled', 'MarkerEdgeColor', 'black');
            end
        end
    end
    
    title(sprintf('Particle Trajectories - Time = %.1f', current_time));
    xlabel('X position (pixels)');
    ylabel('Y position (pixels)');
    axis equal;
    grid on;
    
    % Keep axis limits fixed
    if frame == 1
        axis_limits = axis;
    else
        axis(axis_limits);
    end
    
    pause(0.1); % Pause between frames
    
    if mod(frame, 10) == 0
        fprintf('Frame %d/%d (t = %.1f)\n', frame, num_frames, current_time);
    end
end

fprintf('Animation complete!\n');
%%
%% TRAJECTORIES + TRIANGULATION ANIMATION
fprintf('Creating trajectory + triangulation animation...\n');

% Find time range
all_times = [];
for i = 1:num_tracks
    all_times = [all_times; T{i}];
end
time_min = min(all_times);
time_max = max(all_times);

fprintf('Time range: %.1f to %.1f\n', time_min, time_max);

% Create time vector for animation
num_frames = 100;
time_vector = linspace(time_min, time_max, num_frames);

% Decide which trajectories to show
num_trajectories_to_show = min(50, num_tracks);
fprintf('Showing %d trajectories\n', num_trajectories_to_show);

% Create figure
figure('Position', [100 100 800 600]);

% Animation loop
for frame = 1:num_frames
    current_time = time_vector(frame);
    
    clf; % Clear figure
    hold on;
    
    % Collect current positions of the particles we're showing
    current_positions_x = [];
    current_positions_y = [];
    
    % For each trajectory we're showing
    for i = 1:num_trajectories_to_show
        
        % Find points up to current time
        time_mask = T{i} <= current_time;
        
        if sum(time_mask) > 0
            % Plot trajectory path up to current time
            plot(X{i}(time_mask), Y{i}(time_mask), 'b-', 'LineWidth', 1, 'Color', [0.6 0.6 0.9]);
            
            % Get current particle position
            current_x_traj = X{i}(time_mask);
            current_y_traj = Y{i}(time_mask);
            current_pos_x = current_x_traj(end);
            current_pos_y = current_y_traj(end);
            
            % Store current position for triangulation
            current_positions_x = [current_positions_x; current_pos_x];
            current_positions_y = [current_positions_y; current_pos_y];
        end
    end
    
    % Create triangulation on current positions (if we have enough particles)
    if length(current_positions_x) >= 3
        % Remove any duplicate positions
        positions = [current_positions_x, current_positions_y];
        [unique_positions, ~] = unique(positions, 'rows');
        triangulation_x = unique_positions(:, 1);
        triangulation_y = unique_positions(:, 2);
        
        % Make triangulation
        DT = delaunayTriangulation(triangulation_x, triangulation_y);
        triangles = DT.ConnectivityList;
        
        % Draw triangles FIRST (so they appear behind particles)
        for tri = 1:size(triangles, 1)
            triangle_x = triangulation_x(triangles(tri, [1,2,3,1]));
            triangle_y = triangulation_y(triangles(tri, [1,2,3,1]));
            plot(triangle_x, triangle_y, 'k-', 'LineWidth', 0.5, 'Color', [0.3, 0.3, 0.3]);
        end
    end
    
    % Draw current particle positions ON TOP of triangles
    if ~isempty(current_positions_x)
        scatter(current_positions_x, current_positions_y, 60, 'red', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
    end
    
    % Title and formatting
    num_particles_shown = length(current_positions_x);
    num_triangles_shown = 0;
    if length(current_positions_x) >= 3
        num_triangles_shown = size(triangles, 1);
    end
    
    title(sprintf('Trajectories + Triangulation - Time = %.1f\n%d particles, %d triangles', ...
          current_time, num_particles_shown, num_triangles_shown));
    xlabel('X position (pixels)');
    ylabel('Y position (pixels)');
    axis equal;
    grid on;
    
    % Keep axis limits fixed
    if frame == 1
        axis_limits = axis;
    else
        axis(axis_limits);
    end
    
    pause(0.1);
    
    if mod(frame, 20) == 0
        fprintf('Frame %d/%d (t = %.1f): %d particles, %d triangles\n', ...
                frame, num_frames, current_time, num_particles_shown, num_triangles_shown);
    end
end

fprintf('Animation complete!\n');

%%
%% TRAJECTORIES + TRIANGULATION ANIMATION - ALL PARTICLES
fprintf('Creating trajectory + triangulation animation for ALL particles...\n');

% Find time range
all_times = [];
for i = 1:num_tracks
    all_times = [all_times; T{i}];
end
time_min = min(all_times);
time_max = max(all_times);

fprintf('Time range: %.1f to %.1f\n', time_min, time_max);
fprintf('Showing ALL %d trajectories\n', num_tracks);

% Create time vector for animation
num_frames = 100;
time_vector = linspace(time_min, time_max, num_frames);

% Create figure
figure('Position', [100 100 800 600]);

% Animation loop
for frame = 1:num_frames
    current_time = time_vector(frame);
    
    clf; % Clear figure
    hold on;
    
    % Collect current positions of ALL particles
    current_positions_x = [];
    current_positions_y = [];
    
    % For EVERY trajectory
    for i = 1:num_tracks
        
        % Find points up to current time
        time_mask = T{i} <= current_time;
        
        if sum(time_mask) > 0
            % Plot trajectory path up to current time
            % plot(X{i}(time_mask), Y{i}(time_mask), 'b-', 'LineWidth', 0.8, 'Color', [0.6 0.6 0.9]);
            
            % Get current particle position
            current_x_traj = X{i}(time_mask);
            current_y_traj = Y{i}(time_mask);
            current_pos_x = current_x_traj(end);
            current_pos_y = current_y_traj(end);
            
            % Store current position for triangulation
            current_positions_x = [current_positions_x; current_pos_x];
            current_positions_y = [current_positions_y; current_pos_y];
        end
    end
    
    % Create triangulation on current positions (if we have enough particles)
    if length(current_positions_x) >= 3
        % Remove any duplicate positions
        positions = [current_positions_x, current_positions_y];
        [unique_positions, ~] = unique(positions, 'rows');
        triangulation_x = unique_positions(:, 1);
        triangulation_y = unique_positions(:, 2);
        
        % Make triangulation
        DT = delaunayTriangulation(triangulation_x, triangulation_y);
        triangles = DT.ConnectivityList;
        
        % Draw triangles FIRST (so they appear behind particles)
        for tri = 1:size(triangles, 1)
            triangle_x = triangulation_x(triangles(tri, [1,2,3,1]));
            triangle_y = triangulation_y(triangles(tri, [1,2,3,1]));
            plot(triangle_x, triangle_y, 'k-', 'LineWidth', 0.5, 'Color', [0.3, 0.3, 0.3]);
        end
    end
    
    % Draw current particle positions ON TOP of triangles
    if ~isempty(current_positions_x)
        scatter(current_positions_x, current_positions_y, 10, 'red', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
    end
    
    % Title and formatting
    num_particles_shown = length(current_positions_x);
    num_triangles_shown = 0;
    if length(current_positions_x) >= 3
        num_triangles_shown = size(triangles, 1);
    end
    
    title(sprintf('ALL Trajectories + Triangulation - Time = %.1f\n%d particles, %d triangles', ...
          current_time, num_particles_shown, num_triangles_shown));
    xlabel('X position (pixels)');
    ylabel('Y position (pixels)');
    axis equal;
    grid on;
    
    % Keep axis limits fixed
    if frame == 1
        axis_limits = axis;
    else
        axis(axis_limits);
    end
    
    pause(0.1);
    
    if mod(frame, 20) == 0
        fprintf('Frame %d/%d (t = %.1f): %d particles, %d triangles\n', ...
                frame, num_frames, current_time, num_particles_shown, num_triangles_shown);
    end
end

fprintf('Animation complete!\n');