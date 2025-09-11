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
    
    time_step = 5; 
    analysis_times = unique_times(1:time_step:end);
    
    %% Parameters
    calib = 1.0;
    distance_threshold = 5;
    
    all_x = []; all_y = [];
    for i = 1:length(X)
        all_x = [all_x; X{i}];
        all_y = [all_y; Y{i}];
    end
    all_x = all_x * calib;
    all_y = all_y * calib;
    
    x_padding = 0.05 * (max(all_x) - min(all_x));
    y_padding = 0.05 * (max(all_y) - min(all_y));
    x_limits = [min(all_x) - x_padding, max(all_x) + x_padding];
    y_limits = [min(all_y) - y_padding, max(all_y) + y_padding];
    
    %%
    figure('Position', [100, 100, 1200, 500]);
    
    % Animation loop
    for t_idx = 1:length(analysis_times)
        current_time = analysis_times(t_idx);
        
        % Extract positions at this time point
        positions = [];
        particle_ids = [];
        
        for i = 1:length(X)
            timeIdx = find(T{i} == current_time, 1);
            if ~isempty(timeIdx)
                positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
                particle_ids = [particle_ids; i];
            end
        end
        
        positions = positions * calib;
        N_particles = size(positions, 1);
        
        if N_particles < 3
            continue;
        end
        
        % Create Delaunay triangulation as ben code
        DT = delaunayTriangulation(positions(:,1), positions(:,2));
        triangles = DT.ConnectivityList;
        n_triangles = size(triangles, 1);
        
        valid_triangles = [];
        for i = 1:n_triangles
            p1 = triangles(i, 1);
            p2 = triangles(i, 2);
            p3 = triangles(i, 3);
            
            % Calculate the three side lengths
            side1 = norm(positions(p1,:) - positions(p2,:));
            side2 = norm(positions(p2,:) - positions(p3,:));
            side3 = norm(positions(p3,:) - positions(p1,:));
            
            max_side = max([side1, side2, side3]);
            
            % Keep triangle only if  sides are below threshold
            if max_side <= distance_threshold
                valid_triangles = [valid_triangles; triangles(i,:)];
            end
        end
        
    
    
        
        n_valid_triangles = size(valid_triangles, 1);
            clf;
        
        subplot(1, 2, 1);
        triplot(DT, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5);
        hold on;
        scatter(positions(:,1), positions(:,2), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.8);
        axis equal; grid on;
        xlim(x_limits); ylim(y_limits);  % Fix axis limits
        % title(sprintf('Full Delaunay Triangulation\n(%d triangles)', n_triangles));
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
        xlim(x_limits); ylim(y_limits);  % Fix axis limits
        % title(sprintf('Distance-Filtered Triangles\n(%d triangles, threshold = %.1f μm)', ...
        %     n_valid_triangles, distance_threshold));
        xlabel('x (μm)'); ylabel('y (μm)');
        
        % sgtitle(sprintf('Time = %d | Particles = %d | Valid/Total = %d/%d', ...
        %     current_time, N_particles, n_valid_triangles, n_triangles));
            drawnow;
        % pause(0.1);  % 
    end