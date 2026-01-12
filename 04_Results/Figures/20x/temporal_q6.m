clear all; clc; close all;

%% Calibration
pixels_per_micron = 5.5;  % 5.5 pixels/μm
particle_diameter_um = 1.4;  % μm
particle_diameter_px = particle_diameter_um * pixels_per_micron;
particle_radius_px = particle_diameter_px / 2;

%% Image files (4 different time steps)
img_files = {'1.jpg', '2.jpg', '3.jpg', '4.jpg'};
n_images = length(img_files);

% Time values in seconds
times = [0, 60, 300, 600];  % seconds

%% Create detection check figure first
fig_check = figure('Position', [100, 100, 1600, 400], 'Color', 'white');

for img_idx = 1:n_images
    
    fprintf('\n=== Processing %s ===\n', img_files{img_idx});
    
    %% Load image
    img = imread(img_files{img_idx});
    fprintf('Image size: %d x %d pixels\n', size(img, 1), size(img, 2));
    
    %% Detect particles
    fprintf('Detecting particles...\n');
    radius_range = round([particle_radius_px*0.7, particle_radius_px*1.3]);
    sensitivity = 0.95;
    edge_threshold = 0.1;
    
    [centers_px, radii_px] = imfindcircles(img, radius_range, ...
        'ObjectPolarity', 'bright', ...
        'Sensitivity', sensitivity, ...
        'EdgeThreshold', edge_threshold);
    
    fprintf('Detected %d particles\n', length(radii_px));
    fprintf('Radius range searched: [%.1f, %.1f] pixels\n', radius_range(1), radius_range(2));
    
    %% Plot detection result in subplot
    figure(fig_check);
    subplot(1, n_images, img_idx);
    imshow(img);
    hold on;
    
    % Draw circles around detected particles
    if ~isempty(centers_px)
        viscircles(centers_px, radii_px, 'Color', 'r', 'LineWidth', 1.5);
        
        % Also mark centers
        scatter(centers_px(:,1), centers_px(:,2), 30, 'g', 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1);
    end
    
    % Title
    title(sprintf('t = %d s (%d particles)', times(img_idx), length(radii_px)), ...
        'FontSize', 14, 'Color', 'w');
    
    hold off;
end

%% Save detection check figure
figure(fig_check);
export_fig('detection_check_time_series', '-png', '-r300');
fprintf('\n\nDetection check saved as detection_check_time_series.png\n');
fprintf('Red circles = detected particles, Green dots = centers\n\n');

%% Now process each image for Q6 analysis
for img_idx = 1:n_images
    
    fprintf('\n=== Creating Q6 figure for t=%d s ===\n', times(img_idx));
    
    %% Load image
    img = imread(img_files{img_idx});
    
    %% Detect particles (re-run detection)
    radius_range = round([particle_radius_px*0.7, particle_radius_px*1.3]);
    sensitivity = 0.96;
    edge_threshold = 0.1;
    
    [centers_px, radii_px] = imfindcircles(img, radius_range, ...
        'ObjectPolarity', 'bright', ...
        'Sensitivity', sensitivity, ...
        'EdgeThreshold', edge_threshold);
    
    if isempty(centers_px)
        warning('No particles detected in %s!', img_files{img_idx});
        
        % Create figure with just the image (no particles)
        fig = figure('Position', [100, 100, 800, 800], 'Color', 'white');
        imagesc(img); 
        colormap(gca, gray);
        axis image;
        axis off;
        
        % Title with time
        title(sprintf('$t = %d$ s', times(img_idx)), ...
            'FontSize', 18, 'Interpreter', 'latex', 'Color', 'k');
        
        % Save
        set(gcf, 'Units', 'inches');
        set(gcf, 'Position', [1, 1, 6, 6]);
        output_filename = sprintf('Q6_time_%d_s', times(img_idx));
        export_fig(output_filename, '-pdf', '-painters', '-r300', '-transparent');
        fprintf('Saved: %s.pdf\n', output_filename);
        close(gcf);
        
        continue;
    end
    
    %% Convert to microns
    positions_um = centers_px / pixels_per_micron;
    
    %% Compute Delaunay triangulation
    DT = delaunayTriangulation(positions_um(:,1), positions_um(:,2));
    nPoints = size(DT.Points, 1);
    
    %% Calculate Q6
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
    
    fprintf('Q6: Mean=%.3f, Std=%.3f, Max=%.3f\n', ...
        mean(Q6), std(Q6), max(Q6));
    
    %% Map Q6 to colors
    cmap = internet(256);
    Q6_normalized = max(0, min(1, Q6));
    indices = round(Q6_normalized * 255) + 1;
    indices = max(1, min(256, indices));
    colors = cmap(indices, :);
    
    %% Create individual figure for this image
    fig = figure('Position', [100, 100, 800, 800], 'Color', 'white');
    
    % Display image with imagesc and gray colormap
    imagesc(img); 
    colormap(gca, gray);
    axis image;
    axis off;
    
    hold on;
    
    % Plot particles with Q6 colors
    for i = 1:nPoints
        scatter(centers_px(i,1), centers_px(i,2), 10, colors(i,:), 'filled', ...
            'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
    end
    
    % Add scale bar ONLY to the last image
    [img_height, img_width] = size(img);
    if img_idx == n_images
        scale_bar_length_um = 20;
        scale_bar_length_px = scale_bar_length_um * pixels_per_micron;
        x_bar = [50, 50 + scale_bar_length_px];
        y_bar = [img_height - 30, img_height - 30];
        plot(x_bar, y_bar, 'w-', 'LineWidth', 4);
        text(50 + scale_bar_length_px/2, img_height - 50, ...
            sprintf('%d $\\mu$m', scale_bar_length_um), ...
            'Color', 'w', 'FontSize', 14, 'Interpreter', 'latex', ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    % Title with time and Q6 stats
    title(sprintf('$t = %d$ s ($\\langle Q_6 \\rangle$ = %.2f)', ...
        times(img_idx), mean(Q6)), ...
        'FontSize', 18, 'Interpreter', 'latex', 'Color', 'k');
    
    hold off;
    
    % Save individual image using export_fig
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [1, 1, 6, 6]);
    output_filename = sprintf('Q6_time_%d_s', times(img_idx));
    export_fig(output_filename, '-pdf', '-painters', '-r300', '-transparent');
    fprintf('Saved: %s.pdf\n', output_filename);
    
    close(gcf);
end

%% Create separate colorbar figure
fig_cb = figure('Position', [100, 100, 200, 600], 'Color', 'white');
ax_cb = axes('Position', [0.3, 0.1, 0.3, 0.8]);
imagesc(ax_cb, linspace(0, 1, 256)');
ax_cb.YDir = 'normal';
colormap(ax_cb, internet);
ax_cb.XTick = [];
ax_cb.YTick = [1, 64, 128, 192, 256];
ax_cb.YTickLabel = {'0.0', '0.25', '0.5', '0.75', '1.0'};
ax_cb.TickLabelInterpreter = 'latex';
ax_cb.FontSize = 16;
ylabel(ax_cb, '$Q_6$', 'Interpreter', 'latex', 'FontSize', 20);

% Save colorbar using export_fig
export_fig('Q6_colorbar', '-pdf', '-painters', '-r300', '-transparent');
fprintf('\nSaved: Q6_colorbar.pdf\n');

fprintf('\n=== All images saved! ===\n');

%% Helper function
function connected = isConnected(DT, i, j)
    triangles = vertexAttachments(DT, i);
    triangles = triangles{1};
    connected = false;
    for t = 1:length(triangles)
        if any(DT.ConnectivityList(triangles(t), :) == j)
            connected = true;
            break;
        end
    end
end