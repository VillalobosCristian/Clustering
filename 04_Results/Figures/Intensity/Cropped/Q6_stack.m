clear all; clc; close all;

%% Calibration
pixels_per_micron = 11.5;  % 16 pixels/μm for 60x
particle_diameter_um = 1.90;  % μm
particle_diameter_px = particle_diameter_um * pixels_per_micron;
particle_radius_px = particle_diameter_px / 2;

%% Image files (4 different intensities)
img_files = {'exp10000.tif', 'exp10001.tif', 'exp10002.tif', 'exp10003.tif'};
n_images = length(img_files);

% You can specify intensities/powers for labels
powers = [20, 40, 60, 100];  % Adjust these values to match your experiments

%% Process each image and save separately
for img_idx = 1:n_images
    
    fprintf('\n=== Processing %s ===\n', img_files{img_idx});
    
    %% Load image
    img = imread(img_files{img_idx});
    fprintf('Image size: %d x %d pixels\n', size(img, 1), size(img, 2));
    
    %% Detect particles
    fprintf('Detecting particles...\n');
    radius_range = round([particle_radius_px*0.7, particle_radius_px*1.3]);
    sensitivity = 0.92;
    edge_threshold = 0.1;
    
    [centers_px, radii_px] = imfindcircles(img, radius_range, ...
        'ObjectPolarity', 'bright', ...
        'Sensitivity', sensitivity, ...
        'EdgeThreshold', edge_threshold);
    
    fprintf('Detected %d particles\n', length(radii_px));
    
    if isempty(centers_px)
        warning('No particles detected in %s!', img_files{img_idx});
        continue;
    end
    
    %% Convert to microns
    positions_um = centers_px / pixels_per_micron;
    
    %% Compute Delaunay triangulation
    DT = delaunayTriangulation(positions_um(:,1), positions_um(:,2));
    nPoints = size(DT.Points, 1);
    
    %% Calculate Q6
    fprintf('Calculating Q6...\n');
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
    axis image;  % Maintain aspect ratio
    axis off;    % Hide axes
    
    hold on;
    
    % Plot particles with Q6 colors
    for i = 1:nPoints
        scatter(centers_px(i,1), centers_px(i,2), 100, colors(i,:), 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
    end
    
    % Add scale bar ONLY to the last image (img_idx == n_images)
    [img_height, img_width] = size(img);
    if img_idx == n_images
        scale_bar_length_um = 10;
        scale_bar_length_px = scale_bar_length_um * pixels_per_micron;
        x_bar = [50, 50 + scale_bar_length_px];
        y_bar = [img_height - 30, img_height - 30];
        plot(x_bar, y_bar, 'w-', 'LineWidth', 4);
        text(50 + scale_bar_length_px/2, img_height - 50, ...
            sprintf('%d $\\mu$m', scale_bar_length_um), ...
            'Color', 'w', 'FontSize', 14, 'Interpreter', 'latex', ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    % Title with power and Q6 stats
    title(sprintf('%d\\%% ($\\langle Q_6 \\rangle$ = %.2f)', ...
        powers(img_idx), mean(Q6)), ...
        'FontSize', 18, 'Interpreter', 'latex', 'Color', 'k');
    
    hold off;
    
    % Save individual image using export_fig
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [1, 1, 6, 6]);
    output_filename = sprintf('Q6_power_%d_percent', powers(img_idx));
    export_fig(output_filename, '-pdf', '-painters', '-r300', '-transparent');
    fprintf('Saved: %s.pdf\n', output_filename);
    
    close(gcf);  % Close figure after saving
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