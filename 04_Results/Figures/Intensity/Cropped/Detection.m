clear all; clc; close all;

%% Calibration
pixels_per_micron = 16.0;  % 16 pixels/μm for 60x
particle_diameter_um = 1.90;  % μm
particle_diameter_px = particle_diameter_um * pixels_per_micron;
particle_radius_px = particle_diameter_px / 2;

%% Image files
img_files = {'exp10000.tif', 'exp10001.tif', 'exp10002.tif', 'exp10003.tif'};
n_images = length(img_files);

%% Create figure
fig = figure('Position', [100, 100, 1600, 400], 'Color', 'white');

%% Process each image
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
    fprintf('Radius range searched: [%.1f, %.1f] pixels\n', radius_range(1), radius_range(2));
    
    %% Plot detection result
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
    title(sprintf('Image %d: %d particles', img_idx, length(radii_px)), ...
        'FontSize', 14, 'Color', 'w');
    
    hold off;
end

%% Save
exportgraphics(gcf, 'detection_check.png', 'Resolution', 300);
fprintf('\n\nFigure saved as detection_check.png\n');
fprintf('\nRed circles = detected particles, Green dots = centers\n');