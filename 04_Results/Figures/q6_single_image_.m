clear all; clc; close all;

%% Calibration
pixels_per_micron = 5.5;  % 15.5 pixels/μm
particle_diameter_um = 1.4;  % μm
particle_diameter_px = particle_diameter_um * pixels_per_micron;
particle_radius_px = particle_diameter_px / 2;

%% Load the last image
img_file = 'ps_2um_20x_b100_10fps_5_croped.jpg';
img = imread(img_file);

fprintf('Loaded image: %s\n', img_file);
fprintf('Image size: %d x %d pixels\n', size(img, 1), size(img, 2));

%% Detect particles using circular Hough transform
fprintf('\nDetecting particles...\n');

% Parameters for circle detection
radius_range = round([particle_radius_px*0.7, particle_radius_px*1.3]);
sensitivity = 0.955;
edge_threshold = 0.1;

[centers_px, radii_px] = imfindcircles(img, radius_range, ...
    'ObjectPolarity', 'bright', ...
    'Sensitivity', sensitivity, ...
    'EdgeThreshold', edge_threshold);

fprintf('Detected %d particles\n', length(radii_px));

if isempty(centers_px)
    error('No particles detected! Try adjusting sensitivity or edge_threshold');
end

%% Convert positions to microns for Q6 calculation
positions_um = centers_px / pixels_per_micron;

%% Compute Delaunay triangulation
DT = delaunayTriangulation(positions_um(:,1), positions_um(:,2));
nPoints = size(DT.Points, 1);

%% Calculate hexatic order parameter Q6
fprintf('\nCalculating Q6...\n');
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

% Calculate statistics
fprintf('\nQ6 statistics:\n');
fprintf('  Mean: %.3f\n', mean(Q6));
fprintf('  Std:  %.3f\n', std(Q6));
fprintf('  Max:  %.3f\n', max(Q6));
fprintf('  Min:  %.3f\n', min(Q6));

%% Map Q6 values to RGB colors manually
% Create internet colormap
cmap = internet(256);

% Normalize Q6 to [0, 1] range (already in this range)
Q6_normalized = Q6;
Q6_normalized(Q6_normalized < 0) = 0;
Q6_normalized(Q6_normalized > 1) = 1;

% Map to colormap indices
indices = round(Q6_normalized * 255) + 1;  % 1 to 256
indices(indices < 1) = 1;
indices(indices > 256) = 256;

% Get RGB colors for each particle
colors = cmap(indices, :);

%% Create figure with Q6 overlay
fig = figure('Position', [100, 100, 1200, 900], 'Color', 'white');

% Show grayscale image
imshow(img);
hold on;

% Overlay particles with explicit RGB colors
for i = 1:nPoints
    scatter(centers_px(i,1), centers_px(i,2), 110, colors(i,:), 'filled', ...
        'MarkerEdgeColor', 'none', 'LineWidth', 2);
end

% Create manual colorbar
ax_main = gca;
ax_cb = axes('Position', [0.92, 0.3, 0.02, 0.4]);
imagesc(ax_cb, linspace(0, 1, 256)');
ax_cb.YDir = 'normal';
colormap(ax_cb, internet);
ax_cb.XTick = [];
ax_cb.YTick = [1, 64, 128, 192, 256];
ax_cb.YTickLabel = {'0.0', '0.25', '0.5', '0.75', '1.0'};
ax_cb.TickLabelInterpreter = 'latex';
ax_cb.FontSize = 14;
ylabel(ax_cb, '$\Psi_6$', 'Interpreter', 'latex', 'FontSize', 18);

% Return to main axes
axes(ax_main);

% Add scale bar
[img_height, img_width] = size(img);
scale_bar_length_um = 10;
scale_bar_length_px = scale_bar_length_um * pixels_per_micron;
x_bar = [50, 50 + scale_bar_length_px];
y_bar = [img_height - 30, img_height - 30];
plot(x_bar, y_bar, 'w-', 'LineWidth', 5);
text(50 + scale_bar_length_px/2, img_height - 50, ...
    sprintf('%d $\\mu$m', scale_bar_length_um), ...
    'Color', 'w', 'FontSize', 16, 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Title
% title(sprintf('Hexatic Order Parameter (N=%d, $\\langle Q_6 \\rangle$=%.2f)', ...
%     nPoints, mean(Q6)), ...
%     'FontSize', 18, 'Interpreter', 'latex');

hold off;

% Save
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 10, 8]);
export_fig( 'Q6_dots_overlay.pdf');
% 
% fprintf('\nFigure saved as Q6_dots_overlay.png\n');


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