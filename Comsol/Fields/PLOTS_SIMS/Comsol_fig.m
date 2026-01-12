clear all; clc; close all;

%% FIGURE PARAMETERS
FIGURE_WIDTH = 8;
FIGURE_HEIGHT = 5;
FONT_SIZE_LABELS = 28;
FONT_SIZE_AXES = 28;
FONT_SIZE_LEGEND = 14;
LINE_WIDTH = 2.5;
AXES_WIDTH = 1.5;

%% SETUP
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');

%% LOAD COMSOL DATA
[file_name, path_name] = uigetfile('*.csv', 'Select COMSOL CSV file');
if isequal(file_name, 0), error('No file selected'); end

filename = fullfile(path_name, file_name);
fprintf('Loading: %s\n', file_name);

data = readmatrix(filename, 'NumHeaderLines', 9);

% Extract and clean data immediately
r = real(data(:, 1));     z = real(data(:, 2));     T = real(data(:, 3));
F_TPr = real(data(:, 4)); F_Dr = real(data(:, 5));  F_TPz = real(data(:, 6));
F_Dz = real(data(:, 7));  F_g = real(data(:, 8));   u = real(data(:, 9));     
w = real(data(:, 10));

fprintf('Data points: %d\n', length(r));
fprintf('Domain: r = 0-%.0f μm, z = 0-%.0f μm\n', max(r)*1e6, max(z)*1e6);

% Calculate total forces
F_total_r = F_TPr + F_Dr;
F_total_z = F_TPz + F_Dz + F_g;

% Get illumination radius
prompt = {'Illumination radius R_{illum} (μm):'};
dlgtitle = 'Experimental Parameter';
dims = [1 50];
defaultans = {'75'};
answer = inputdlg(prompt, dlgtitle, dims, defaultans);
R_illu = str2double(answer{1});

fprintf('R_illum = %.0f μm\n', R_illu);
fprintf('Temperature: %.2f - %.2f K (ΔT = %.2f K)\n\n', ...
    min(T), max(T), max(T)-min(T));

%% OUTPUT DIRECTORY
output_dir = fullfile(path_name, 'figures_comsol_analysis');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% CREATE HIGH-RESOLUTION INTERPOLATION GRID
grid_resolution = 300;  % Reduced from 400 for speed
r_grid = linspace(0, max(r), grid_resolution);
z_grid = linspace(0, max(z), grid_resolution);
[R_grid, Z_grid] = meshgrid(r_grid, z_grid);

fprintf('Interpolating fields...\n');
warning('off', 'MATLAB:griddata:DuplicateDataPoints');
T_grid = griddata(r, z, T, R_grid, Z_grid, 'cubic');
u_grid = griddata(r, z, u, R_grid, Z_grid, 'cubic');
w_grid = griddata(r, z, w, R_grid, Z_grid, 'cubic');
F_TPr_grid = griddata(r, z, F_TPr, R_grid, Z_grid, 'cubic');
F_Dr_grid = griddata(r, z, F_Dr, R_grid, Z_grid, 'cubic');
F_total_r_grid = griddata(r, z, F_total_r, R_grid, Z_grid, 'cubic');
F_total_z_grid = griddata(r, z, F_total_z, R_grid, Z_grid, 'cubic');
warning('on', 'MATLAB:griddata:DuplicateDataPoints');

%% COMPUTE TEMPERATURE GRADIENTS
fprintf('Computing temperature gradients...\n');

% Calculate ∂T/∂r and ∂T/∂z
[dT_dz, dT_dr] = gradient(T_grid, z_grid, r_grid);

% Calculate gradient magnitude |∇T|
grad_T_mag = sqrt(dT_dr.^2 + dT_dz.^2);

% Convert to K/μm
dT_dr_um = dT_dr * 1e-6;
dT_dz_um = dT_dz * 1e-6;
grad_T_mag_um = grad_T_mag * 1e-6;

fprintf('Gradient statistics:\n');
fprintf('  |∇T| max = %.4f K/μm\n', max(grad_T_mag_um(:), [], 'omitnan'));
fprintf('  |∇T| mean = %.4f K/μm\n', mean(grad_T_mag_um(:), 'omitnan'));
fprintf('  ∂T/∂r max = %.4f K/μm\n', max(abs(dT_dr_um(:)), [], 'omitnan'));
fprintf('  ∂T/∂z max = %.4f K/μm\n\n', max(abs(dT_dz_um(:)), [], 'omitnan'));

%% EXTRACT NEAR-SURFACE LAYERS
fprintf('Extracting layers...\n');

z_layers_um = [0,3, 5, 10];  % 6 well-spaced layers
z_tolerance_um = 2;  % ±2 μm
n_layers = length(z_layers_um);
colors_layers = turbo(n_layers);

layer_data = cell(n_layers, 1);

for i_layer = 1:n_layers
    z_target = z_layers_um(i_layer) * 1e-6;
    z_tol = z_tolerance_um * 1e-6;
    idx_layer = abs(z - z_target) < z_tol;
    
    if sum(idx_layer) > 10
        layer_data{i_layer}.r = r(idx_layer);
        layer_data{i_layer}.T = T(idx_layer);
        layer_data{i_layer}.u = u(idx_layer);
        layer_data{i_layer}.w = w(idx_layer);
        layer_data{i_layer}.F_r = F_total_r(idx_layer);
        layer_data{i_layer}.F_TPr = F_TPr(idx_layer);
        layer_data{i_layer}.F_Dr = F_Dr(idx_layer);
        layer_data{i_layer}.z_actual = mean(z(idx_layer))*1e6;
        layer_data{i_layer}.n_points = sum(idx_layer);
        
        fprintf('  z=%.0fμm: %d points\n', layer_data{i_layer}.z_actual, ...
                layer_data{i_layer}.n_points);
    end
end

%% COMPUTE RADIAL PROFILES
n_bins = 25;  % Good resolution
r_max = max(r);
r_edges = linspace(0, r_max, n_bins+1);
r_centers = (r_edges(1:end-1) + r_edges(2:end))/2 * 1e6;
r_centers = r_centers(:);

profiles = struct();
profiles.r_centers = r_centers;

for i_layer = 1:n_layers
    if isempty(layer_data{i_layer}), continue; end
    
    data_layer = layer_data{i_layer};
    r_layer = data_layer.r;
    u_layer = data_layer.u;
    F_layer = data_layer.F_r;
    F_TP_layer = data_layer.F_TPr;
    F_D_layer = data_layer.F_Dr;
    
    u_mean = nan(n_bins, 1);  u_std = nan(n_bins, 1);
    F_mean = nan(n_bins, 1);  F_std = nan(n_bins, 1);
    F_TP_mean = nan(n_bins, 1);
    F_D_mean = nan(n_bins, 1);
    
    for i_bin = 1:n_bins
        idx_bin = r_layer >= r_edges(i_bin) & r_layer < r_edges(i_bin+1);
        if sum(idx_bin) >= 2
            u_mean(i_bin) = mean(u_layer(idx_bin)) * 1e6;
            u_std(i_bin) = std(u_layer(idx_bin)) * 1e6;
            F_mean(i_bin) = mean(F_layer(idx_bin)) * 1e15;
            F_std(i_bin) = std(F_layer(idx_bin)) * 1e15;
            F_TP_mean(i_bin) = mean(F_TP_layer(idx_bin)) * 1e15;
            F_D_mean(i_bin) = mean(F_D_layer(idx_bin)) * 1e15;
        end
    end
    
    % Light smoothing
    u_mean = movmean(u_mean, 2, 'omitnan');
    F_mean = movmean(F_mean, 2, 'omitnan');
    F_TP_mean = movmean(F_TP_mean, 2, 'omitnan');
    F_D_mean = movmean(F_D_mean, 2, 'omitnan');
    
    % Consistent field names: z000, z050, z100, etc.
    field_name = sprintf('z%03d', round(z_layers_um(i_layer)*10));
    profiles.(field_name).u_mean = u_mean;
    profiles.(field_name).u_std = u_std;
    profiles.(field_name).F_mean = F_mean;
    profiles.(field_name).F_std = F_std;
    profiles.(field_name).F_TP_mean = F_TP_mean;
    profiles.(field_name).F_D_mean = F_D_mean;
    profiles.(field_name).z_actual = data_layer.z_actual;
end

%% KEY METRICS
Delta_T = max(T) - min(T);
fprintf('\n=== Key Results ===\n');
fprintf('ΔT = %.2f K\n', Delta_T);

for i_layer = 1:n_layers
    field_name = sprintf('z%03d', round(z_layers_um(i_layer)*10));
    if ~isfield(profiles, field_name), continue; end
    
    prof = profiles.(field_name);
    valid_u = ~isnan(prof.u_mean);
    
    if sum(valid_u) > 0
        u_vals = prof.u_mean(valid_u);
        min_u = min(u_vals);
        if min_u < -0.01
            fprintf('z=%.0fμm: INWARD v_r = %.3f μm/s\n', prof.z_actual, min_u);
        end
    end
end

%% ========================================================================
%% FIGURE 1: 2D Temperature Field with Flow
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');

contourf(R_grid*1e6, Z_grid*1e6, T_grid, 150, 'LineStyle', 'none');
% colormap(flipud(spectral));
colormap(thermal)
hold on;

skip_vec = 15;
quiver(R_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       Z_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       u_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       w_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       1, 'k', 'LineWidth', 1,'ShowArrowHead','on');

xline(R_illu, 'r--', 'LineWidth', 2);
plot([0, max(r)*1e6], [30, 30], 'g--', 'LineWidth', 2.5);

cb = colorbar;
cb.Label.String = 'T (K)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = FONT_SIZE_LABELS;
cb.FontSize = FONT_SIZE_AXES;
cb.TickLabelInterpreter='latex';
xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);

% xlim([0, min(max(r)*1e6, R_illu*3)]);
% ylim([0, 50]);
 ylim([0, 120]);

box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig1_Temperature_Flow');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');

%% ========================================================================
%% FIGURE 2: Temperature Gradient + Flow Field
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Plot gradient magnitude with proper color scaling
% pcolor(R_grid*1e6, Z_grid*1e6, grad_T_mag_um);
contourf(R_grid*1e6, Z_grid*1e6, grad_T_mag_um, 150, 'LineStyle', 'none');

shading interp;
hold on;

% Overlay velocity field
skip_vec = 15;
quiver(R_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       Z_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       u_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       w_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       1, 'w', 'LineWidth', 1, 'MaxHeadSize', 1);

xline(R_illu, 'r--', 'LineWidth', 2);
% plot([0, max(r_grid)*1e6], [30, 30], 'g--', 'LineWidth', 2.5);

cb = colorbar;
cb.Label.String = '$|\nabla T|$ (K/$\mu$m)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = FONT_SIZE_LABELS;
cb.FontSize = FONT_SIZE_AXES;
cb.TickLabelInterpreter='latex'

% Set proper colormap and limits
colormap(thermal);
% Fix color scale - gradient is small!
caxis([0, max(grad_T_mag_um(:), [], 'omitnan')]);

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);

% xlim([0, min(max(r_grid)*1e6, R_illu*3)]);
% ylim([0, 50]);
box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig2_Gradient_Flow');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');

%% ========================================================================
%% FIGURE 3: Surface Temperature Gradient (1D)
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Extract surface data
z_surf_tolerance = 0.5e-6;
surf_idx = abs(z) < z_surf_tolerance;
r_surf = r(surf_idx) * 1e6;
T_surf = T(surf_idx);

% Bin temperature
n_bins_surf = 50;
r_edges_surf = linspace(0, max(r), n_bins_surf+1);
r_centers_surf = (r_edges_surf(1:end-1) + r_edges_surf(2:end))/2 * 1e6;

T_mean_surf = nan(n_bins_surf, 1);
for i = 1:n_bins_surf
    idx = r_surf >= r_edges_surf(i)*1e6 & r_surf < r_edges_surf(i+1)*1e6;
    if sum(idx) >= 2
        T_mean_surf(i) = mean(T_surf(idx));
    end
end

T_mean_surf = movmean(T_mean_surf, 2, 'omitnan');
valid_T_surf = ~isnan(T_mean_surf);

% Compute gradient
grad_T_r_surf = gradient(T_mean_surf(valid_T_surf), r_centers_surf(valid_T_surf));

% Plot
plot(r_centers_surf(valid_T_surf), grad_T_r_surf*1e-6, 'b-', 'LineWidth', LINE_WIDTH);
hold on;
yline(0, 'k--', 'LineWidth', 1.5);
xline(R_illu, 'r--', 'LineWidth', 2);

% % Mark minimum (steepest negative)
% [min_grad, min_idx] = min(grad_T_r_surf);
% plot(r_centers_surf(valid_T_surf(min_idx)), min_grad*1e-6, 'ro', ...
%      'MarkerSize', 15, 'MarkerFaceColor', 'r', 'LineWidth', 2);


xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$\partial T/\partial r$ (K/$\mu$m)', 'FontSize', FONT_SIZE_LABELS);

grid on; box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig3_Surface_Gradient');
% export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
% fprintf('Saved: Fig3\n');

%% ========================================================================
%% FIGURE 4: Velocity Profiles at Multiple Heights
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');
hold on;
colors_layers = viridis(n_layers);

for i_layer = 2:n_layers
    field_name = sprintf('z%03d', round(z_layers_um(i_layer)*10));
    if ~isfield(profiles, field_name), continue; end
    
    prof = profiles.(field_name);
    valid = ~isnan(prof.u_mean);
    
    if sum(valid) > 0
        plot(profiles.r_centers(valid), prof.u_mean(valid), '-', ...
             'Color', colors_layers(i_layer,:), 'LineWidth', LINE_WIDTH, ...
             'DisplayName', sprintf('$z = %.0f$ $\\mu$m', prof.z_actual));
    end
end

yline(0, 'k--', 'LineWidth', 1.5);
xline(R_illu, 'r--', 'LineWidth', 1.5, 'Alpha', 0.5);

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$v_r$ ($\mu$m/s)', 'FontSize', FONT_SIZE_LABELS);
% legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'NumColumns', 2);
grid off; box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 5, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'F2ig4_Velocity_vs_Height');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');

%% ========================================================================
%% FIGURE 5: Force Balance Components

%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');
hold on;

% Focus on a specific height (z=15 μm - particle observation)
field_name = 'z150';
if isfield(profiles, field_name)
    prof = profiles.(field_name);
    valid = ~isnan(prof.F_TP_mean) & ~isnan(prof.F_D_mean);
    
    if sum(valid) > 0
        % Plot force components
        plot(profiles.r_centers(valid), prof.F_TP_mean(valid), 'r-', ...
             'LineWidth', LINE_WIDTH, 'DisplayName', 'Thermophoresis $F_{TP}$');
        plot(profiles.r_centers(valid), prof.F_D_mean(valid), 'b-', ...
             'LineWidth', LINE_WIDTH, 'DisplayName', 'Drag $F_D$ (convection)');
        plot(profiles.r_centers(valid), prof.F_mean(valid), 'k-', ...
             'LineWidth', LINE_WIDTH+0.5, 'DisplayName', 'Total $F_{total}$');
    end
end

yline(0, 'k--', 'LineWidth', 1.5);
% yline(4.3, 'k:', 'LineWidth', 1.5, 'DisplayName', '$k_BT/a$');
xline(R_illu, 'r--', 'LineWidth', 1.5, 'Alpha', 0.5, 'HandleVisibility', 'off');

% Add annotations explaining force directions
y_range = ylim;
% text(0.05, 0.95, '$F_{TP} > 0$ (outward)', 'Units', 'normalized', ...
%      'FontSize', FONT_SIZE_LEGEND-2, 'Color', 'r', 'VerticalAlignment', 'top');
% text(0.05, 0.88, '$F_D < 0$ (inward)', 'Units', 'normalized', ...
%      'FontSize', FONT_SIZE_LEGEND-2, 'Color', 'b', 'VerticalAlignment', 'top');
% text(0.05, 0.81, '$F_{total} < 0$ (net inward)', 'Units', 'normalized', ...
%      'FontSize', FONT_SIZE_LEGEND-2, 'Color', 'k', 'VerticalAlignment', 'top');

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('Force (fN)', 'FontSize', FONT_SIZE_LABELS);
% title(sprintf('Force Balance at $z = %.0f$ $\\mu$m', profiles.z150.z_actual), ...
      % 'FontSize', FONT_SIZE_LABELS-2);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-2);
% grid on;
box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH/2 ...
    , FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig5_Force_Balance');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');

%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');
hold on;

for i_layer = 1:n_layers
    field_name = sprintf('z%03d', round(z_layers_um(i_layer)*10));
    if ~isfield(profiles, field_name), continue; end
    
    prof = profiles.(field_name);
    valid = ~isnan(prof.F_mean);
    
    if sum(valid) > 0
        plot(profiles.r_centers(valid), prof.F_mean(valid), '-', ...
             'Color', colors_layers(i_layer,:), 'LineWidth', LINE_WIDTH, ...
             'DisplayName', sprintf('$z = %.0f$ $\\mu$m', prof.z_actual));
    end
end

yline(0, 'k--', 'LineWidth', 1.5);
yline(4.3, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
text(max(profiles.r_centers)*0.9, 4.3, '$k_BT/a$', 'FontSize', FONT_SIZE_LEGEND);
xline(R_illu, 'r--', 'LineWidth', 1.5, 'Alpha', 0.5, 'HandleVisibility', 'off');

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$F_r$ (fN)', 'FontSize', FONT_SIZE_LABELS);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-4, 'NumColumns', 2);
grid on; box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig6_Force_vs_Height');
% export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');

%% ========================================================================
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');

contourf(R_grid*1e6, Z_grid*1e6, T_grid, 50, 'LineStyle', 'none');
colormap(thermal);
hold on;

% Overlay force vectors (sparse)
skip_force = 8;
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_total_r_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_total_z_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       2.5, 'w', 'LineWidth', 1, 'MaxHeadSize', 0.8);

xline(R_illu, 'r--', 'LineWidth', 2);
plot([0, max(r)*1e6], [30, 30], 'g--', 'LineWidth', 2.5);

cb = colorbar;
cb.Label.String = 'Temperature (K)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = FONT_SIZE_LABELS;
cb.FontSize = FONT_SIZE_AXES;

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);

xlim([0, min(max(r)*1e6, R_illu*3)]);
ylim([0, 40]);
box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig7_Forces_on_Temperature');
% export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');


%% force and gradt

figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Plot gradient magnitude with proper color scaling
% pcolor(R_grid*1e6, Z_grid*1e6, grad_T_mag_um);
contourf(R_grid*1e6, Z_grid*1e6, grad_T_mag_um, 150, 'LineStyle', 'none');

shading interp;
hold on;

% Overlay velocity field
skip_force = 15;
% quiver(R_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        Z_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        u_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        w_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        2, 'w', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_total_r_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_total_z_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       1, 'w ', 'LineWidth', 1, 'MaxHeadSize', 1);
xline(R_illu, 'r--', 'LineWidth', 2);
% plot([0, max(r_grid)*1e6], [30, 30], 'g--', 'LineWidth', 2.5);
% quiver(R_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        Z_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        u_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        w_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        1, 'w', 'LineWidth', 1.5,'ShowArrowHead','on');


% cb = colorbar;
% cb.Label.String = '$|\nabla T|$ (K/$\mu$m)';
% cb.Label.Interpreter = 'latex';
% cb.Label.FontSize = FONT_SIZE_LABELS;
% cb.FontSize = FONT_SIZE_AXES;
% cb.TickLabelInterpreter='latex'

% Set proper colormap and limits
colormap(thermal);
% colormap(flipud(spectral))
% Fix color scale - gradient is small!
caxis([0, max(grad_T_mag_um(:), [], 'omitnan')]);

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);

% xlim([0, min(max(r_grid)*1e6, R_illu*3)]);
% ylim([0, 50]);
ylim([0 120])

box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);


set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig2_FGradient_Flow');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
% fprintf('Saved: Fig2\n');

%%
fig = figure('Position', [100, 100, 1200, 500], 'Color', 'w');

% Panel A: Temperature
subplot(1,2,1);
contourf(R_grid*1e6, Z_grid*1e6, T_grid, 50, 'LineStyle', 'none');
colormap(thermal);
hold on;

skip_vec = 8;
quiver(R_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       Z_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       u_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       w_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       2, 'w', 'LineWidth', 1.5);

xline(R_illu, 'r--', 'LineWidth', 2);
plot([0, max(r)*1e6], [30, 30], 'g--', 'LineWidth', 2.5);

cb = colorbar;
cb.Label.String = 'Temperature (K)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = FONT_SIZE_LABELS;
cb.FontSize = FONT_SIZE_AXES;
cb.TickLabelInterpreter='latex';
xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);

% xlim([0, min(max(r)*1e6, R_illu*3)]);
% ylim([0, 50]);
box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig1_Temperature_Flow');
% Panel B: Gradient + Forces
subplot(1,2,2);
pcolor(R_grid*1e6, Z_grid*1e6, grad_T_mag_um);
shading interp;
hold on;

% Overlay velocity field
skip_vec = 8;
% quiver(R_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        Z_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        u_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        w_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
%        2, 'w', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_total_r_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_total_z_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       2.5, 'k ', 'LineWidth', 1, 'MaxHeadSize', 0.8);
xline(R_illu, 'r--', 'LineWidth', 2);
% plot([0, max(r_grid)*1e6], [30, 30], 'g--', 'LineWidth', 2.5);

cb = colorbar;
cb.Label.String = '$|\nabla T|$ (K/$\mu$m)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = FONT_SIZE_LABELS;
cb.FontSize = FONT_SIZE_AXES;
cb.TickLabelInterpreter='latex'

% Set proper colormap and limits
colormap(thermal);
% Fix color scale - gradient is small!
caxis([0, max(grad_T_mag_um(:), [], 'omitnan')]);

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);

% xlim([0, min(max(r_grid)*1e6, R_illu*3)]);
% ylim([0, 50]);
box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig2_Gradient_Flow');

set(gcf, 'Units', 'inches', 'Position', [1, 1, 12, 5]);
% export_fig(fullfile(output_dir, 'Fig6_Overview_2Panel'), ...
%     '-pdf', '-painters', '-r300', '-transparent');











%%
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% --- Compute dot product F · u (before scaling used in quiver) ---
FdotU = F_total_r_grid .* u_grid + F_total_z_grid .* w_grid;

% --- Normalize the dot product to [-1, 1] ---
FdotU_norm = FdotU ./ max(abs(FdotU(:)), [], 'omitnan');

% Plot DOT PRODUCT (normalized) instead of gradient magnitude
pcolor(R_grid*1e6, Z_grid*1e6, FdotU_norm);
shading interp;
hold on;

% Overlay velocity field
skip_vec = 2;

% Force field (unchanged)
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_total_r_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_total_z_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       2.5, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.8);

xline(R_illu, 'r--', 'LineWidth', 2);

% Flow field (unchanged)
quiver(R_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       Z_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       u_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       w_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       1, 'w', 'LineWidth', 1.5,'ShowArrowHead','on');

% Colormap (diverging)
colormap(balanced);

% Now the range is exactly [-1, 1]
caxis([-1 1]);

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);

box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

%%
%% ========================================================================
%% SURFACE TEMPERATURE GRADIENT AT MULTIPLE HEIGHTS
%% Compute ∂T/∂r at different z-heights (not just z=0)
%% ========================================================================

fprintf('\n=== Computing temperature gradients at multiple heights ===\n');

% Define heights to analyze
z_heights_um = [0, 5, 10, 15, 20, 25, 30];  % μm
n_heights = length(z_heights_um);
z_tolerance_um = 2;  % ±2 μm tolerance

% Binning parameters for radial profiles
n_bins_grad = 60;  % Good resolution for gradient calculation
r_max = max(r);
r_edges_grad = linspace(0, r_max, n_bins_grad+1);
r_centers_grad = (r_edges_grad(1:end-1) + r_edges_grad(2:end))/2 * 1e6;  % in μm

% Storage for results
gradT_results = struct();
gradT_results.r_centers = r_centers_grad;
gradT_results.z_heights = z_heights_um;

% Colors for plotting
colors_grad = jet(n_heights);

%% Extract and compute gradient for each height
for i_h = 1:n_heights
    z_target = z_heights_um(i_h) * 1e-6;
    z_tol = z_tolerance_um * 1e-6;
    
    % Extract data at this height
    idx_height = abs(z - z_target) < z_tol;
    
    if sum(idx_height) < 10
        fprintf('  z=%.0fμm: Insufficient points, skipping\n', z_heights_um(i_h));
        continue;
    end
    
    r_h = r(idx_height);
    T_h = T(idx_height);
    z_actual = mean(z(idx_height)) * 1e6;
    
    fprintf('  z=%.0fμm: %d points (actual: %.1fμm)\n', ...
            z_heights_um(i_h), sum(idx_height), z_actual);
    
    % Bin temperature by radius
    T_mean = nan(n_bins_grad, 1);
    T_std = nan(n_bins_grad, 1);
    
    for i_bin = 1:n_bins_grad
        idx_bin = r_h >= r_edges_grad(i_bin) & r_h < r_edges_grad(i_bin+1);
        if sum(idx_bin) >= 2
            T_mean(i_bin) = mean(T_h(idx_bin));
            T_std(i_bin) = std(T_h(idx_bin));
        end
    end
    
    % Smooth temperature profile slightly
    T_mean_smooth = movmean(T_mean, 3, 'omitnan');
    
    % Compute gradient ∂T/∂r
    valid = ~isnan(T_mean_smooth);
    if sum(valid) > 5
        % Use central differences for gradient
        gradT_r = gradient(T_mean_smooth(valid), r_centers_grad(valid));
        
        % Convert to K/μm
        gradT_r_um = gradT_r * 1e-6;
        
        % Store results
        field_name = sprintf('z%03d', round(z_heights_um(i_h)*10));
        gradT_results.(field_name).r = r_centers_grad(valid);
        gradT_results.(field_name).T_mean = T_mean_smooth(valid);
        gradT_results.(field_name).gradT_r = gradT_r_um;
        gradT_results.(field_name).z_actual = z_actual;
        
        % Find minimum (steepest negative gradient)
        [min_grad, idx_min] = min(gradT_r_um);
        r_min = r_centers_grad(valid);
        r_min = r_min(idx_min);
        
        gradT_results.(field_name).min_grad = min_grad;
        gradT_results.(field_name).r_min = r_min;
        
        fprintf('    → Min gradient: %.4f K/μm at r=%.1f μm\n', min_grad, r_min);
    end
end

%% FIGURE: Temperature Gradient ∂T/∂r at Multiple Heights
figure('Position', [100, 100, 900, 700], 'Color', 'w');
hold on;

for i_h = 1:n_heights
    field_name = sprintf('z%03d', round(z_heights_um(i_h)*10));
    
    if ~isfield(gradT_results, field_name)
        continue;
    end
    
    data_h = gradT_results.(field_name);
    
    plot(data_h.r, data_h.gradT_r, '-', ...
         'Color', colors_grad(i_h,:), 'LineWidth', LINE_WIDTH, ...
         'DisplayName', sprintf('$z = %.0f$ $\\mu$m', data_h.z_actual));
    % 
    % % Mark minimum
    % plot(data_h.r_min, data_h.min_grad, 'o', ...
    %      'Color', colors_grad(i_h,:), 'MarkerSize', 8, ...
    %      'MarkerFaceColor', colors_grad(i_h,:), 'HandleVisibility', 'off');
end

% Reference lines
yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(R_illu, 'r--', 'LineWidth', 2, 'Alpha', 0.5, ...
      'DisplayName', 'Illumination edge');

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$\partial T/\partial r$ (K/$\mu$m)', 'FontSize', FONT_SIZE_LABELS);

% legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-2, 'NumColumns', 2);
% grid on;
box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);
xlim([0, max(r_centers_grad)]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH/2, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig_GradT_MultiHeight');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');

%%
