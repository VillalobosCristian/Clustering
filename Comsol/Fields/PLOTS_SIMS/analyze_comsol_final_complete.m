clear all; clc; close all;

%% FIGURE PARAMETERS
FIGURE_WIDTH = 10;
FIGURE_HEIGHT = 5;
FONT_SIZE_LABELS = 24;
FONT_SIZE_AXES = 24;
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

z_layers_um = [0, 5, 10, 15, 20, 30];  % 6 well-spaced layers
z_tolerance_um = 2;  % ±2 μm
n_layers = length(z_layers_um);
colors_layers = jet(n_layers);

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

contourf(R_grid*1e6, Z_grid*1e6, T_grid, 50, 'LineStyle', 'none');
colormap(thermal);
hold on;

skip_vec = 15;
quiver(R_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       Z_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       u_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       w_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       1, 'k', 'LineWidth', 1,'MaxHeadSize',1.5);

xline(R_illu, 'r--', 'LineWidth', 2);
% plot([0, max(r)*1e6], [30, 30], 'g--', 'LineWidth', 2.5);

% cb = colorbar;
% cb.Label.String = 'T (K)';
% cb.Label.Interpreter = 'latex';
% cb.Label.FontSize = FONT_SIZE_LABELS;
% cb.FontSize = FONT_SIZE_AXES;
% cb.TickLabelInterpreter='latex';
xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
% 
% xlim([0, min(max(r)*1e6, R_illu*3)]);
% ylim([0, 50]);
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
pcolor(R_grid*1e6, Z_grid*1e6, grad_T_mag_um);
shading interp;
hold on;

% Overlay velocity field
skip_vec = 15;
quiver(R_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       Z_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       u_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       w_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       1, 'k', 'LineWidth', 1,'MaxHeadSize',1.5);


xline(R_illu, 'r--', 'LineWidth', 2);
% plot([0, max(r_grid)*1e6], [30, 30], 'g--', 'LineWidth', 2.5);

% cb = colorbar;
% cb.Label.String = '$|\nabla T|$ (K/$\mu$m)';
% cb.Label.Interpreter = 'latex';
% cb.Label.FontSize = FONT_SIZE_LABELS;
% cb.FontSize = FONT_SIZE_AXES;
% cb.TickLabelInterpreter='latex'

% Set proper colormap and limits
colormap(thermal);
% Fix color scale - gradient is small!
caxis([0, max(grad_T_mag_um(:), [], 'omitnan')]);

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
% 
% xlim([0, min(max(r_grid)*1e6, R_illu*3)]);
% ylim([0, 50]);
box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig2_Gradient_Flow');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Fig2\n');

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
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Fig3\n');

%% ========================================================================
%% FIGURE 4: Velocity Profiles at Multiple Heights
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');
hold on;

for i_layer = 1:n_layers
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
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-4, 'NumColumns', 2);
grid on; box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig4_Velocity_vs_Height');
% export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Fig4\n');

%% ========================================================================
%% FIGURE 5: Force Balance Components (NEW!)
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
yline(4.3, 'k:', 'LineWidth', 1.5, 'DisplayName', '$k_BT/a$');
xline(R_illu, 'r--', 'LineWidth', 1.5, 'Alpha', 0.5, 'HandleVisibility', 'off');

% Add annotations explaining force directions
y_range = ylim;
text(0.05, 0.95, '$F_{TP} > 0$ (outward)', 'Units', 'normalized', ...
     'FontSize', FONT_SIZE_LEGEND-2, 'Color', 'r', 'VerticalAlignment', 'top');
text(0.05, 0.88, '$F_D < 0$ (inward)', 'Units', 'normalized', ...
     'FontSize', FONT_SIZE_LEGEND-2, 'Color', 'b', 'VerticalAlignment', 'top');
text(0.05, 0.81, '$F_{total} < 0$ (net inward)', 'Units', 'normalized', ...
     'FontSize', FONT_SIZE_LEGEND-2, 'Color', 'k', 'VerticalAlignment', 'top');

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('Force (fN)', 'FontSize', FONT_SIZE_LABELS);
title(sprintf('Force Balance at $z = %.0f$ $\\mu$m', profiles.z150.z_actual), ...
      'FontSize', FONT_SIZE_LABELS-2);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-2);
grid on; box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig5_Force_Balance');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Fig5\n');

%% ========================================================================
%% FIGURE 6: Force Components at Multiple Heights
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
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Fig6\n');

%% ========================================================================
%% FIGURE 7: Forces on Temperature Field (2D)
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');

contourf(R_grid*1e6, Z_grid*1e6, T_grid, 50, 'LineStyle', 'none');
colormap(balanced);
hold on;

% Overlay force vectors (sparse)
skip_force = 15;
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_total_r_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_total_z_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       1, 'w', 'LineWidth', 1, 'MaxHeadSize', 1);

xline(R_illu, 'r--', 'LineWidth', 2);
% plot([0, max(r)*1e6], [30, 30], 'g--', 'LineWidth', 2.5);

cb = colorbar;
cb.Label.String = 'Temperature (K)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = FONT_SIZE_LABELS;
cb.FontSize = FONT_SIZE_AXES;
cb.TickLabelInterpreter = 'latex';

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);

% xlim([0, min(max(r)*1e6, R_illu*3)]);
% ylim([0, 50]);
box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig7_Forces_on_Temperature');
% export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
% fprintf('Saved: Fig7\n');

%% FINAL SUMMARY
fprintf('\n========================================\n');
fprintf('Analysis Complete\n');
fprintf('========================================\n');
fprintf('Temperature: ΔT = %.2f K\n', Delta_T);
fprintf('Max |∇T| = %.4f K/μm\n', max(grad_T_mag_um(:), [], 'omitnan'));
fprintf('R_illum = %.0f μm\n', R_illu);
fprintf('Layers analyzed: %d (z = 0-30 μm)\n', n_layers);
fprintf('Radial bins: %d\n', n_bins);
fprintf('\n7 figures saved to:\n%s\n', output_dir);
fprintf('========================================\n');

%% Print Force Balance Summary
fprintf('\n=== FORCE BALANCE SUMMARY (z=15 μm) ===\n');
if isfield(profiles, 'z150')
    prof = profiles.z150;
    valid = ~isnan(prof.F_TP_mean) & ~isnan(prof.F_D_mean);
    
    if sum(valid) > 0
        fprintf('Thermophoresis (F_TP): %.2f to %.2f fN\n', ...
            min(prof.F_TP_mean(valid)), max(prof.F_TP_mean(valid)));
        fprintf('Drag/Convection (F_D): %.2f to %.2f fN\n', ...
            min(prof.F_D_mean(valid)), max(prof.F_D_mean(valid)));
        fprintf('Total (F_total): %.2f to %.2f fN\n', ...
            min(prof.F_mean(valid)), max(prof.F_mean(valid)));
        fprintf('Thermal force (kT/a): 4.3 fN\n');
    end
end
fprintf('========================================\n');


%%

figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Plot gradient magnitude with proper color scaling
pcolor(R_grid*1e6, Z_grid*1e6, grad_T_mag_um);
shading interp;
hold on;

% Overlay velocity field
skip_vec = 15;
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
% 
% cb = colorbar;
% cb.Label.String = '$|\nabla T|$ (K/$\mu$m)';
% cb.Label.Interpreter = 'latex';
% cb.Label.FontSize = FONT_SIZE_LABELS;
% cb.FontSize = FONT_SIZE_AXES;
% cb.TickLabelInterpreter='latex'

% Set proper colormap and limits
colormap(thermal);
% Fix color scale - gradient is small!
caxis([0, max(grad_T_mag_um(:), [], 'omitnan')]);

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);

% xlim([0, min(max(r_grid)*1e6, R_illu*3)]);
% ylim([0, 50]);
% xlim([0 150])
box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig2_Gradient_FFlow');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');

%%
% For each height:
%% ========================================================================
%% COMPUTE PARTICLE VELOCITY (including thermophoresis)
%% Add this after your existing velocity profile analysis
%% ========================================================================

%% This computes v_particle = u_fluid - D_T × ∂T/∂r
%% Currently you plot u_fluid, but particles experience v_particle!

fprintf('\n=== Computing Particle Velocity (with thermophoresis) ===\n');

% % Physical parameter
% D_T = 0.2e-12;  % m²/(s·K) - thermophoretic mobility for PS in water
% 
% % For each layer, compute particle velocity
% for i_layer = 1:n_layers
%     field_name = sprintf('z%03d', round(z_layers_um(i_layer)*10));
%     if ~isfield(profiles, field_name), continue; end
% 
%     z_val = z_layers_um(i_layer);
%     z_target = z_val * 1e-6;
%     z_tol = z_tolerance_um * 1e-6;
% 
%     % Extract temperature at this height for gradient calculation
%     idx_layer = abs(z - z_target) < z_tol;
% 
%     if sum(idx_layer) < 10, continue; end
% 
%     r_layer = r(idx_layer);
%     T_layer = T(idx_layer);
% 
%     % Bin temperature (same bins as velocity)
%     T_binned = nan(n_bins, 1);
%     for i_bin = 1:n_bins
%         idx_bin = r_layer >= r_edges(i_bin) & r_layer < r_edges(i_bin+1);
%         if sum(idx_bin) >= 2
%             T_binned(i_bin) = mean(T_layer(idx_bin));
%         end
%     end
% 
%     % Smooth temperature
%     T_binned = movmean(T_binned, 2, 'omitnan');
% 
%     % Compute radial gradient ∂T/∂r
%     valid = ~isnan(T_binned) & ~isnan(profiles.(field_name).u_mean);
%     if sum(valid) < 3, continue; end
% 
%     r_vals = profiles.r_centers(valid) * 1e-6;  % Convert to meters
%     T_vals = T_binned(valid);
% 
%     % Compute gradient using central differences
%     dT_dr = gradient(T_vals, r_vals);  % K/m
% 
%     % Thermophoretic velocity contribution: -D_T × ∂T/∂r
%     % Note: ∂T/∂r < 0 near R_illum, so -D_T × ∂T/∂r > 0 (outward)
%     v_thermophoresis = -D_T * dT_dr * 1e6;  % Convert to μm/s
% 
%     % Particle velocity = fluid velocity + thermophoretic contribution
%     u_fluid = profiles.(field_name).u_mean(valid);
%     v_particle = u_fluid + v_thermophoresis;
% 
%     % Store in profiles structure
%     profiles.(field_name).v_particle = nan(n_bins, 1);
%     profiles.(field_name).v_particle(valid) = v_particle;
%     profiles.(field_name).v_thermophoresis = nan(n_bins, 1);
%     profiles.(field_name).v_thermophoresis(valid) = v_thermophoresis;
%     profiles.(field_name).dT_dr = nan(n_bins, 1);
%     profiles.(field_name).dT_dr(valid) = dT_dr * 1e-6;  % K/μm for display
% 
%     % Print statistics
%     fprintf('z=%.0fμm:\n', profiles.(field_name).z_actual);
%     fprintf('  Fluid velocity (u_fluid): %.3f to %.3f μm/s\n', ...
%         min(u_fluid), max(u_fluid));
%     fprintf('  Thermophoresis (v_thermo): %.3f to %.3f μm/s (outward)\n', ...
%         min(v_thermophoresis), max(v_thermophoresis));
%     fprintf('  Particle velocity (v_particle): %.3f to %.3f μm/s\n', ...
%         min(v_particle), max(v_particle));
%     fprintf('  Difference: %.1f%% reduction in inward velocity\n', ...
%         (1 - min(v_particle)/min(u_fluid))*100);
% end
% 
% fprintf('\n');

%% ========================================================================
%% FIGURE: Fluid vs Particle Velocity Comparison
%% ========================================================================
% 
% figure('Position', [100, 100, 900, 700], 'Color', 'w');
% hold on;
% 
% % Focus on z=15 μm (particle observation height)
% field_name = 'z150';
% if isfield(profiles, field_name)
%     prof = profiles.(field_name);
%     valid = ~isnan(prof.u_mean) & ~isnan(prof.v_particle);
% 
%     if sum(valid) > 0
%         % Plot fluid velocity
%         plot(profiles.r_centers(valid), prof.u_mean(valid), 'b-', ...
%              'LineWidth', LINE_WIDTH, 'DisplayName', ...
%              '$u_{\\mathrm{fluid}}$ (what you currently plot)');
% 
%         % Plot particle velocity
%         plot(profiles.r_centers(valid), prof.v_particle(valid), 'r-', ...
%              'LineWidth', LINE_WIDTH, 'DisplayName', ...
%              '$v_{\\mathrm{particle}}$ (actual particle motion)');
% 
%         % Plot thermophoretic contribution
%         plot(profiles.r_centers(valid), prof.v_thermophoresis(valid), 'g--', ...
%              'LineWidth', LINE_WIDTH-0.5, 'DisplayName', ...
%              '$-D_T \\partial T/\\partial r$ (thermophoresis)');
% 
%         % Shade region showing the difference
%         x_fill = profiles.r_centers(valid);
%         y1 = prof.u_mean(valid);
%         y2 = prof.v_particle(valid);
%         fill([x_fill; flipud(x_fill)], [y1; flipud(y2)], 'y', ...
%              'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
%              'DisplayName', 'Thermophoretic reduction');
%     end
% end
% 
% yline(0, 'k--', 'LineWidth', 1.5);
% xline(R_illu, 'r--', 'LineWidth', 1.5, 'Alpha', 0.5, 'HandleVisibility', 'off');
% 
% xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
% ylabel('Velocity ($\mu$m/s)', 'FontSize', FONT_SIZE_LABELS);
% title('Fluid vs Particle Velocity at $z = 15$ $\mu$m', 'FontSize', FONT_SIZE_LABELS-2);
% legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-2);
% 
% % Add annotation
% text(0.05, 0.95, {'$v_{\\mathrm{particle}} = u_{\\mathrm{fluid}} - D_T \\frac{\\partial T}{\\partial r}$', ...
%                    '', ...
%                    'Thermophoresis opposes inward flow', ...
%                    'but convection dominates'}, ...
%      'Units', 'normalized', 'FontSize', FONT_SIZE_LEGEND-2, ...
%      'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', 'k');
% 
% grid on; box on;
% set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);
% 
% set(gcf, 'Units', 'inches');
% set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);
% 
% fileName = fullfile(output_dir, 'Fig_Fluid_vs_Particle_Velocity');
% export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
% fprintf('Saved: Fluid vs Particle Velocity comparison\n');

%% ========================================================================
%% FIGURE: All Heights - Particle Velocity
% %% ========================================================================
% 
% figure('Position', [100, 100, 900, 700], 'Color', 'w');
% hold on;
% 
% for i_layer = 1:n_layers
%     field_name = sprintf('z%03d', round(z_layers_um(i_layer)*10));
%     if ~isfield(profiles, field_name), continue; end
%     if ~isfield(profiles.(field_name), 'v_particle'), continue; end
% 
%     prof = profiles.(field_name);
%     valid = ~isnan(prof.v_particle);
% 
%     if sum(valid) > 0
%         plot(profiles.r_centers(valid), prof.v_particle(valid), '-', ...
%              'Color', colors_layers(i_layer,:), 'LineWidth', LINE_WIDTH, ...
%              'DisplayName', sprintf('$z = %.0f$ $\\mu$m', prof.z_actual));
%     end
% end
% 
% yline(0, 'k--', 'LineWidth', 1.5);
% xline(R_illu, 'r--', 'LineWidth', 1.5, 'Alpha', 0.5);
% 
% xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
% ylabel('$v_{\\mathrm{particle}}$ ($\mu$m/s)', 'FontSize', FONT_SIZE_LABELS);
% title('Particle Velocity (with thermophoresis) vs Height', 'FontSize', FONT_SIZE_LABELS-2);
% legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-4, 'NumColumns', 2);
% grid on; box on;
% set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);
% 
% set(gcf, 'Units', 'inches');
% set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);
% 
% fileName = fullfile(output_dir, 'Fig_Particle_Velocity_vs_Height');
% export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
% fprintf('Saved: Particle velocity at multiple heights\n');
% 
% %% ========================================================================
% %% FIGURE: Temperature Gradient (shows why thermophoresis is outward)
% %% ========================================================================
% 
% figure('Position', [100, 100, 900, 700], 'Color', 'w');
% 
% field_name = 'z150';
% if isfield(profiles, field_name) && isfield(profiles.(field_name), 'dT_dr')
%     prof = profiles.(field_name);
%     valid = ~isnan(prof.dT_dr);
% 
%     if sum(valid) > 0
%         plot(profiles.r_centers(valid), prof.dT_dr(valid), 'b-', ...
%              'LineWidth', LINE_WIDTH);
%         hold on;
% 
%         yline(0, 'k--', 'LineWidth', 1.5);
%         xline(R_illu, 'r--', 'LineWidth', 2);
% 
%         % Mark region where ∂T/∂r < 0 (thermophoresis outward)
%         idx_neg = prof.dT_dr(valid) < 0;
%         r_neg = profiles.r_centers(valid);
%         if sum(idx_neg) > 0
%             y_range = ylim;
%             fill([r_neg(idx_neg); flipud(r_neg(idx_neg))], ...
%                  [ones(sum(idx_neg),1)*y_range(1); ones(sum(idx_neg),1)*0], ...
%                  'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
%             text(mean(r_neg(idx_neg)), y_range(1)*0.8, ...
%                  '$\\partial T/\\partial r < 0$ $\\rightarrow$ thermophoresis outward', ...
%                  'FontSize', FONT_SIZE_LEGEND, 'Color', 'r', ...
%                  'HorizontalAlignment', 'center');
%         end
% 
%         xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
%         ylabel('$\\partial T/\\partial r$ (K/$\mu$m)', 'FontSize', FONT_SIZE_LABELS);
%         title(sprintf('Temperature Gradient at $z = %.0f$ $\\mu$m', prof.z_actual), ...
%               'FontSize', FONT_SIZE_LABELS-2);
% 
%         grid on; box on;
%         set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);
%     end
% end
% 
% set(gcf, 'Units', 'inches');
% set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);
% 
% fileName = fullfile(output_dir, 'Fig_Temperature_Gradient_z15');
% export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
% fprintf('Saved: Temperature gradient at z=15μm\n');
% 
% %% Print Summary
% fprintf('\n========================================\n');
% fprintf('PARTICLE VS FLUID VELOCITY SUMMARY\n');
% fprintf('========================================\n');
% 
% if isfield(profiles, 'z150') && isfield(profiles.z150, 'v_particle')
%     prof = profiles.z150;
%     valid = ~isnan(prof.v_particle) & ~isnan(prof.u_mean);
% 
%     if sum(valid) > 0
%         u_min = min(prof.u_mean(valid));
%         v_min = min(prof.v_particle(valid));
%         thermo_contrib = mean(prof.v_thermophoresis(valid));
% 
%         fprintf('At z = %.0f μm:\n', prof.z_actual);
%         fprintf('  Min fluid velocity:    %.3f μm/s (inward)\n', u_min);
%         fprintf('  Thermophoresis:        %.3f μm/s (outward)\n', thermo_contrib);
%         fprintf('  Min particle velocity: %.3f μm/s (inward)\n', v_min);
%         fprintf('\n');
%         fprintf('Thermophoresis reduces inward velocity by %.1f%%\n', ...
%             (1 - v_min/u_min)*100);
%         fprintf('But particle still accumulates inward! ✓\n');
%         fprintf('Convection/Drag dominates thermophoresis\n');
%     end
% end
% 
% fprintf('========================================\n\n');