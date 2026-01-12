clear all; clc; close all;

%% ========================================================================
%% COMSOL SIMULATION ANALYSIS
%% ========================================================================
% Analysis of optothermal flow simulations
% Matches style and metrics from experimental analysis (CSV_analysis_colloids.m)
%
% Input: Fields.csv or Fields_superGaussian.csv from COMSOL
% Output: Flow fields, velocity profiles, gradients, forces
%% ========================================================================

set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');

%% ========================================================================
%% SECTION 1: LOAD DATA
%% ========================================================================

fprintf('\n========================================================================\n');
fprintf('COMSOL SIMULATION ANALYSIS\n');
fprintf('========================================================================\n\n');

% Select file to analyze
[file_name, path_name] = uigetfile('*.csv', 'Select COMSOL Fields CSV file');
if isequal(file_name, 0)
    error('No file selected');
end

% Extract simulation name
[~, sim_name, ~] = fileparts(file_name);
fprintf('Analyzing: %s\n\n', sim_name);

% Ask for simulation parameters
prompt = {'Profile type (Gaussian/SuperGaussian):', ...
          'Characteristic width w (μm):', ...
          'Shape exponent d:'};
dlgtitle = 'Simulation Parameters';
dims = [1 50];
defaultans = {'SuperGaussian', '22.4', '5.5'};
answer = inputdlg(prompt, dlgtitle, dims, defaultans);

if isempty(answer)
    error('Parameters not provided');
end

profile_type = answer{1};
w_um = str2double(answer{2});
d_shape = str2double(answer{3});

fprintf('Profile: %s\n', profile_type);
fprintf('Width w: %.1f μm\n', w_um);
fprintf('Shape exponent d: %.1f\n\n', d_shape);

%% Load CSV data
fprintf('Loading data from: %s\n', file_name);
data = readtable(fullfile(path_name, file_name), 'HeaderLines', 8);

% Assign column names
data.Properties.VariableNames = {'r', 'z', 'T', 'F_TPr', 'F_Dr', ...
                                  'F_TPz', 'F_Dz', 'F_g', 'u', 'w'};

% Convert to micrometers and micrometers/second
data.r_um = data.r * 1e6;
data.z_um = data.z * 1e6;
data.u_um = data.u * 1e6;  % radial velocity
data.w_um = data.w * 1e6;  % vertical velocity

fprintf('Loaded %d data points\n', height(data));
fprintf('Domain: r = [%.1f, %.1f] μm, z = [%.1f, %.1f] μm\n\n', ...
        min(data.r_um), max(data.r_um), min(data.z_um), max(data.z_um));

%% Create output directories
results_base = fullfile(pwd, 'COMSOL_results', sim_name);
fig_dir = fullfile(results_base, 'figures');
data_dir = fullfile(results_base, 'data');

if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

fprintf('Output folder: %s\n\n', results_base);

%% ========================================================================
%% SECTION 2: TEMPERATURE FIELD ANALYSIS
%% ========================================================================

fprintf('=== TEMPERATURE FIELD ===\n');

T_ref = 293.15;  % K (20°C)
Delta_T = max(data.T) - T_ref;

fprintf('T_min = %.2f K (%.1f °C)\n', min(data.T), min(data.T)-273.15);
fprintf('T_max = %.2f K (%.1f °C)\n', max(data.T), max(data.T)-273.15);
fprintf('ΔT = %.2f K\n\n', Delta_T);

%% FIGURE 1: Temperature field (2D map)
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Subsample for plotting if too many points
if height(data) > 50000
    idx_plot = 1:5:height(data);
else
    idx_plot = 1:height(data);
end

scatter(data.r_um(idx_plot), data.z_um(idx_plot), 8, data.T(idx_plot), 'filled');
colormap(hot);
cb = colorbar;
cb.Label.String = '$T$ (K)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$z\ (\mu\mathrm{m})$', 'FontSize', 18);
title('Temperature field $T(r,z)$', 'FontSize', 18);

axis equal tight;
grid on;
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'temperature_field';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% FIGURE 2: Temperature profile at surface (z < 5 μm)
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Extract near-surface data
z_surface_max = 5;  % μm
surf_data = data(data.z_um < z_surface_max, :);
surf_data = sortrows(surf_data, 'r_um');

% Bin and average
n_bins = 50;
r_edges = linspace(0, max(surf_data.r_um), n_bins+1);
r_centers = (r_edges(1:end-1) + r_edges(2:end))/2;

T_mean = zeros(n_bins, 1);
T_std = zeros(n_bins, 1);

for k = 1:n_bins
    idx = surf_data.r_um >= r_edges(k) & surf_data.r_um < r_edges(k+1);
    if sum(idx) > 0
        T_mean(k) = mean(surf_data.T(idx));
        T_std(k) = std(surf_data.T(idx));
    else
        T_mean(k) = NaN;
        T_std(k) = NaN;
    end
end

% Plot
plot(r_centers, T_mean - T_ref, 'b-', 'LineWidth', 3);
hold on;
% fill([r_centers, fliplr(r_centers)], ...
%      [T_mean - T_ref + T_std, fliplr(T_mean - T_ref - T_std)], ...
%      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$\Delta T$ (K)', 'FontSize', 18);
title(sprintf('Temperature profile at surface ($z < %.0f$ $\\mu$m)', z_surface_max), 'FontSize', 18);

grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;
xlim([0, max(r_centers)]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'temperature_profile_surface';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% ========================================================================
%% SECTION 3: TEMPERATURE GRADIENTS
%% ========================================================================

fprintf('\n=== TEMPERATURE GRADIENTS ===\n');

% Compute gradients using gridded data
% Create interpolant for smooth gradient calculation
r_grid = linspace(0, max(data.r_um), 100);
z_grid = linspace(0, max(data.z_um), 100);
[R_grid, Z_grid] = meshgrid(r_grid, z_grid);

% Interpolate temperature field
F_T = scatteredInterpolant(data.r_um, data.z_um, data.T, 'linear', 'none');
T_grid = F_T(R_grid, Z_grid);

% Calculate gradients
[dT_dr, dT_dz] = gradient(T_grid, r_grid(2)-r_grid(1), z_grid(2)-z_grid(1));
grad_T_mag = sqrt(dT_dr.^2 + dT_dz.^2);

fprintf('Max |∇T|: %.3f K/μm\n', max(grad_T_mag(:)));

%% FIGURE 3: Temperature gradient magnitude
figure('Position', [100, 100, 900, 700], 'Color', 'w');

contourf(R_grid, Z_grid, grad_T_mag, 20, 'LineStyle', 'none');
colormap(jet);
cb = colorbar;
cb.Label.String = '$|\nabla T|$ (K/$\mu$m)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$z\ (\mu\mathrm{m})$', 'FontSize', 18);
title('Temperature gradient magnitude $|\nabla T|$', 'FontSize', 18);

axis equal tight;
grid on;
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'gradient_magnitude';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% FIGURE 4: Radial gradient at surface
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Extract gradient at surface
grad_T_surf = zeros(size(r_centers));
for k = 1:length(r_centers)
    % Find closest r in grid
    [~, ir] = min(abs(r_grid - r_centers(k)));
    % Average over near-surface z points
    iz_surf = find(z_grid < z_surface_max);
    if ~isempty(iz_surf)
        grad_T_surf(k) = mean(grad_T_mag(iz_surf, ir), 'omitnan');
    end
end

plot(r_centers, grad_T_surf, 'r-', 'LineWidth', 3);

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$|\nabla T|$ (K/$\mu$m)', 'FontSize', 18);
title(sprintf('Temperature gradient at surface ($z < %.0f$ $\\mu$m)', z_surface_max), 'FontSize', 18);

grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;
xlim([0, max(r_centers)]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'gradient_radial_surface';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% ========================================================================
%% SECTION 4: VELOCITY FIELD ANALYSIS
%% ========================================================================

fprintf('\n=== VELOCITY FIELD ===\n');

% Velocity magnitude
data.v_mag = sqrt(data.u_um.^2 + data.w_um.^2);

fprintf('Radial velocity u_r:\n');
fprintf('  Min: %.3f μm/s (max inflow)\n', min(data.u_um));
fprintf('  Max: %.3f μm/s\n', max(data.u_um));

fprintf('Vertical velocity u_z:\n');
fprintf('  Min: %.3f μm/s\n', min(data.w_um));
fprintf('  Max: %.3f μm/s (max upwelling)\n', max(data.w_um));

fprintf('Velocity magnitude:\n');
fprintf('  Max: %.3f μm/s\n\n', max(data.v_mag));

%% FIGURE 5: Velocity magnitude field
figure('Position', [100, 100, 900, 700], 'Color', 'w');

scatter(data.r_um(idx_plot), data.z_um(idx_plot), 8, data.v_mag(idx_plot), 'filled');
colormap(parula);
cb = colorbar;
cb.Label.String = '$|\mathbf{v}|$ ($\mu$m/s)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$z\ (\mu\mathrm{m})$', 'FontSize', 18);
title('Velocity magnitude $|\mathbf{v}|$', 'FontSize', 18);

axis equal tight;
grid on;
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'velocity_magnitude';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% FIGURE 6: Radial velocity field (color-coded)
figure('Position', [100, 100, 900, 700], 'Color', 'w');

scatter(data.r_um(idx_plot), data.z_um(idx_plot), 8, data.u_um(idx_plot), 'filled');
colormap(redblue_colormap);
caxis([-max(abs(data.u_um)), max(abs(data.u_um))]);
cb = colorbar;
cb.Label.String = '$v_r$ ($\mu$m/s)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$z\ (\mu\mathrm{m})$', 'FontSize', 18);
title('Radial velocity $v_r$ (negative = inflow)', 'FontSize', 18);

axis equal tight;
grid on;
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'velocity_radial_field';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% ========================================================================
%% SECTION 5: SPATIALLY-AVERAGED RADIAL VELOCITY <v_r>(r)
%% ========================================================================

fprintf('\n=== RADIAL VELOCITY PROFILE ===\n');

% Bin by radius and average (similar to experimental analysis)
vr_mean = zeros(n_bins, 1);
vr_sem = zeros(n_bins, 1);

for k = 1:n_bins
    idx = surf_data.r_um >= r_edges(k) & surf_data.r_um < r_edges(k+1);
    if sum(idx) > 5
        vr_mean(k) = mean(surf_data.u_um(idx));
        vr_sem(k) = std(surf_data.u_um(idx)) / sqrt(sum(idx));
    else
        vr_mean(k) = NaN;
        vr_sem(k) = NaN;
    end
end

% Find minimum (maximum inflow)
[vr_min, idx_min] = min(vr_mean);
r_min = r_centers(idx_min);

fprintf('Maximum inflow: %.3f μm/s at r = %.1f μm\n', vr_min, r_min);

%% FIGURE 7: Spatially-averaged radial velocity
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Plot individual cells
scatter(surf_data.r_um, surf_data.u_um, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.2);

hold on;

% Plot binned average with error bars
errorbar(r_centers, vr_mean, vr_sem, 'ro-', 'LineWidth', 2.5, ...
         'MarkerSize', 10, 'MarkerFaceColor', 'r', 'CapSize', 10);

yline(0, 'k--', 'LineWidth', 1.5);

% Mark minimum location
plot(r_min, vr_min, 'gs', 'MarkerSize', 15, 'LineWidth', 3);
text(r_min, vr_min, sprintf('  Min: %.2f $\\mu$m/s\n  at $r=%.1f$ $\\mu$m', vr_min, r_min), ...
     'FontSize', 12, 'VerticalAlignment', 'top');

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$\langle v_r \rangle\ (\mu\mathrm{m/s})$', 'FontSize', 18);
title('Spatially-averaged radial velocity', 'FontSize', 18);

legend('Grid cells', 'Binned average', 'Zero line', 'Maximum inflow', ...
       'Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);

grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;
xlim([0, max(r_centers)]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'vr_radial_averaged';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% FIGURE 8: Vertical velocity profile
figure('Position', [100, 100, 900, 700], 'Color', 'w');

wz_mean = zeros(n_bins, 1);
wz_sem = zeros(n_bins, 1);

for k = 1:n_bins
    idx = surf_data.r_um >= r_edges(k) & surf_data.r_um < r_edges(k+1);
    if sum(idx) > 5
        wz_mean(k) = mean(surf_data.w_um(idx));
        wz_sem(k) = std(surf_data.w_um(idx)) / sqrt(sum(idx));
    else
        wz_mean(k) = NaN;
        wz_sem(k) = NaN;
    end
end

scatter(surf_data.r_um, surf_data.w_um, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.2);
hold on;
errorbar(r_centers, wz_mean, wz_sem, 'ro-', 'LineWidth', 2.5, ...
         'MarkerSize', 10, 'MarkerFaceColor', 'r', 'CapSize', 10);
yline(0, 'k--', 'LineWidth', 1.5);

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$\langle v_z \rangle\ (\mu\mathrm{m/s})$', 'FontSize', 18);
title('Spatially-averaged vertical velocity', 'FontSize', 18);

legend('Grid cells', 'Binned average', 'Zero line', ...
       'Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);

grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;
xlim([0, max(r_centers)]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'vz_radial_averaged';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% ========================================================================
%% SECTION 6: VELOCITY VECTOR FIELD & STREAMLINES
%% ========================================================================

fprintf('\n=== VELOCITY VECTORS AND STREAMLINES ===\n');

%% FIGURE 9: Vector field (quiver)
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Subsample for vectors
step = max(1, floor(height(data)/2000));
quiver_data = data(1:step:end, :);

quiver(quiver_data.r_um, quiver_data.z_um, ...
       quiver_data.u_um, quiver_data.w_um, 2, ...
       'LineWidth', 1.2, 'Color', 'b', 'MaxHeadSize', 0.5);

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$z\ (\mu\mathrm{m})$', 'FontSize', 18);
title('Velocity field $\mathbf{v}(r,z)$', 'FontSize', 18);

axis equal tight;
grid on;
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'velocity_vectors';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% FIGURE 10: Streamlines
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Interpolate velocities on grid
F_u = scatteredInterpolant(data.r_um, data.z_um, data.u_um, 'linear', 'none');
F_w = scatteredInterpolant(data.r_um, data.z_um, data.w_um, 'linear', 'none');

U_grid = F_u(R_grid, Z_grid);
W_grid = F_w(R_grid, Z_grid);

% Speed for background color
speed_grid = sqrt(U_grid.^2 + W_grid.^2);

% Plot streamlines
contourf(R_grid, Z_grid, speed_grid, 20, 'LineStyle', 'none', 'FaceAlpha', 0.6);
hold on;
streamslice(R_grid, Z_grid, U_grid, W_grid, 2);

colormap(parula);
cb = colorbar;
cb.Label.String = 'Speed ($\mu$m/s)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$z\ (\mu\mathrm{m})$', 'FontSize', 18);
title('Velocity streamlines', 'FontSize', 18);

axis equal tight;
grid on;
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'streamlines';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% ========================================================================
%% SECTION 7: FORCE ANALYSIS
%% ========================================================================

fprintf('\n=== FORCE ANALYSIS ===\n');

% Force magnitudes
data.F_r_total = data.F_TPr + data.F_Dr;
data.F_z_total = data.F_TPz + data.F_Dz + data.F_g;
data.F_total_mag = sqrt(data.F_r_total.^2 + data.F_z_total.^2);

fprintf('Force components (N):\n');
fprintf('  Thermophoresis (radial): %.2e to %.2e\n', min(data.F_TPr), max(data.F_TPr));
fprintf('  Drag (radial): %.2e to %.2e\n', min(data.F_Dr), max(data.F_Dr));
fprintf('  Thermophoresis (vertical): %.2e to %.2e\n', min(data.F_TPz), max(data.F_TPz));
fprintf('  Drag (vertical): %.2e to %.2e\n', min(data.F_Dz), max(data.F_Dz));
fprintf('  Gravity: %.2e N (constant)\n', data.F_g(1));

fprintf('\nTotal forces:\n');
fprintf('  F_r (net radial): %.2e to %.2e N\n', min(data.F_r_total), max(data.F_r_total));
fprintf('  F_z (net vertical): %.2e to %.2e N\n', min(data.F_z_total), max(data.F_z_total));

% Check near-surface behavior
near_surface = data(data.z_um < 10, :);
fprintf('\nNear surface (z < 10 μm):\n');
fprintf('  Mean F_z: %.2e N\n', mean(near_surface.F_z_total));
if mean(near_surface.F_z_total) < 0
    fprintf('  → Net force is DOWNWARD (particles pushed to substrate)\n');
else
    fprintf('  → Net force is UPWARD\n');
end

%% FIGURE 11: Force components vs height (at center)
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Extract center data (r < 5 μm)
center_data = data(data.r_um < 5, :);
center_data = sortrows(center_data, 'z_um');

plot(center_data.F_TPz * 1e15, center_data.z_um, 'o-', 'LineWidth', 2.5, ...
     'MarkerSize', 6, 'DisplayName', 'Thermophoresis');
hold on;
plot(center_data.F_Dz * 1e15, center_data.z_um, 's-', 'LineWidth', 2.5, ...
     'MarkerSize', 6, 'DisplayName', 'Drag');
plot(center_data.F_g(1) * 1e15 * ones(size(center_data.z_um)), center_data.z_um, ...
     '^-', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Gravity');
plot(center_data.F_z_total * 1e15, center_data.z_um, 'k-', 'LineWidth', 3, ...
     'DisplayName', 'Net force');

xline(0, 'k--', 'LineWidth', 1.5);
yline(10, 'r:', 'LineWidth', 2, 'Alpha', 0.5);

xlabel('$F_z$ (fN)', 'FontSize', 18);
ylabel('$z$ ($\mu$m)', 'FontSize', 18);
title('Vertical force components (center, $r < 5$ $\mu$m)', 'FontSize', 18);

legend('Location', 'southeast', 'FontSize', 14);
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;
ylim([0, min(50, max(center_data.z_um))]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'forces_vertical_center';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% FIGURE 12: Net force spatial distribution (near surface)
figure('Position', [100, 100, 900, 700], 'Color', 'w');

near_surf_plot = data(data.z_um < 20 & data.z_um > 0, :);

scatter(near_surf_plot.r_um, near_surf_plot.z_um, 15, ...
        near_surf_plot.F_z_total * 1e15, 'filled');
colormap(redblue_colormap);
max_F = max(abs(near_surf_plot.F_z_total * 1e15));
caxis([-max_F, max_F]);
cb = colorbar;
cb.Label.String = '$F_z$ (fN)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$z\ (\mu\mathrm{m})$', 'FontSize', 18);
title('Net vertical force (near surface)', 'FontSize', 18);

axis equal tight;
grid on;
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'forces_net_spatial';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% FIGURE 13: Force magnitude total
figure('Position', [100, 100, 900, 700], 'Color', 'w');

scatter(data.r_um(idx_plot), data.z_um(idx_plot), 8, ...
        data.F_total_mag(idx_plot) * 1e15, 'filled');
colormap(hot);
cb = colorbar;
cb.Label.String = '$|\mathbf{F}|$ (fN)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$z\ (\mu\mathrm{m})$', 'FontSize', 18);
title('Total force magnitude', 'FontSize', 18);

axis equal tight;
grid on;
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'forces_total_magnitude';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);
%% ========================================================================
%% ADDITIONAL FIGURES: FORCE FIELD ANALYSIS
%% ========================================================================
%% Add these after SECTION 7 in your COMSOL_analysis.m script
%% ========================================================================

%% FIGURE 14: Total Force Vector Field (Quiver)
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Subsample for force vectors
step_force = max(1, floor(height(data)/2000));
force_data = data(1:step_force:end, :);

% Calculate force components
force_data.F_r_total = force_data.F_TPr + force_data.F_Dr;
force_data.F_z_total = force_data.F_TPz + force_data.F_Dz + force_data.F_g;

% Convert to fN for better visualization
force_data.F_r_fN = force_data.F_r_total * 1e15;
force_data.F_z_fN = force_data.F_z_total * 1e15;

% Background: force magnitude
scatter(data.r_um(idx_plot), data.z_um(idx_plot), 8, ...
        data.F_total_mag(idx_plot) * 1e15, 'filled', 'MarkerFaceAlpha', 0.4);
hold on;

% Overlay force vectors
quiver(force_data.r_um, force_data.z_um, ...
       force_data.F_r_fN, force_data.F_z_fN, 2, ...
       'LineWidth', 1.5, 'Color', 'k', 'MaxHeadSize', 0.5);

colormap(hot);
cb = colorbar;
cb.Label.String = '$|\mathbf{F}_{\mathrm{total}}|$ (fN)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$z\ (\mu\mathrm{m})$', 'FontSize', 18);
title('Total force field $\mathbf{F}_{\mathrm{total}} = \mathbf{F}_{\mathrm{TP}} + \mathbf{F}_{\mathrm{D}} + \mathbf{F}_g$', 'FontSize', 18);

axis equal tight;
grid on;
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;

% Add annotation
text(0.05, 0.95, sprintf('Gravity: %.2f fN (downward)', abs(data.F_g(1)*1e15)), ...
     'Units', 'normalized', 'FontSize', 12, 'BackgroundColor', 'w', ...
     'EdgeColor', 'k', 'VerticalAlignment', 'top');

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);

fileName = 'force_field_total_quiver';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% FIGURE 15: Vertical Force F_z vs r (compared with v_z)
figure('Position', [100, 100, 1400, 600], 'Color', 'w');

%% Panel A: F_z (radially averaged at surface)
subplot(1,2,1);

% Calculate F_z_total for surface data if not already present
if ~ismember('F_z_total', surf_data.Properties.VariableNames)
    surf_data.F_z_total = surf_data.F_TPz + surf_data.F_Dz + surf_data.F_g;
end

% Extract surface data and bin by radius
Fz_mean = zeros(n_bins, 1);
Fz_sem = zeros(n_bins, 1);

for k = 1:n_bins
    idx = surf_data.r_um >= r_edges(k) & surf_data.r_um < r_edges(k+1);
    if sum(idx) > 5
        Fz_mean(k) = mean(real(surf_data.F_z_total(idx))) * 1e15;  % Convert to fN, use real part
        Fz_sem(k) = std(real(surf_data.F_z_total(idx))) * 1e15 / sqrt(sum(idx));
    else
        Fz_mean(k) = NaN;
        Fz_sem(k) = NaN;
    end
end

% Plot individual cells (use real parts to avoid warnings)
scatter(surf_data.r_um, real(surf_data.F_z_total) * 1e15, 15, 'b', ...
        'filled', 'MarkerFaceAlpha', 0.2);
hold on;

% Plot binned average with error bars
errorbar(r_centers, Fz_mean, Fz_sem, 'ro-', 'LineWidth', 2.5, ...
         'MarkerSize', 10, 'MarkerFaceColor', 'r', 'CapSize', 10);

yline(0, 'k--', 'LineWidth', 1.5);
yline(-abs(real(data.F_g(1))*1e15), 'g--', 'LineWidth', 1.5, ...
      'DisplayName', sprintf('Gravity (%.2f fN)', abs(real(data.F_g(1))*1e15)));

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$\langle F_z \rangle$ (fN)', 'FontSize', 18);
title('(a) Net vertical force', 'FontSize', 18);

legend('Grid cells', 'Binned average', 'Zero line', 'Gravity', ...
       'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;
xlim([0, max(r_centers)]);

% Add shaded region for negative (downward) forces
if any(Fz_mean < 0)
    y_fill = [min(Fz_mean-Fz_sem), 0];
    fill([0, max(r_centers), max(r_centers), 0], ...
         [y_fill(1), y_fill(1), 0, 0], ...
         'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

%% Panel B: v_z (from earlier calculation, for comparison)
subplot(1,2,2);

% Re-plot v_z for comparison
scatter(surf_data.r_um, surf_data.w_um, 15, 'b', ...
        'filled', 'MarkerFaceAlpha', 0.2);
hold on;
errorbar(r_centers, wz_mean, wz_sem, 'ro-', 'LineWidth', 2.5, ...
         'MarkerSize', 10, 'MarkerFaceColor', 'r', 'CapSize', 10);
yline(0, 'k--', 'LineWidth', 1.5);

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 18);
ylabel('$\langle v_z \rangle$ ($\mu$m/s)', 'FontSize', 18);
title('(b) Vertical velocity', 'FontSize', 18);

legend('Grid cells', 'Binned average', 'Zero line', ...
       'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 16);
box on;
xlim([0, max(r_centers)]);

% Add shaded region for upward velocities
if any(wz_mean > 0)
    y_fill = [0, max(wz_mean+wz_sem)];
    fill([0, max(r_centers), max(r_centers), 0], ...
         [0, 0, y_fill(2), y_fill(2)], ...
         'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 14, 6]);

fileName = 'Fz_vs_vz_comparison';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% FIGURE 16: Force Components Decomposition (Radial Profiles at Surface)
figure('Position', [100, 100, 1400, 1000], 'Color', 'w');

%% Calculate all force components vs radius
F_TPr_mean = zeros(n_bins, 1);
F_Dr_mean = zeros(n_bins, 1);
F_TPz_mean = zeros(n_bins, 1);
F_Dz_mean = zeros(n_bins, 1);

for k = 1:n_bins
    idx = surf_data.r_um >= r_edges(k) & surf_data.r_um < r_edges(k+1);
    if sum(idx) > 5
        F_TPr_mean(k) = mean(real(surf_data.F_TPr(idx))) * 1e15;
        F_Dr_mean(k) = mean(real(surf_data.F_Dr(idx))) * 1e15;
        F_TPz_mean(k) = mean(real(surf_data.F_TPz(idx))) * 1e15;
        F_Dz_mean(k) = mean(real(surf_data.F_Dz(idx))) * 1e15;
    else
        F_TPr_mean(k) = NaN;
        F_Dr_mean(k) = NaN;
        F_TPz_mean(k) = NaN;
        F_Dz_mean(k) = NaN;
    end
end

%% Panel A: Radial force components
subplot(2,2,1);
hold on;
plot(r_centers, F_TPr_mean, 'o-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
     'DisplayName', 'Thermophoresis');
plot(r_centers, F_Dr_mean, 's-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
     'DisplayName', 'Drag');
plot(r_centers, F_TPr_mean + F_Dr_mean, 'k-', 'LineWidth', 3, ...
     'DisplayName', 'Total $F_r$');
yline(0, 'k--', 'LineWidth', 1.5);

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 16);
ylabel('$F_r$ (fN)', 'FontSize', 16);
title('(a) Radial force components', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);
xlim([0, max(r_centers)]);

%% Panel B: Vertical force components
subplot(2,2,2);
hold on;
plot(r_centers, F_TPz_mean, 'o-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
     'DisplayName', 'Thermophoresis');
plot(r_centers, F_Dz_mean, 's-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
     'DisplayName', 'Drag');
plot(r_centers, ones(size(r_centers)) * real(data.F_g(1)) * 1e15, '^-', ...
     'LineWidth', 2.5, 'MarkerSize', 8, 'DisplayName', 'Gravity');
plot(r_centers, Fz_mean, 'k-', 'LineWidth', 3, 'DisplayName', 'Total $F_z$');
yline(0, 'k--', 'LineWidth', 1.5);

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 16);
ylabel('$F_z$ (fN)', 'FontSize', 16);
title('(b) Vertical force components', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);
xlim([0, max(r_centers)]);

%% Panel C: Ratio Thermophoresis/Drag (radial)
subplot(2,2,3);
ratio_r = abs(F_TPr_mean ./ F_Dr_mean);
plot(r_centers, ratio_r, 'o-', 'LineWidth', 2.5, 'MarkerSize', 8, 'Color', [0.8, 0.2, 0.2]);
yline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Equal contribution');

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 16);
ylabel('$|F_{\mathrm{TP},r}| / |F_{\mathrm{D},r}|$', 'FontSize', 16);
title('(c) Thermophoresis/Drag ratio (radial)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);
xlim([0, max(r_centers)]);
% ylim([0, max(ratio_r)*1.1]);

%% Panel D: Ratio Thermophoresis/Drag (vertical)
subplot(2,2,4);
ratio_z = abs(F_TPz_mean ./ F_Dz_mean);
plot(r_centers, ratio_z, 's-', 'LineWidth', 2.5, 'MarkerSize', 8, 'Color', [0.2, 0.4, 0.8]);
yline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Equal contribution');

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 16);
ylabel('$|F_{\mathrm{TP},z}| / |F_{\mathrm{D},z}|$', 'FontSize', 16);
title('(d) Thermophoresis/Drag ratio (vertical)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);
xlim([0, max(r_centers)]);
ylim([0, max(ratio_z)*1.1]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 14, 10]);

fileName = 'force_components_decomposition';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

%% Update summary with force ratios
fprintf('\n=== FORCE RATIOS ===\n');
fprintf('Radial:\n');
fprintf('  Mean |F_TP| / |F_D|: %.2f\n', nanmean(ratio_r));
fprintf('Vertical:\n');
fprintf('  Mean |F_TP| / |F_D|: %.2f\n', nanmean(ratio_z));

if nanmean(ratio_r) < 0.1
    fprintf('  → Drag dominates radially (convection-driven)\n');
elseif nanmean(ratio_r) > 10
    fprintf('  → Thermophoresis dominates radially\n');
else
    fprintf('  → Mixed regime radially\n');
end

if nanmean(ratio_z) < 0.1
    fprintf('  → Drag dominates vertically (convection-driven)\n');
elseif nanmean(ratio_z) > 10
    fprintf('  → Thermophoresis dominates vertically\n');
else
    fprintf('  → Mixed regime vertically\n');
end

fprintf('\n');

%% FIGURE 17: Correlation between ∇T and v_r (Mechanism Validation)
figure('Position', [100, 100, 1400, 600], 'Color', 'w');

%% Panel A: Both profiles on same plot (normalized)
subplot(1,2,1);
yyaxis left
plot(r_centers, grad_T_surf / max(grad_T_surf), 'r-', 'LineWidth', 3);
ylabel('$|\nabla T| / |\nabla T|_{\max}$', 'FontSize', 16);
set(gca, 'YColor', 'r');

yyaxis right
plot(r_centers, abs(vr_mean) / max(abs(vr_mean)), 'b-', 'LineWidth', 3);
ylabel('$|v_r| / |v_r|_{\max}$', 'FontSize', 16);
set(gca, 'YColor', 'b');

xlabel('$r\ (\mu\mathrm{m})$', 'FontSize', 16);
title('(a) Normalized profiles: $v_r \propto \nabla T$', 'FontSize', 16, 'FontWeight', 'bold');
legend('$|\nabla T|$', '$|v_r|$', 'Location', 'best', 'FontSize', 12);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);
xlim([0, max(r_centers)]);

%% Panel B: Direct correlation plot
subplot(1,2,2);

% Remove NaN values for correlation
valid_idx = ~isnan(grad_T_surf) & ~isnan(vr_mean);
grad_T_valid = grad_T_surf(valid_idx);
vr_valid = abs(vr_mean(valid_idx));

% Scatter plot
scatter(grad_T_valid, vr_valid, 80, r_centers(valid_idx), 'filled', 'MarkerEdgeColor', 'k');
hold on;

% Linear fit
p = polyfit(grad_T_valid, vr_valid, 1);
grad_T_fit = linspace(min(grad_T_valid), max(grad_T_valid), 100);
vr_fit = polyval(p, grad_T_fit);
plot(grad_T_fit, vr_fit, 'k--', 'LineWidth', 2);

% Calculate R²
vr_pred = polyval(p, grad_T_valid);
SS_res = sum((vr_valid - vr_pred).^2);
SS_tot = sum((vr_valid - mean(vr_valid)).^2);
R_squared = 1 - SS_res/SS_tot;

xlabel('$|\nabla T|$ (K/$\mu$m)', 'FontSize', 16);
ylabel('$|v_r|$ ($\mu$m/s)', 'FontSize', 16);
title('(b) Correlation: $v_r \sim (\\beta g h^3 / \\nu) \\nabla T$', 'FontSize', 16, 'FontWeight', 'bold');

% Add fit info
text(0.05, 0.95, sprintf('Linear fit:\n$v_r = %.3f \\times |\\nabla T| + %.3f$\n$R^2 = %.3f$', ...
     p(1), p(2), R_squared), ...
     'Units', 'normalized', 'FontSize', 12, 'BackgroundColor', 'w', ...
     'EdgeColor', 'k', 'VerticalAlignment', 'top', 'Interpreter', 'latex');

colormap(parula);
cb = colorbar;
cb.Label.String = '$r$ ($\mu$m)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;

grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 14, 6]);

fileName = 'correlation_gradT_vr';
saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
saveas(gcf, fullfile(fig_dir, [fileName '.png']));
fprintf('Saved: %s\n', fileName);

fprintf('\n=== CORRELATION ANALYSIS ===\n');
fprintf('Linear fit: v_r = %.3f × |∇T| + %.3f\n', p(1), p(2));
fprintf('R² = %.3f\n', R_squared);
fprintf('Slope interpretation:\n');
fprintf('  Predicted: v ~ (βgh³/ν) × ∇T\n');
fprintf('  With β=2.1e-4 K⁻¹, g=9.81 m/s², h=100μm, ν=1e-6 m²/s:\n');
fprintf('  Expected slope: ~20 μm/s per (K/μm)\n');
fprintf('  Measured slope: %.1f μm/s per (K/μm)\n', p(1));
if abs(p(1) - 20) < 10
    fprintf('  → EXCELLENT agreement with convection theory! ✓\n');
elseif abs(p(1) - 20) < 50
    fprintf('  → Reasonable agreement with convection theory\n');
else
    fprintf('  → Consider temperature scaling adjustment\n');
end

fprintf('\n');
%% ========================================================================
%% SECTION 8: SUMMARY METRICS
%% ========================================================================

fprintf('\n=== SUMMARY METRICS ===\n');

summary = struct();
summary.simulation_name = sim_name;
summary.profile_type = profile_type;
summary.w_um = w_um;
summary.d_shape = d_shape;

summary.temperature.Delta_T = Delta_T;
summary.temperature.T_min = min(data.T);
summary.temperature.T_max = max(data.T);
summary.temperature.grad_max = max(grad_T_mag(:));

summary.velocity.u_r_min = min(data.u_um);
summary.velocity.u_r_max = max(data.u_um);
summary.velocity.u_z_min = min(data.w_um);
summary.velocity.u_z_max = max(data.w_um);
summary.velocity.v_mag_max = max(data.v_mag);
summary.velocity.r_max_inflow = r_min;

summary.forces.F_TP_r_range = [min(data.F_TPr), max(data.F_TPr)];
summary.forces.F_D_r_range = [min(data.F_Dr), max(data.F_Dr)];
summary.forces.F_TP_z_range = [min(data.F_TPz), max(data.F_TPz)];
summary.forces.F_D_z_range = [min(data.F_Dz), max(data.F_Dz)];
summary.forces.F_g = data.F_g(1);
summary.forces.F_z_mean_near_surface = mean(near_surface.F_z_total);

% Péclet number
D_Brownian = 0.43;  % μm²/s for 1 μm particle at 20°C
Pe_max = max(data.v_mag) * 2 / D_Brownian;  % using 2 μm as length scale
summary.transport.Pe_max = Pe_max;
summary.transport.D_Brownian = D_Brownian;

fprintf('\nSimulation: %s\n', sim_name);
fprintf('Profile: %s (w=%.1f μm, d=%.1f)\n', profile_type, w_um, d_shape);
fprintf('\nTemperature:\n');
fprintf('  ΔT = %.2f K\n', Delta_T);
fprintf('  Max gradient: %.3f K/μm\n', max(grad_T_mag(:)));
fprintf('\nVelocity:\n');
fprintf('  Max inflow: %.3f μm/s at r = %.1f μm\n', vr_min, r_min);
fprintf('  Max upwelling: %.3f μm/s\n', max(data.w_um));
fprintf('  Max speed: %.3f μm/s\n', max(data.v_mag));
fprintf('\nForces:\n');
fprintf('  Mean F_z (z<10μm): %.2e N\n', summary.forces.F_z_mean_near_surface);
% fprintf('  Direction: %s\n', char("DOWNWARD" * (summary.forces.F_z_mean_near_surface < 0) + ...
                                   % "UPWARD" * (summary.forces.F_z_mean_near_surface >= 0)));
% fprintf('\nTransport:\n');
% fprintf('  Max Péclet: %.1f\n', Pe_max);
% fprintf('  Regime: %s\n', char("Advection-dominated" * (Pe_max > 10) + ...
%                                 "Diffusion-dominated" * (Pe_max <= 10)));

%% Save data
save(fullfile(data_dir, 'analysis_results.mat'), 'summary', 'data', ...
     'r_centers', 'vr_mean', 'vr_sem', 'wz_mean', 'wz_sem', ...
     'T_mean', 'T_std', 'grad_T_surf');

fprintf('\n========================================================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('Results saved to: %s\n', results_base);
fprintf('Generated %d figures\n', 13);
fprintf('========================================================================\n');

%% ========================================================================
%% HELPER FUNCTION: Red-Blue Colormap
%% ========================================================================
function c = redblue_colormap(m)
    if nargin < 1
        m = size(get(gcf,'colormap'),1);
    end
    r = (0:m-1)'/max(m-1,1);
    c = [r.^0.5, (1-r).^0.5, 1-r];
    c = flipud(c);
end