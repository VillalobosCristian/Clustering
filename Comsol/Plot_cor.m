%% ========================================================================
%% PUBLICATION QUALITY COMSOL DATA VISUALIZATION
%% Velocity field, Force fields, and Radial velocity profile
%% ========================================================================
%
% Generates publication-ready figures for optothermal colloidal assembly
%
% Output figures:
%   - Fig1: Temperature field with velocity streamlines
%   - Fig2: Force vector fields (Thermophoretic, Drag, Total)
%   - Fig3: Radial velocity profile <v_r>
%   - Fig4: 2D Force magnitude maps
%
% Author: Optothermal Colloidal Assembly Research
% Date: December 2025
%% ========================================================================

clear; close all; clc;

%% ========================================================================
%% CONFIGURATION
%% ========================================================================

% File selection
filename = 'Fields_35um.csv';  % Change to: Fields_25um.csv, Fields_35um.csv, Fields_75um.csv
R_illum = 35;                  % Illumination radius [μm] - match to file!

% Figure export settings
export_figures = true;
export_format = 'png';         % 'png', 'pdf', 'eps', 'svg'
export_dpi = 300;

% Plotting parameters
plot_params.r_max = 80;        % Max r for plots [μm]
plot_params.z_max = 25;        % Max z for plots [μm]
plot_params.grid_res = 150;    % Grid resolution
plot_params.skip_vec = 8;      % Vector skip for quiver plots
plot_params.z_profile = 5;     % Height for radial profiles [μm]

% Set default figure properties for publication quality
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesLineWidth', 1.2);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultAxesTickDir', 'out');
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');

%% ========================================================================
%% LOAD DATA
%% ========================================================================

fprintf('Loading data from: %s\n', filename);

if ~exist(filename, 'file')
    error('File not found: %s\nPlease check the filename and path.', filename);
end

% Read header info
fid = fopen(filename, 'r');
for i = 1:9
    fgetl(fid);  % Skip header lines
end
fclose(fid);

% Load numerical data
data = readmatrix(filename, 'NumHeaderLines', 9);

% Extract columns
r = data(:, 1);        % r coordinate [m]
z = data(:, 2);        % z coordinate [m]
T = data(:, 3);        % Temperature [K]
F_TPr = data(:, 4);    % Thermophoretic force r [N]
F_Dr = data(:, 5);     % Drag force r [N]
F_TPz = data(:, 6);    % Thermophoretic force z [N]
F_Dz = data(:, 7);     % Drag force z [N]
F_g = data(:, 8);      % Gravitational force [N]
u = data(:, 9);        % Velocity r [m/s]
w = data(:, 10);       % Velocity z [m/s]

% Handle complex numbers (take real part)
F_TPr = real(F_TPr);
F_TPz = real(F_TPz);
F_Dr = real(F_Dr);
F_Dz = real(F_Dz);
F_g = real(F_g);

% Convert units
r_um = r * 1e6;        % [μm]
z_um = z * 1e6;        % [μm]
u_um = u * 1e6;        % [μm/s]
w_um = w * 1e6;        % [μm/s]
Delta_T = T - 293.15;  % Temperature rise [K]

% Calculate derived quantities
vel_mag = sqrt(u.^2 + w.^2);           % Velocity magnitude [m/s]
F_TP_mag = sqrt(F_TPr.^2 + F_TPz.^2);  % Thermophoretic force magnitude [N]
F_D_mag = sqrt(F_Dr.^2 + F_Dz.^2);     % Drag force magnitude [N]

% Total force (sum of all components)
F_total_r = F_TPr + F_Dr;              % Total r-component [N]
F_total_z = F_TPz + F_Dz + F_g;        % Total z-component [N]
F_total_mag = sqrt(F_total_r.^2 + F_total_z.^2);

fprintf('  Loaded %d data points\n', length(r));
fprintf('  r: [%.1f, %.1f] μm\n', min(r_um), max(r_um));
fprintf('  z: [%.1f, %.1f] μm\n', min(z_um), max(z_um));
fprintf('  ΔT_max: %.2f K\n', max(Delta_T));

%% ========================================================================
%% CREATE INTERPOLATION GRID
%% ========================================================================

fprintf('Creating interpolation grid...\n');

% Define grid boundaries
r_min = 0;
r_max_grid = min(max(r_um), plot_params.r_max);
z_min = 0;
z_max_grid = min(max(z_um), plot_params.z_max);

% Create regular grid
r_grid = linspace(r_min, r_max_grid, plot_params.grid_res);
z_grid = linspace(z_min + 0.5, z_max_grid, plot_params.grid_res);
[R_grid, Z_grid] = meshgrid(r_grid, z_grid);

% Interpolate all fields
T_grid = griddata(r_um, z_um, T, R_grid, Z_grid, 'natural');
Delta_T_grid = griddata(r_um, z_um, Delta_T, R_grid, Z_grid, 'natural');
u_grid = griddata(r_um, z_um, u_um, R_grid, Z_grid, 'natural');
w_grid = griddata(r_um, z_um, w_um, R_grid, Z_grid, 'natural');
vel_mag_grid = griddata(r_um, z_um, vel_mag*1e6, R_grid, Z_grid, 'natural');

% Forces in fN
F_TPr_grid = griddata(r_um, z_um, F_TPr*1e15, R_grid, Z_grid, 'natural');
F_TPz_grid = griddata(r_um, z_um, F_TPz*1e15, R_grid, Z_grid, 'natural');
F_Dr_grid = griddata(r_um, z_um, F_Dr*1e15, R_grid, Z_grid, 'natural');
F_Dz_grid = griddata(r_um, z_um, F_Dz*1e15, R_grid, Z_grid, 'natural');
F_total_r_grid = griddata(r_um, z_um, F_total_r*1e15, R_grid, Z_grid, 'natural');
F_total_z_grid = griddata(r_um, z_um, F_total_z*1e15, R_grid, Z_grid, 'natural');
F_total_mag_grid = griddata(r_um, z_um, F_total_mag*1e15, R_grid, Z_grid, 'natural');

%% ========================================================================
%% FIGURE 1: TEMPERATURE FIELD WITH VELOCITY STREAMLINES
%% ========================================================================

fprintf('Generating Figure 1: Temperature & Velocity...\n');

fig1 = figure('Name', 'Temperature & Velocity', 'Position', [50, 50, 700, 550]);
set(fig1, 'Color', 'white');

% Temperature contour
contourf(R_grid, Z_grid, Delta_T_grid, 30, 'LineStyle', 'none');
hold on;

% Colormap and colorbar
colormap(gca, parula);
cb = colorbar;
ylabel(cb, '$\Delta T$ (K)', 'Interpreter', 'latex', 'FontSize', 13);
caxis([0, max(Delta_T)]);

% Velocity streamlines
h_stream = streamslice(R_grid, Z_grid, u_grid, w_grid, 2);
set(h_stream, 'Color', 'w', 'LineWidth', 0.6);

% Illumination zone marker
xline(R_illum, 'r--', 'LineWidth', 2);
text(R_illum + 2, z_max_grid - 3, sprintf('$R = %d$ $\\mu$m', R_illum), ...
    'Color', 'r', 'FontSize', 11, 'Interpreter', 'latex');

% Labels
xlabel('$r$ ($\mu$m)', 'FontSize', 14);
ylabel('$z$ ($\mu$m)', 'FontSize', 14);

axis equal;
xlim([0, plot_params.r_max]);
ylim([0, plot_params.z_max]);
set(gca, 'Layer', 'top');

if export_figures
    exportgraphics(fig1, ['Fig1_Temperature_Velocity.' export_format], 'Resolution', export_dpi);
    fprintf('  Saved: Fig1_Temperature_Velocity.%s\n', export_format);
end

%% ========================================================================
%% FIGURE 2: FORCE VECTOR FIELDS (3 PANELS)
%% ========================================================================

fprintf('Generating Figure 2: Force Fields...\n');

fig2 = figure('Name', 'Force Fields', 'Position', [100, 100, 1500, 450]);
set(fig2, 'Color', 'white');
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

skip = plot_params.skip_vec;

% Colors
color_TP = [0.8, 0.2, 0.2];    % Red for thermophoretic
color_D = [0.2, 0.4, 0.8];     % Blue for drag
color_tot = [0.1, 0.1, 0.1];   % Black for total

% --- Panel (a): Thermophoretic Force ---
nexttile;
quiver(R_grid(1:skip:end, 1:skip:end), Z_grid(1:skip:end, 1:skip:end), ...
       F_TPr_grid(1:skip:end, 1:skip:end), F_TPz_grid(1:skip:end, 1:skip:end), ...
       2, 'Color', color_TP, 'LineWidth', 1.0);
hold on;
xline(R_illum, 'k--', 'LineWidth', 1.5);

xlabel('$r$ ($\mu$m)', 'FontSize', 14);
ylabel('$z$ ($\mu$m)', 'FontSize', 14);
title('(a) $\mathbf{F}_{\mathrm{TP}}$', 'FontSize', 14, 'FontWeight', 'normal');
axis equal;
xlim([0, plot_params.r_max]);
ylim([0, plot_params.z_max]);
set(gca, 'Layer', 'top');

% --- Panel (b): Drag Force ---
nexttile;
quiver(R_grid(1:skip:end, 1:skip:end), Z_grid(1:skip:end, 1:skip:end), ...
       F_Dr_grid(1:skip:end, 1:skip:end), F_Dz_grid(1:skip:end, 1:skip:end), ...
       2, 'Color', color_D, 'LineWidth', 1.0);
hold on;
xline(R_illum, 'k--', 'LineWidth', 1.5);

xlabel('$r$ ($\mu$m)', 'FontSize', 14);
ylabel('$z$ ($\mu$m)', 'FontSize', 14);
title('(b) $\mathbf{F}_{\mathrm{D}}$', 'FontSize', 14, 'FontWeight', 'normal');
axis equal;
xlim([0, plot_params.r_max]);
ylim([0, plot_params.z_max]);
set(gca, 'Layer', 'top');

% --- Panel (c): Total Force ---
nexttile;
quiver(R_grid(1:skip:end, 1:skip:end), Z_grid(1:skip:end, 1:skip:end), ...
       F_total_r_grid(1:skip:end, 1:skip:end), F_total_z_grid(1:skip:end, 1:skip:end), ...
       2, 'Color', color_tot, 'LineWidth', 1.0);
hold on;
xline(R_illum, 'r--', 'LineWidth', 1.5);

xlabel('$r$ ($\mu$m)', 'FontSize', 14);
ylabel('$z$ ($\mu$m)', 'FontSize', 14);
title('(c) $\mathbf{F}_{\mathrm{total}} = \mathbf{F}_{\mathrm{TP}} + \mathbf{F}_{\mathrm{D}} + \mathbf{F}_g$', ...
      'FontSize', 14, 'FontWeight', 'normal');
axis equal;
xlim([0, plot_params.r_max]);
ylim([0, plot_params.z_max]);
set(gca, 'Layer', 'top');

if export_figures
    exportgraphics(fig2, ['Fig2_Force_Fields.' export_format], 'Resolution', export_dpi);
    fprintf('  Saved: Fig2_Force_Fields.%s\n', export_format);
end

%% ========================================================================
%% FIGURE 3: RADIAL VELOCITY PROFILE <v_r>
%% ========================================================================

fprintf('Generating Figure 3: Radial Velocity Profile...\n');

% Extract near-surface data
z_thresh = plot_params.z_profile;
idx_surf = z_um < z_thresh;

% Bin by r
r_bins = 0:2:plot_params.r_max;
r_centers = (r_bins(1:end-1) + r_bins(2:end)) / 2;
n_bins = length(r_centers);

v_r_mean = zeros(1, n_bins);
v_r_std = zeros(1, n_bins);

for i = 1:n_bins
    idx_bin = idx_surf & (r_um >= r_bins(i)) & (r_um < r_bins(i+1));
    if sum(idx_bin) > 0
        v_r_mean(i) = mean(u_um(idx_bin));
        v_r_std(i) = std(u_um(idx_bin));
    else
        v_r_mean(i) = NaN;
        v_r_std(i) = NaN;
    end
end

% Create figure
fig3 = figure('Name', 'Radial Velocity Profile', 'Position', [150, 150, 600, 450]);
set(fig3, 'Color', 'white');

hold on;

% Shaded regions for inflow/outflow
fill([r_centers, fliplr(r_centers)], ...
     [max(v_r_mean, 0), zeros(1, n_bins)], ...
     [0.9, 0.7, 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
fill([r_centers, fliplr(r_centers)], ...
     [min(v_r_mean, 0), zeros(1, n_bins)], ...
     [0.7, 0.8, 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Main curve
plot(r_centers, v_r_mean, 'b-', 'LineWidth', 2.5);
plot(r_centers, v_r_mean, 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');

% Reference lines
yline(0, 'k-', 'LineWidth', 0.8);
xline(R_illum, 'r--', 'LineWidth', 2);

% Labels and annotations
xlabel('$r$ ($\mu$m)', 'FontSize', 14);
ylabel('$\langle v_r \rangle$ ($\mu$m/s)', 'FontSize', 14);
title(sprintf('Radial velocity profile ($z < %d$ $\\mu$m)', z_thresh), ...
      'FontSize', 14, 'FontWeight', 'normal');

% Legend
text(R_illum + 2, max(v_r_mean)*0.8, sprintf('$R = %d$ $\\mu$m', R_illum), ...
    'Color', 'r', 'FontSize', 11, 'Interpreter', 'latex');

% Find zero crossing
idx_cross = find(diff(sign(v_r_mean)) ~= 0, 1);
if ~isempty(idx_cross)
    r_cross = r_centers(idx_cross);
    text(r_cross, 0.002, sprintf('$r_0 \\approx %.0f$ $\\mu$m', r_cross), ...
        'FontSize', 10, 'Interpreter', 'latex', 'VerticalAlignment', 'bottom');
end

grid on;
xlim([0, plot_params.r_max]);
set(gca, 'Layer', 'top');

if export_figures
    exportgraphics(fig3, ['Fig3_Radial_Velocity.' export_format], 'Resolution', export_dpi);
    fprintf('  Saved: Fig3_Radial_Velocity.%s\n', export_format);
end

%% ========================================================================
%% FIGURE 4: TOTAL FORCE ON TEMPERATURE BACKGROUND
%% ========================================================================

fprintf('Generating Figure 4: Total Force with Temperature...\n');

fig4 = figure('Name', 'Total Force on Temperature', 'Position', [200, 200, 700, 550]);
set(fig4, 'Color', 'white');

% Temperature background
contourf(R_grid, Z_grid, Delta_T_grid, 30, 'LineStyle', 'none');
hold on;

% Colormap
colormap(gca, hot);
cb = colorbar;
ylabel(cb, '$\Delta T$ (K)', 'Interpreter', 'latex', 'FontSize', 13);
caxis([0, max(Delta_T)]);

% Force vectors (white for visibility)
quiver(R_grid(1:skip:end, 1:skip:end), Z_grid(1:skip:end, 1:skip:end), ...
       F_total_r_grid(1:skip:end, 1:skip:end), F_total_z_grid(1:skip:end, 1:skip:end), ...
       2, 'w', 'LineWidth', 1.0);

% Illumination zone
xline(R_illum, 'c--', 'LineWidth', 2);

xlabel('$r$ ($\mu$m)', 'FontSize', 14);
ylabel('$z$ ($\mu$m)', 'FontSize', 14);

axis equal;
xlim([0, plot_params.r_max]);
ylim([0, plot_params.z_max]);
set(gca, 'Layer', 'top');

if export_figures
    exportgraphics(fig4, ['Fig4_Force_Temperature.' export_format], 'Resolution', export_dpi);
    fprintf('  Saved: Fig4_Force_Temperature.%s\n', export_format);
end

%% ========================================================================
%% FIGURE 5: COMBINED OVERVIEW (2x2)
%% ========================================================================

fprintf('Generating Figure 5: Combined Overview...\n');

fig5 = figure('Name', 'Combined Overview', 'Position', [250, 50, 1200, 900]);
set(fig5, 'Color', 'white');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel (a): Temperature + Streamlines ---
nexttile;
contourf(R_grid, Z_grid, Delta_T_grid, 30, 'LineStyle', 'none');
hold on;
colormap(gca, parula);
cb = colorbar;
ylabel(cb, '$\Delta T$ (K)', 'Interpreter', 'latex', 'FontSize', 11);
caxis([0, max(Delta_T)]);
h_stream = streamslice(R_grid, Z_grid, u_grid, w_grid, 2);
set(h_stream, 'Color', 'w', 'LineWidth', 0.5);
xline(R_illum, 'r--', 'LineWidth', 1.5);
xlabel('$r$ ($\mu$m)', 'FontSize', 12);
ylabel('$z$ ($\mu$m)', 'FontSize', 12);
title('(a) Temperature \& flow', 'FontSize', 13, 'FontWeight', 'normal');
axis equal;
xlim([0, plot_params.r_max]);
ylim([0, plot_params.z_max]);

% --- Panel (b): Velocity magnitude ---
nexttile;
contourf(R_grid, Z_grid, vel_mag_grid, 30, 'LineStyle', 'none');
hold on;
colormap(gca, viridis_like());
cb = colorbar;
ylabel(cb, '$|\mathbf{v}|$ ($\mu$m/s)', 'Interpreter', 'latex', 'FontSize', 11);
xline(R_illum, 'w--', 'LineWidth', 1.5);
xlabel('$r$ ($\mu$m)', 'FontSize', 12);
ylabel('$z$ ($\mu$m)', 'FontSize', 12);
title('(b) Velocity magnitude', 'FontSize', 13, 'FontWeight', 'normal');
axis equal;
xlim([0, plot_params.r_max]);
ylim([0, plot_params.z_max]);

% --- Panel (c): Total force vectors ---
nexttile;
quiver(R_grid(1:skip:end, 1:skip:end), Z_grid(1:skip:end, 1:skip:end), ...
       F_total_r_grid(1:skip:end, 1:skip:end), F_total_z_grid(1:skip:end, 1:skip:end), ...
       2, 'k', 'LineWidth', 0.8);
hold on;
xline(R_illum, 'r--', 'LineWidth', 1.5);
xlabel('$r$ ($\mu$m)', 'FontSize', 12);
ylabel('$z$ ($\mu$m)', 'FontSize', 12);
title('(c) Total force $\mathbf{F}_{\mathrm{total}}$', 'FontSize', 13, 'FontWeight', 'normal');
axis equal;
xlim([0, plot_params.r_max]);
ylim([0, plot_params.z_max]);
grid on;

% --- Panel (d): Radial velocity profile ---
nexttile;
hold on;
fill([r_centers, fliplr(r_centers)], [max(v_r_mean, 0), zeros(1, n_bins)], ...
     [0.9, 0.7, 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
fill([r_centers, fliplr(r_centers)], [min(v_r_mean, 0), zeros(1, n_bins)], ...
     [0.7, 0.8, 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(r_centers, v_r_mean, 'b-', 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 0.8);
xline(R_illum, 'r--', 'LineWidth', 1.5);
xlabel('$r$ ($\mu$m)', 'FontSize', 12);
ylabel('$\langle v_r \rangle$ ($\mu$m/s)', 'FontSize', 12);
title(sprintf('(d) Radial velocity ($z < %d$ $\\mu$m)', z_thresh), ...
      'FontSize', 13, 'FontWeight', 'normal');
xlim([0, plot_params.r_max]);
grid on;

if export_figures
    exportgraphics(fig5, ['Fig5_Combined_Overview.' export_format], 'Resolution', export_dpi);
    fprintf('  Saved: Fig5_Combined_Overview.%s\n', export_format);
end

%% ========================================================================
%% SUMMARY OUTPUT
%% ========================================================================

fprintf('\n========================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('========================================\n');
fprintf('File: %s\n', filename);
fprintf('R_illum: %d μm\n', R_illum);
fprintf('ΔT_max: %.2f K\n', max(Delta_T));
fprintf('v_r max (outflow): %.4f μm/s\n', max(v_r_mean));
fprintf('v_r min (inflow): %.4f μm/s\n', min(v_r_mean));
fprintf('----------------------------------------\n');
if export_figures
    fprintf('Figures exported as %s at %d dpi\n', upper(export_format), export_dpi);
end
fprintf('========================================\n');

%% ========================================================================
%% HELPER FUNCTION: VIRIDIS-LIKE COLORMAP
%% ========================================================================

function cmap = viridis_like()
    % Approximate viridis colormap
    cmap = [
        0.267004, 0.004874, 0.329415;
        0.282327, 0.140926, 0.457517;
        0.253935, 0.265254, 0.529983;
        0.206756, 0.371758, 0.553117;
        0.163625, 0.471133, 0.558148;
        0.127568, 0.566949, 0.550556;
        0.134692, 0.658636, 0.517649;
        0.266941, 0.748751, 0.440573;
        0.477504, 0.821444, 0.318195;
        0.741388, 0.873449, 0.149561;
        0.993248, 0.906157, 0.143936
    ];
    % Interpolate to 256 colors
    x = linspace(1, size(cmap, 1), 256);
    cmap = interp1(1:size(cmap, 1), cmap, x);
end