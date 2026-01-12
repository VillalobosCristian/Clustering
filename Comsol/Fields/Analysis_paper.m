clear all; clc; close all;

%% ========================================================================
%% COMSOL ANALYSIS - ESSENTIAL FIGURES FOR PAPER
%% ========================================================================
% Generates only the key publication-quality figures
% ~7 figures total, matching experimental analysis style
%% ========================================================================

set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');

%% LOAD DATA
fprintf('========================================================================\n');
fprintf('COMSOL ANALYSIS - PAPER FIGURES\n');
fprintf('========================================================================\n\n');

[file_name, path_name] = uigetfile('*.csv', 'Select COMSOL Fields CSV');
if isequal(file_name, 0), error('No file selected'); end

[~, sim_name, ~] = fileparts(file_name);
fprintf('Analyzing: %s\n', sim_name);

% Parameters
prompt = {'Profile type:', 'Width w (μm):', 'Shape d:', 'Illumination radius R (μm):'};
defaultans = {'SuperGaussian', '22.4', '5.5', '35'};
answer = inputdlg(prompt, 'Parameters', [1 50], defaultans);
if isempty(answer), error('Parameters not provided'); end

profile_type = answer{1};
w_um = str2double(answer{2});
d_shape = str2double(answer{3});
R_illum = str2double(answer{4});

fprintf('Profile: %s (w=%.1f μm, d=%.1f)\n', profile_type, w_um, d_shape);
fprintf('R_illum: %.0f μm\n\n', R_illum);

% Load CSV
data = readtable(fullfile(path_name, file_name), 'HeaderLines', 8);
data.Properties.VariableNames = {'r', 'z', 'T', 'F_TPr', 'F_Dr', ...
                                  'F_TPz', 'F_Dz', 'F_g', 'u', 'w'};

% Convert to μm and μm/s
data.r_um = data.r * 1e6;
data.z_um = data.z * 1e6;
data.u_um = real(data.u * 1e6);  % Use real part
data.w_um = real(data.w * 1e6);

fprintf('Loaded %d points\n', height(data));
fprintf('Domain: r=[%.0f,%.0f] μm, z=[%.0f,%.0f] μm\n\n', ...
        min(data.r_um), max(data.r_um), min(data.z_um), max(data.z_um));

% Output directory
fig_dir = fullfile(pwd, 'Paper_Figures', sim_name);
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Extract surface data (z < 5 μm)
T_ref = 293.15;
Delta_T = max(data.T) - T_ref;
surf_data = data(data.z_um < 5, :);
surf_data = sortrows(surf_data, 'r_um');

fprintf('ΔT = %.2f K\n', Delta_T);
fprintf('Surface data: %d points\n\n', height(surf_data));

% Binning
n_bins = 50;
r_edges = linspace(0, max(surf_data.r_um), n_bins+1);
r_centers = (r_edges(1:end-1) + r_edges(2:end))/2;

%% ========================================================================
%% FIGURE 1: Temperature Profile & Gradient (2 panels)
%% ========================================================================
fprintf('Generating Figure 1: Temperature profiles...\n');

figure('Position', [100, 100, 1400, 600], 'Color', 'w');

% Bin temperature
T_mean = zeros(n_bins, 1);
for k = 1:n_bins
    idx = surf_data.r_um >= r_edges(k) & surf_data.r_um < r_edges(k+1);
    if sum(idx) > 0
        T_mean(k) = mean(surf_data.T(idx));
    else
        T_mean(k) = NaN;
    end
end

% Calculate gradient
grad_T = abs(gradient(T_mean, mean(diff(r_centers))));

%% Panel A: Temperature
subplot(1,2,1);
plot(r_centers, T_mean - T_ref, 'b-', 'LineWidth', 3);
hold on;
xline(R_illum, 'r--', 'LineWidth', 2, 'DisplayName', '$R_{\mathrm{illum}}$');

xlabel('$r$ ($\mu$m)', 'FontSize', 18);
ylabel('$\Delta T$ (K)', 'FontSize', 18);
title('(a) Temperature profile', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 14);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 16);
xlim([0, max(r_centers)]);

%% Panel B: Gradient
subplot(1,2,2);
plot(r_centers, grad_T, 'r-', 'LineWidth', 3);
hold on;
xline(R_illum, 'r--', 'LineWidth', 2);

% Shade "flat" region
flat_threshold = 0.05;  % K/μm
flat_idx = grad_T < flat_threshold & r_centers < R_illum;
if any(flat_idx)
    fill([r_centers(flat_idx), fliplr(r_centers(flat_idx))], ...
         [zeros(sum(flat_idx),1)', fliplr(grad_T(flat_idx)')], ...
         'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Flat region');
end

xlabel('$r$ ($\mu$m)', 'FontSize', 18);
ylabel('$|\nabla T|$ (K/$\mu$m)', 'FontSize', 18);
title('(b) Temperature gradient', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 16);
xlim([0, max(r_centers)]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 14, 6]);
saveas(gcf, fullfile(fig_dir, 'Fig1_temperature_gradient.png'));
saveas(gcf, fullfile(fig_dir, 'Fig1_temperature_gradient.fig'));

%% ========================================================================
%% FIGURE 2: Velocity Profiles (2 panels: v_r and v_z)
%% ========================================================================
fprintf('Generating Figure 2: Velocity profiles...\n');

figure('Position', [100, 100, 1400, 600], 'Color', 'w');

% Bin velocities
vr_mean = zeros(n_bins, 1);
vr_sem = zeros(n_bins, 1);
vz_mean = zeros(n_bins, 1);
vz_sem = zeros(n_bins, 1);

for k = 1:n_bins
    idx = surf_data.r_um >= r_edges(k) & surf_data.r_um < r_edges(k+1);
    if sum(idx) > 5
        vr_mean(k) = mean(surf_data.u_um(idx));
        vr_sem(k) = std(surf_data.u_um(idx)) / sqrt(sum(idx));
        vz_mean(k) = mean(surf_data.w_um(idx));
        vz_sem(k) = std(surf_data.w_um(idx)) / sqrt(sum(idx));
    else
        vr_mean(k) = NaN;
        vr_sem(k) = NaN;
        vz_mean(k) = NaN;
        vz_sem(k) = NaN;
    end
end

%% Panel A: Radial velocity
subplot(1,2,1);
scatter(surf_data.r_um, surf_data.u_um, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.2);
hold on;
errorbar(r_centers, vr_mean, vr_sem, 'ro-', 'LineWidth', 2.5, ...
         'MarkerSize', 10, 'MarkerFaceColor', 'r', 'CapSize', 10);
yline(0, 'k--', 'LineWidth', 1.5);
xline(R_illum, 'r--', 'LineWidth', 2);

[vr_min, idx_min] = min(vr_mean);
r_min = r_centers(idx_min);
plot(r_min, vr_min, 'gs', 'MarkerSize', 15, 'LineWidth', 3);

xlabel('$r$ ($\mu$m)', 'FontSize', 18);
ylabel('$\langle v_r \rangle$ ($\mu$m/s)', 'FontSize', 18);
title(sprintf('(a) Radial velocity (min: %.2f $\\mu$m/s)', vr_min), ...
      'FontSize', 18, 'FontWeight', 'bold');
legend('Grid cells', 'Binned', 'Zero', '$R_{\mathrm{illum}}$', 'Max inflow', ...
       'Location', 'southwest', 'FontSize', 12);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 16);
xlim([0, max(r_centers)]);

%% Panel B: Vertical velocity
subplot(1,2,2);
scatter(surf_data.r_um, surf_data.w_um, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.2);
hold on;
errorbar(r_centers, vz_mean, vz_sem, 'ro-', 'LineWidth', 2.5, ...
         'MarkerSize', 10, 'MarkerFaceColor', 'r', 'CapSize', 10);
yline(0, 'k--', 'LineWidth', 1.5);
xline(R_illum, 'r--', 'LineWidth', 2);

xlabel('$r$ ($\mu$m)', 'FontSize', 18);
ylabel('$\langle v_z \rangle$ ($\mu$m/s)', 'FontSize', 18);
title(sprintf('(b) Vertical velocity (max: %.2f $\\mu$m/s)', max(vz_mean)), ...
      'FontSize', 18, 'FontWeight', 'bold');
legend('Grid cells', 'Binned', 'Zero', '$R_{\mathrm{illum}}$', ...
       'Location', 'best', 'FontSize', 12);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 16);
xlim([0, max(r_centers)]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 14, 6]);
saveas(gcf, fullfile(fig_dir, 'Fig2_velocity_profiles.png'));
saveas(gcf, fullfile(fig_dir, 'Fig2_velocity_profiles.fig'));

%% ========================================================================
%% FIGURE 3: Streamlines (Flow Pattern)
%% ========================================================================
fprintf('Generating Figure 3: Streamlines...\n');

figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Create grid
r_grid = linspace(0, max(data.r_um), 100);
z_grid = linspace(0, max(data.z_um), 100);
[R_grid, Z_grid] = meshgrid(r_grid, z_grid);

% Interpolate velocities
F_u = scatteredInterpolant(data.r_um, data.z_um, data.u_um, 'linear', 'none');
F_w = scatteredInterpolant(data.r_um, data.z_um, data.w_um, 'linear', 'none');
U_grid = F_u(R_grid, Z_grid);
W_grid = F_w(R_grid, Z_grid);
speed_grid = sqrt(U_grid.^2 + W_grid.^2);

% Plot
contourf(R_grid, Z_grid, speed_grid, 20, 'LineStyle', 'none', 'FaceAlpha', 0.6);
hold on;
streamslice(R_grid, Z_grid, U_grid, W_grid, 2);

% Mark illumination zone
plot([R_illum, R_illum], [0, max(z_grid)], 'r--', 'LineWidth', 3);

colormap(parula);
cb = colorbar;
cb.Label.String = 'Speed ($\mu$m/s)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.FontSize = 14;

xlabel('$r$ ($\mu$m)', 'FontSize', 18);
ylabel('$z$ ($\mu$m)', 'FontSize', 18);
title('Velocity streamlines', 'FontSize', 18, 'FontWeight', 'bold');

axis equal tight;
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 16);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);
saveas(gcf, fullfile(fig_dir, 'Fig3_streamlines.png'));
saveas(gcf, fullfile(fig_dir, 'Fig3_streamlines.fig'));

%% ========================================================================
%% FIGURE 4: Force Analysis - F_z vs Height (Center)
%% ========================================================================
fprintf('Generating Figure 4: Force analysis...\n');

figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Extract center data
center_data = data(data.r_um < 5, :);
center_data = sortrows(center_data, 'z_um');

% Use real parts
F_TPz = real(center_data.F_TPz) * 1e15;
F_Dz = real(center_data.F_Dz) * 1e15;
F_g = real(center_data.F_g(1)) * 1e15;
F_z_net = F_TPz + F_Dz + F_g;

plot(F_TPz, center_data.z_um, 'o-', 'LineWidth', 2.5, 'MarkerSize', 6, ...
     'DisplayName', 'Thermophoresis');
hold on;
plot(F_Dz, center_data.z_um, 's-', 'LineWidth', 2.5, 'MarkerSize', 6, ...
     'DisplayName', 'Drag');
plot(F_g * ones(size(center_data.z_um)), center_data.z_um, '^-', ...
     'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Gravity');
plot(F_z_net, center_data.z_um, 'k-', 'LineWidth', 3, ...
     'DisplayName', 'Net $F_z$');

xline(0, 'k--', 'LineWidth', 1.5);
yline(10, 'r:', 'LineWidth', 2, 'Alpha', 0.5);

xlabel('$F_z$ (fN)', 'FontSize', 18);
ylabel('$z$ ($\mu$m)', 'FontSize', 18);
title('Vertical forces at center ($r < 5$ $\mu$m)', 'FontSize', 18, 'FontWeight', 'bold');

legend('Location', 'southeast', 'FontSize', 14);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 16);
ylim([0, min(50, max(center_data.z_um))]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 8, 6]);
saveas(gcf, fullfile(fig_dir, 'Fig4_forces_vertical.png'));
saveas(gcf, fullfile(fig_dir, 'Fig4_forces_vertical.fig'));

%% ========================================================================
%% FIGURE 5: Mechanism Validation - ∇T vs v_r Correlation
%% ========================================================================
fprintf('Generating Figure 5: Mechanism validation...\n');

figure('Position', [100, 100, 1400, 600], 'Color', 'w');

%% Panel A: Normalized profiles overlay
subplot(1,2,1);
yyaxis left
plot(r_centers, grad_T / max(grad_T), 'r-', 'LineWidth', 3);
ylabel('$|\nabla T| / |\nabla T|_{\max}$', 'FontSize', 16);
set(gca, 'YColor', 'r');

yyaxis right
plot(r_centers, abs(vr_mean) / max(abs(vr_mean)), 'b-', 'LineWidth', 3);
ylabel('$|v_r| / |v_r|_{\max}$', 'FontSize', 16);
set(gca, 'YColor', 'b');

xlabel('$r$ ($\mu$m)', 'FontSize', 16);
title('(a) Profiles: $v_r \propto \nabla T$', 'FontSize', 16, 'FontWeight', 'bold');
legend('$|\nabla T|$', '$|v_r|$', 'Location', 'best', 'FontSize', 12);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);
xlim([0, max(r_centers)]);

%% Panel B: Direct correlation
subplot(1,2,2);

valid_idx = ~isnan(grad_T) & ~isnan(vr_mean);
grad_T_valid = grad_T(valid_idx);
vr_valid = abs(vr_mean(valid_idx));

scatter(grad_T_valid, vr_valid, 80, r_centers(valid_idx), 'filled', ...
        'MarkerEdgeColor', 'k');
hold on;

% Linear fit
p = polyfit(grad_T_valid, vr_valid, 1);
grad_fit = linspace(min(grad_T_valid), max(grad_T_valid), 100);
vr_fit = polyval(p, grad_fit);
plot(grad_fit, vr_fit, 'k--', 'LineWidth', 2);

% R²
vr_pred = polyval(p, grad_T_valid);
R2 = 1 - sum((vr_valid - vr_pred).^2) / sum((vr_valid - mean(vr_valid)).^2);

xlabel('$|\nabla T|$ (K/$\mu$m)', 'FontSize', 16);
ylabel('$|v_r|$ ($\mu$m/s)', 'FontSize', 16);
title('(b) Linear correlation', 'FontSize', 16, 'FontWeight', 'bold');

text(0.05, 0.95, sprintf('$v_r = %.2f \\times |\\nabla T| + %.3f$\n$R^2 = %.3f$', ...
     p(1), p(2), R2), 'Units', 'normalized', 'FontSize', 12, ...
     'BackgroundColor', 'w', 'EdgeColor', 'k', 'VerticalAlignment', 'top', ...
     'Interpreter', 'latex');

colormap(parula);
cb = colorbar;
cb.Label.String = '$r$ ($\mu$m)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;

grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 14, 6]);
saveas(gcf, fullfile(fig_dir, 'Fig5_correlation_mechanism.png'));
saveas(gcf, fullfile(fig_dir, 'Fig5_correlation_mechanism.fig'));

%% ========================================================================
%% SUMMARY & SAVE DATA
%% ========================================================================

fprintf('\n========================================================================\n');
fprintf('SUMMARY FOR PAPER\n');
fprintf('========================================================================\n\n');

fprintf('SIMULATION: %s\n', sim_name);
fprintf('Profile: %s (w=%.1f μm, d=%.1f, R=%.0f μm)\n\n', ...
        profile_type, w_um, d_shape, R_illum);

fprintf('TEMPERATURE:\n');
fprintf('  ΔT = %.2f K\n', Delta_T);
fprintf('  Max |∇T| = %.3f K/μm\n', max(grad_T));
fprintf('  Center |∇T| (r<5μm) = %.4f K/μm\n', mean(grad_T(r_centers<5)));
fprintf('  Edge |∇T| (r≈R) = %.3f K/μm\n', ...
        mean(grad_T(abs(r_centers-R_illum)<5)));

fprintf('\nVELOCITY:\n');
fprintf('  Max inflow: %.3f μm/s at r = %.1f μm\n', vr_min, r_min);
fprintf('  Max upwelling: %.3f μm/s\n', max(vz_mean));
fprintf('  Center velocity (r<5μm): %.3f μm/s\n', mean(abs(vr_mean(r_centers<5))));

fprintf('\nFORCES (near surface z<10μm):\n');
near_surf = data(data.z_um < 10, :);
F_z_mean = mean(real(near_surf.F_TPz + near_surf.F_Dz + near_surf.F_g)) * 1e15;
fprintf('  Mean F_z = %.2f fN\n', F_z_mean);
if F_z_mean < 0
    fprintf('  Direction: DOWNWARD ✓\n');
else
    fprintf('  Direction: UPWARD\n');
end

fprintf('\nMECHANISM:\n');
fprintf('  v_r vs ∇T slope: %.2f μm/s per (K/μm)\n', p(1));
fprintf('  R² = %.3f\n', R2);
fprintf('  Expected (theory): ~20 μm/s per (K/μm)\n');
if abs(p(1)-20) < 10
    fprintf('  Agreement: EXCELLENT ✓\n');
else
    fprintf('  Agreement: Check scaling\n');
end

% Save summary
summary = struct('sim_name', sim_name, 'profile', profile_type, ...
                 'w_um', w_um, 'd', d_shape, 'R_illum', R_illum, ...
                 'Delta_T', Delta_T, 'grad_max', max(grad_T), ...
                 'vr_min', vr_min, 'vz_max', max(vz_mean), ...
                 'slope', p(1), 'R2', R2, 'F_z_mean', F_z_mean);

save(fullfile(fig_dir, 'summary.mat'), 'summary', 'r_centers', ...
     'T_mean', 'grad_T', 'vr_mean', 'vr_sem', 'vz_mean', 'vz_sem');

fprintf('\n========================================================================\n');
fprintf('COMPLETE! Generated 5 figures in:\n%s\n', fig_dir);
fprintf('========================================================================\n');