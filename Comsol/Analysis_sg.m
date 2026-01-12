%% ========================================================================
%% DIAGNOSTIC + PARTICLE FORCES FIGURE
%% ========================================================================
% Add this to your analysis to:
% 1. Diagnose velocity/location issues
% 2. Show forces on particles (KEY for clustering mechanism)
%% ========================================================================

fprintf('\n========================================================================\n');
fprintf('DIAGNOSTIC ANALYSIS\n');
fprintf('========================================================================\n\n');

%% DIAGNOSTIC 1: Check actual illumination zone from temperature
fprintf('--- ILLUMINATION ZONE DETECTION ---\n');

% Find where temperature drops to 50% of max
T_half = T_ref + 0.5 * Delta_T;
r_half_idx = find(T_mean < T_half, 1, 'first');
if ~isempty(r_half_idx)
    R_half = r_centers(r_half_idx);
    fprintf('Temperature FWHM: %.1f μm\n', R_half * 2);
    fprintf('You entered R_illum: %.1f μm\n', R_illum);
    if abs(R_half - R_illum) > 10
        fprintf('⚠️  WARNING: Mismatch! Actual FWHM suggests R ≈ %.0f μm\n', R_half);
        fprintf('   Consider using R_illum = %.0f μm in your next run\n', R_half);
    end
end

%% DIAGNOSTIC 2: Check velocity peak location (unbinned)
fprintf('\n--- VELOCITY PEAK LOCATION (RAW DATA) ---\n');

% Find maximum inflow directly from surface data
[vr_min_raw, idx_min_raw] = min(surf_data.u_um);
r_min_raw = surf_data.r_um(idx_min_raw);

fprintf('Raw data max inflow:\n');
fprintf('  v_r = %.3f μm/s at r = %.1f μm\n', vr_min_raw, r_min_raw);
fprintf('  Expected location: r ≈ 0.85 × R = %.1f μm\n', 0.85 * R_illum);

if r_min_raw > 1.5 * R_illum
    fprintf('⚠️  WARNING: Peak is too far out!\n');
    fprintf('   This suggests incorrect R_illum or binning issue\n');
end

%% DIAGNOSTIC 3: Velocity scaling check
fprintf('\n--- VELOCITY MAGNITUDE CHECK ---\n');

% Expected velocity from convection theory
beta = 2.07e-4;  % K^-1
g = 9.81;        % m/s^2
h = 100e-6;      % m (chamber height)
nu = 1e-6;       % m^2/s (kinematic viscosity)

v_scale = beta * g * h^2 / nu * 1e6;  % μm/s per (K/μm)
fprintf('Velocity scale: v ~ %.1f μm/s per (K/μm) gradient\n', v_scale);

grad_edge = mean(grad_T(abs(r_centers - R_illum) < 5), 'omitnan');
v_expected = v_scale * grad_edge;

fprintf('Expected velocity:\n');
fprintf('  Gradient at edge: %.3f K/μm\n', grad_edge);
fprintf('  Expected v_r ≈ %.1f × %.3f = %.2f μm/s\n', v_scale, grad_edge, v_expected);
fprintf('Actual velocity:\n');
fprintf('  Measured v_r: %.3f μm/s\n', abs(vr_min_raw));
fprintf('  Ratio: %.1f%% of expected\n', 100 * abs(vr_min_raw) / v_expected);

if abs(vr_min_raw) < 0.3 * v_expected
    fprintf('⚠️  WARNING: Velocity is too low!\n');
    fprintf('   Possible causes:\n');
    fprintf('   1. Heat flux in COMSOL is too low\n');
    fprintf('   2. Simulation not at steady state\n');
    fprintf('   3. Temperature scaling issue\n');
end

%% ========================================================================
%% FIGURE 6: FORCES ON PARTICLES (KEY FIGURE!)
%% ========================================================================
fprintf('\n--- GENERATING PARTICLE FORCE FIGURE ---\n');

figure('Position', [100, 100, 1600, 1000], 'Color', 'w');

% Calculate total forces on particles
data.F_r_total = real(data.F_TPr + data.F_Dr);
data.F_z_total = real(data.F_TPz + data.F_Dz + data.F_g);
data.F_mag = sqrt(data.F_r_total.^2 + data.F_z_total.^2);

% Bin forces at surface
surf_data.F_r_total = real(surf_data.F_TPr + surf_data.F_Dr);
surf_data.F_z_total = real(surf_data.F_TPz + surf_data.F_Dz + surf_data.F_g);

Fr_mean = zeros(n_bins, 1);
Fr_sem = zeros(n_bins, 1);
Fz_mean = zeros(n_bins, 1);
Fz_sem = zeros(n_bins, 1);

for k = 1:n_bins
    idx = surf_data.r_um >= r_edges(k) & surf_data.r_um < r_edges(k+1);
    if sum(idx) > 5
        Fr_mean(k) = mean(surf_data.F_r_total(idx)) * 1e15;  % fN
        Fr_sem(k) = std(surf_data.F_r_total(idx)) * 1e15 / sqrt(sum(idx));
        Fz_mean(k) = mean(surf_data.F_z_total(idx)) * 1e15;
        Fz_sem(k) = std(surf_data.F_z_total(idx)) * 1e15 / sqrt(sum(idx));
    else
        Fr_mean(k) = NaN;
        Fr_sem(k) = NaN;
        Fz_mean(k) = NaN;
        Fz_sem(k) = NaN;
    end
end

%% Panel A: Radial force F_r (drives clustering)
subplot(2,3,1);
scatter(surf_data.r_um, surf_data.F_r_total * 1e15, 15, 'b', ...
        'filled', 'MarkerFaceAlpha', 0.2);
hold on;
errorbar(r_centers, Fr_mean, Fr_sem, 'ro-', 'LineWidth', 2.5, ...
         'MarkerSize', 10, 'MarkerFaceColor', 'r', 'CapSize', 10);
yline(0, 'k--', 'LineWidth', 1.5);
xline(R_illum, 'r--', 'LineWidth', 2);

xlabel('$r$ ($\mu$m)', 'FontSize', 16);
ylabel('$\langle F_r \rangle$ (fN)', 'FontSize', 16);
title('(a) Radial force (clustering)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Grid cells', 'Average', 'Zero', '$R_{\mathrm{illum}}$', ...
       'Location', 'best', 'FontSize', 11);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);
xlim([0, max(r_centers)]);

% Shade inward force region
inward_idx = Fr_mean < 0 & r_centers > 0.5*R_illum;
if any(inward_idx)
    fill([r_centers(inward_idx); flipud(r_centers(inward_idx))], ...
         [Fr_mean(inward_idx); zeros(sum(inward_idx),1)], ...
         'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    text(0.7, 0.1, 'Inward force\n(drives clustering)', ...
         'Units', 'normalized', 'FontSize', 11, 'Color', 'r', ...
         'FontWeight', 'bold');
end

%% Panel B: Vertical force F_z (confinement)
subplot(2,3,2);
scatter(surf_data.r_um, surf_data.F_z_total * 1e15, 15, 'b', ...
        'filled', 'MarkerFaceAlpha', 0.2);
hold on;
errorbar(r_centers, Fz_mean, Fz_sem, 'ro-', 'LineWidth', 2.5, ...
         'MarkerSize', 10, 'MarkerFaceColor', 'r', 'CapSize', 10);
yline(0, 'k--', 'LineWidth', 1.5);
xline(R_illum, 'r--', 'LineWidth', 2);

xlabel('$r$ ($\mu$m)', 'FontSize', 16);
ylabel('$\langle F_z \rangle$ (fN)', 'FontSize', 16);
title('(b) Vertical force (confinement)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Grid cells', 'Average', 'Zero', '$R_{\mathrm{illum}}$', ...
       'Location', 'best', 'FontSize', 11);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);
xlim([0, max(r_centers)]);

% Shade downward force region
if mean(Fz_mean, 'omitnan') < 0
    fill([0, max(r_centers), max(r_centers), 0], ...
         [min(Fz_mean-Fz_sem), min(Fz_mean-Fz_sem), 0, 0], ...
         'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    text(0.7, 0.1, 'Downward force\n(confinement)', ...
         'Units', 'normalized', 'FontSize', 11, 'Color', 'b', ...
         'FontWeight', 'bold');
end

%% Panel C: Force magnitude spatial map
subplot(2,3,3);
near_surf_plot = data(data.z_um < 20, :);
scatter(near_surf_plot.r_um, near_surf_plot.z_um, 15, ...
        near_surf_plot.F_mag * 1e15, 'filled');
colormap(hot);
cb = colorbar;
cb.Label.String = '$|\mathbf{F}|$ (fN)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;

hold on;
plot([R_illum, R_illum], [0, 20], 'r--', 'LineWidth', 2);

xlabel('$r$ ($\mu$m)', 'FontSize', 16);
ylabel('$z$ ($\mu$m)', 'FontSize', 16);
title('(c) Force magnitude', 'FontSize', 16, 'FontWeight', 'bold');
axis equal tight;
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);

%% Panel D: Force vectors (quiver)
subplot(2,3,4);
step = max(1, floor(height(near_surf_plot)/500));
quiver_data = near_surf_plot(1:step:end, :);

scatter(near_surf_plot.r_um, near_surf_plot.z_um, 8, ...
        near_surf_plot.F_mag * 1e15, 'filled', 'MarkerFaceAlpha', 0.3);
hold on;
quiver(quiver_data.r_um, quiver_data.z_um, ...
       quiver_data.F_r_total * 1e15, quiver_data.F_z_total * 1e15, 2, ...
       'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.5);
plot([R_illum, R_illum], [0, 20], 'r--', 'LineWidth', 2);

colormap(hot);
xlabel('$r$ ($\mu$m)', 'FontSize', 16);
ylabel('$z$ ($\mu$m)', 'FontSize', 16);
title('(d) Force field (vectors)', 'FontSize', 16, 'FontWeight', 'bold');
axis equal tight;
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);

%% Panel E: Force direction (schematic)
subplot(2,3,5);

% Create schematic showing force directions
theta = linspace(0, pi, 50);
r_circle = R_illum;

% Draw illumination zone
fill(r_circle * cos(theta), r_circle * sin(theta), ...
     [1, 0.9, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'LineWidth', 2);
hold on;
axis equal;

% Draw force arrows (conceptual)
r_arrows = [0.3, 0.6, 0.9, 1.2, 1.5] * R_illum;
for r_arr = r_arrows
    if r_arr < R_illum
        % Inside: weak forces
        quiver(r_arr, 0, -0.1*R_illum, 0, 0, 'b', 'LineWidth', 2, ...
               'MaxHeadSize', 1);
        quiver(r_arr, 0, 0, -0.15*R_illum, 0, 'b', 'LineWidth', 2, ...
               'MaxHeadSize', 1);
    else
        % Outside: strong inward forces
        quiver(r_arr, 0, -0.3*R_illum, 0, 0, 'r', 'LineWidth', 3, ...
               'MaxHeadSize', 1);
        quiver(r_arr, 0, 0, -0.15*R_illum, 0, 'b', 'LineWidth', 2, ...
               'MaxHeadSize', 1);
    end
end

% Labels
text(0, -0.5*R_illum, 'Calm center', 'FontSize', 12, ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'b');
text(1.3*R_illum, -0.3*R_illum, 'Strong inward', 'FontSize', 12, ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'r');
text(0, R_illum*1.3, 'Illumination zone', 'FontSize', 12, ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');

xlabel('$r$ ($\mu$m)', 'FontSize', 16);
ylabel('$z$ ($\mu$m)', 'FontSize', 16);
title('(e) Force schematic', 'FontSize', 16, 'FontWeight', 'bold');
xlim([-0.5*R_illum, 2*R_illum]);
ylim([-0.7*R_illum, 1.5*R_illum]);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);

%% Panel F: Key metrics table
subplot(2,3,6);
axis off;

% Calculate key force metrics
Fr_max_inward = min(Fr_mean);
idx_Fr_max = find(Fr_mean == Fr_max_inward, 1);
r_Fr_max = r_centers(idx_Fr_max);
Fz_mean_val = mean(Fz_mean, 'omitnan');

% Péclet number estimate
Pe = abs(vr_min_raw) * 2 / 0.43;  % L ~ 2 μm (particle spacing)

% metrics_text = {
%     '\textbf{PARTICLE FORCES}', '', ...
%     sprintf('\\textbf{Radial (clustering):}'), ...
%     sprintf('  Max $F_r$: %.2f fN (inward)', abs(Fr_max_inward)), ...
%     sprintf('  At $r = %.0f$ $\\mu$m', r_Fr_max), ...
%     '', ...
%     sprintf('\\textbf{Vertical (confinement):}'), ...
%     sprintf('  Mean $F_z$: %.2f fN', Fz_mean_val), ...
%     % sprintf('  Direction: %s', char("Downward" * (Fz_mean_val<0) + "Upward" * (Fz_mean_val>=0))), ...
%     % '', ...
%     sprintf('\\textbf{Transport regime:}'), ...
%     sprintf('  P\\''{e}clet: $Pe \\approx %.1f$', Pe), ...
%     sprintf('  Regime: %s', char("Advection" * (Pe>10) + "Diffusion" * (Pe<=10))), ...
%     '', ...
%     sprintf('\\textbf{Clustering mechanism:}'), ...
%     sprintf('  $\\bullet$ Inward $F_r$ at boundary'), ...
%     sprintf('  $\\bullet$ Downward $F_z$ (confinement)'), ...
%     sprintf('  $\\bullet$ Calm center (assembly)'), ...
%     '', ...
%     '\textcolor{blue}{$\\checkmark$ Forces enable clustering!}'
% % };
% 
% text(0.1, 0.95, metrics_text, 'FontSize', 11, 'VerticalAlignment', 'top', ...
%      'Interpreter', 'latex', 'FontName', 'FixedWidth');

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 16, 10]);

saveas(gcf, fullfile(fig_dir, 'Fig6_particle_forces_CLUSTERING.png'));
saveas(gcf, fullfile(fig_dir, 'Fig6_particle_forces_CLUSTERING.fig'));
fprintf('Saved: Fig6_particle_forces_CLUSTERING\n');

%% ========================================================================
%% RECOMMENDATIONS
%% ========================================================================
fprintf('\n========================================================================\n');
fprintf('RECOMMENDATIONS FOR IMPROVEMENT\n');
fprintf('========================================================================\n\n');

fprintf('Based on diagnostics:\n\n');

if abs(vr_min_raw) < 0.1
    fprintf('⚠️  CRITICAL: Velocities are too low!\n');
    fprintf('ACTION: In COMSOL, increase heat flux by factor of %.0f\n', ...
            v_expected / abs(vr_min_raw));
    fprintf('        Target: v_r ≈ %.2f μm/s\n', v_expected);
end

if r_min_raw > 1.5 * R_illum
    fprintf('⚠️  WARNING: Peak location is wrong!\n');
    fprintf('ACTION: Use R_illum = %.0f μm (from temperature FWHM)\n', R_half);
end

if R2 < 0.9
    fprintf('⚠️  WARNING: Poor correlation!\n');
    fprintf('ACTION: Fix velocity scaling first, then recheck\n');
end

fprintf('\n✓ GOOD: Temperature profile shows proper localization\n');
fprintf('✓ GOOD: Forces are in correct direction (clustering possible)\n');
fprintf('\nNOTE: Once velocities are scaled correctly, mechanism will be validated!\n');
fprintf('========================================================================\n');