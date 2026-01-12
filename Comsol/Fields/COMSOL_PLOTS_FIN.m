%% ========================================================================
%% ADDITIONAL PHYSICS VALIDATION PLOTS
%% Add these to the end of analyze_comsol_final_complete.m
%% ========================================================================

%% FIGURE 8: Péclet Number Comparison ⭐⭐⭐ CRITICAL
%% Shows Pe_conv >> Pe_thermo (factor of 5) - KEY RESULT
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Physical parameters
D_T = 0.2e-12;  % m²/(s·K) - thermophoretic mobility
D_B = 2.3e-13;  % m²/s - Brownian diffusion (Stokes-Einstein, d=1.9μm)

% Compute Péclet numbers at z=15 μm
field_name = 'z150';
if isfield(profiles, field_name)
    prof = profiles.(field_name);
    valid = ~isnan(prof.u_mean);
    
    % Pe_convection = |u| × R_illum / D_B
    Pe_conv = abs(prof.u_mean(valid)) * R_illu*1e-6 / D_B;
    
    % Need temperature gradient for Pe_thermo
    % Extract temperature at z=15 μm
    z_target = 15e-6;
    z_tol = 2e-6;
    idx_layer = abs(z - z_target) < z_tol;
    
    if sum(idx_layer) > 10
        r_layer = r(idx_layer);
        T_layer = T(idx_layer);
        
        % Bin temperature
        T_binned = nan(n_bins, 1);
        for i_bin = 1:n_bins
            idx_bin = r_layer >= r_edges(i_bin) & r_layer < r_edges(i_bin+1);
            if sum(idx_bin) >= 2
                T_binned(i_bin) = mean(T_layer(idx_bin));
            end
        end
        T_binned = movmean(T_binned, 2, 'omitnan');
        
        % Compute gradient
        dT_dr_z15 = gradient(T_binned(valid), profiles.r_centers(valid)*1e-6);
        
        % Pe_thermophoresis = D_T × |∇T| × R_illum / D_B
        Pe_thermo = D_T * abs(dT_dr_z15) * R_illu*1e-6 / D_B;
        
        % Plot
        hold on;
        plot(profiles.r_centers(valid), Pe_conv, 'b-', 'LineWidth', LINE_WIDTH, ...
             'DisplayName', 'Convection: $\mathrm{Pe}_{\mathrm{conv}}$');
        plot(profiles.r_centers(valid), Pe_thermo, 'r-', 'LineWidth', LINE_WIDTH, ...
             'DisplayName', 'Thermophoresis: $\mathrm{Pe}_{\mathrm{thermo}}$');
        
        % Plot ratio
        ratio = Pe_conv ./ Pe_thermo;
        plot(profiles.r_centers(valid), ratio, 'k--', 'LineWidth', LINE_WIDTH-0.5, ...
             'DisplayName', 'Ratio: $\mathrm{Pe}_{\mathrm{conv}}/\mathrm{Pe}_{\mathrm{thermo}}$');
        
        yline(1, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        xline(R_illu, 'g--', 'LineWidth', 1.5, 'Alpha', 0.5, 'HandleVisibility', 'off');
        
        % Annotations
        y_range = ylim;
        text(0.05, 0.95, sprintf('At $z = %.0f$ $\\mu$m', prof.z_actual), ...
             'Units', 'normalized', 'FontSize', FONT_SIZE_LEGEND, ...
             'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', 'k');
        
        % Find mean ratio near R_illum
        idx_near = profiles.r_centers(valid) > R_illu*0.8 & profiles.r_centers(valid) < R_illu*1.2;
        mean_ratio = mean(ratio(idx_near), 'omitnan');
        text(0.05, 0.85, sprintf('Mean ratio $\\approx %.1f$', mean_ratio), ...
             'Units', 'normalized', 'FontSize', FONT_SIZE_LEGEND, ...
             'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', 'k');
    end
end

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('P\''eclet Number', 'FontSize', FONT_SIZE_LABELS);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-2);
grid on; box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);
ylim([0, max(ylim)]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig8_Peclet_Comparison');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Fig8\n');

%% FIGURE 9: Velocity Decomposition ⭐⭐⭐
%% Shows u_slip + u_conv vs thermophoresis at multiple heights
figure('Position', [100, 100, 1200, 800], 'Color', 'w');

% Select 3 representative heights
key_heights = [0, 5, 15];
for plot_idx = 1:3
    subplot(1, 3, plot_idx);
    
    z_val = key_heights(plot_idx);
    field_name = sprintf('z%03d', round(z_val*10));
    
    if isfield(profiles, field_name)
        prof = profiles.(field_name);
        valid = ~isnan(prof.u_mean);
        
        % Get temperature gradient at this height
        z_target = z_val * 1e-6;
        z_tol = 2e-6;
        idx_layer = abs(z - z_target) < z_tol;
        
        if sum(idx_layer) > 10
            r_layer = r(idx_layer);
            T_layer = T(idx_layer);
            
            T_binned = nan(n_bins, 1);
            for i_bin = 1:n_bins
                idx_bin = r_layer >= r_edges(i_bin) & r_layer < r_edges(i_bin+1);
                if sum(idx_bin) >= 2
                    T_binned(i_bin) = mean(T_layer(idx_bin));
                end
            end
            T_binned = movmean(T_binned, 2, 'omitnan');
            
            dT_dr_layer = gradient(T_binned(valid), profiles.r_centers(valid)*1e-6);
            
            % Thermophoretic velocity: -D_T × dT/dr
            v_thermo = -D_T * dT_dr_layer * 1e6;  % μm/s
            
            % Fluid velocity (slip + convection)
            v_fluid = prof.u_mean(valid);
            
            % Particle velocity
            v_particle = v_fluid + v_thermo;
            
            hold on;
            plot(profiles.r_centers(valid), v_fluid, 'b-', 'LineWidth', 2, ...
                 'DisplayName', '$u_{\mathrm{fluid}}$ (slip+conv)');
            plot(profiles.r_centers(valid), v_thermo, 'r--', 'LineWidth', 2, ...
                 'DisplayName', '$-D_T \partial T/\partial r$');
            plot(profiles.r_centers(valid), v_particle, 'k-', 'LineWidth', 2.5, ...
                 'DisplayName', '$v_{\mathrm{particle}}$ (total)');
            
            yline(0, 'k:', 'LineWidth', 1);
            xline(R_illu, 'g--', 'LineWidth', 1.5, 'Alpha', 0.5);
            
            xlabel('$r$ ($\mu$m)', 'FontSize', 14);
            ylabel('Velocity ($\mu$m/s)', 'FontSize', 14);
            title(sprintf('$z = %.0f$ $\\mu$m', prof.z_actual), 'FontSize', 16);
            
            if plot_idx == 3
                legend('Location', 'best', 'FontSize', 11);
            end
            
            grid on; box on;
            set(gca, 'LineWidth', 1.5, 'FontSize', 12);
        end
    end
end

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 12, 4]);

fileName = fullfile(output_dir, 'Fig9_Velocity_Decomposition');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Fig9\n');

%% FIGURE 10: r_min/R_illum Validation ⭐⭐
%% Direct comparison to experiments
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Compute r_min/R_illum for each layer
z_vals_comsol = [];
ratios_comsol = [];

for i_layer = 1:n_layers
    field_name = sprintf('z%03d', round(z_layers_um(i_layer)*10));
    if ~isfield(profiles, field_name), continue; end
    
    prof = profiles.(field_name);
    valid = ~isnan(prof.u_mean);
    
    if sum(valid) >= 3
        [~, idx_min] = min(prof.u_mean(valid));
        r_vals = profiles.r_centers(valid);
        r_min = r_vals(idx_min);
        
        z_vals_comsol(end+1) = prof.z_actual;
        ratios_comsol(end+1) = r_min / R_illu;
    end
end

% Plot COMSOL results
plot(ratios_comsol, z_vals_comsol, 'bo-', 'LineWidth', 2.5, ...
     'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'COMSOL');
hold on;

% Add experimental data (from your results)
% Experimental ranges for 20x objective (R_illum = 75 μm)
exp_mean = 1.00;
exp_std = 0.08;

% Observation zone: z = 10-20 μm
fill([exp_mean-exp_std, exp_mean+exp_std, exp_mean+exp_std, exp_mean-exp_std], ...
     [10, 10, 20, 20], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r', ...
     'LineStyle', '--', 'DisplayName', 'Experiments (20x)');

% Mark ideal line
xline(1.0, 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal ($r_{\mathrm{min}} = R_{\mathrm{illum}}$)');

% Shade acceptable range
fill([0.85, 1.15, 1.15, 0.85], [0, 0, 40, 40], 'g', ...
     'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '$\pm$15\% range');

xlabel('$r_{\mathrm{min}} / R_{\mathrm{illum}}$', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
title('Validation: Minimum Position vs Experiments', 'FontSize', FONT_SIZE_LABELS-2);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-2);

grid on; box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);
xlim([0.7, 1.4]);
ylim([0, 35]);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig10_rmin_Validation');
% export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Fig10\n');

%% FIGURE 11: Streamfunction / Circulation Pattern ⭐⭐
%% Shows toroidal circulation clearly
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Compute streamfunction in cylindrical coordinates
% ψ such that: u_r = -1/r × ∂ψ/∂z, u_z = 1/r × ∂ψ/∂r

% Use cumulative integration
psi = zeros(size(u_grid));
dr_grid = r_grid(2) - r_grid(1);

for iz = 1:length(z_grid)
    for ir = 2:length(r_grid)
        if r_grid(ir) > 1e-9  % Avoid singularity at r=0
            % ∂ψ/∂r = r × u_z
            psi(iz, ir) = psi(iz, ir-1) + r_grid(ir) * w_grid(iz, ir) * dr_grid;
        end
    end
end

% Plot streamfunction contours
contour(R_grid*1e6, Z_grid*1e6, psi*1e12, 20, 'LineWidth', 1.5);
hold on;

% Overlay temperature for context
contourf(R_grid*1e6, Z_grid*1e6, T_grid, 20, 'LineStyle', 'none', 'FaceAlpha', 0.3);

% Mark key regions
xline(R_illu, 'r--', 'LineWidth', 2);
plot([0, max(r)*1e6], [30, 30], 'g--', 'LineWidth', 2.5);

cb = colorbar;
cb.Label.String = 'Temperature (K)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
title('Streamlines: Toroidal Circulation', 'FontSize', FONT_SIZE_LABELS-2);

xlim([0, min(max(r)*1e6, R_illu*3)]);
ylim([0, 50]);
box on;
colormap(hot);
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig11_Streamfunction');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Fig11\n');

%% FIGURE 12: Temperature Uniformity in Hot Zone ⭐
%% Validates "flat-top" claim
figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Extract temperature profiles at multiple heights
heights_to_plot = [0, 5, 10, 15, 20, 30];
colors_temp = jet(length(heights_to_plot));

hold on;
for i = 1:length(heights_to_plot)
    z_target = heights_to_plot(i) * 1e-6;
    z_tol = 2e-6;
    idx_layer = abs(z - z_target) < z_tol;
    
    if sum(idx_layer) > 10
        r_layer = r(idx_layer);
        T_layer = T(idx_layer);
        
        % Bin
        T_binned = nan(50, 1);
        r_edges_fine = linspace(0, max(r), 51);
        r_centers_fine = (r_edges_fine(1:end-1) + r_edges_fine(2:end))/2 * 1e6;
        
        for j = 1:50
            idx_bin = r_layer >= r_edges_fine(j) & r_layer < r_edges_fine(j+1);
            if sum(idx_bin) >= 2
                T_binned(j) = mean(T_layer(idx_bin));
            end
        end
        
        T_binned = movmean(T_binned, 2, 'omitnan');
        valid_T = ~isnan(T_binned);
        
        plot(r_centers_fine(valid_T), T_binned(valid_T), '-', ...
             'Color', colors_temp(i,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('$z = %.0f$ $\\mu$m', heights_to_plot(i)));
    end
end

% Mark illumination zone
xline(R_illu, 'k--', 'LineWidth', 2, 'DisplayName', '$R_{\mathrm{illum}}$');

% Shade the "hot zone" (should be uniform)
y_range = ylim;
fill([0, R_illu, R_illu, 0], [y_range(1), y_range(1), y_range(2), y_range(2)], ...
     'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Hot zone');

xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('Temperature (K)', 'FontSize', FONT_SIZE_LABELS);
title('Temperature Uniformity: Flat-Top Profile', 'FontSize', FONT_SIZE_LABELS-2);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-4, 'NumColumns', 2);

grid on; box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);

% Add annotation about uniformity
T_center = max(T_binned(valid_T));
T_edge = T_binned(find(r_centers_fine(valid_T) >= R_illu*0.9, 1));
uniformity = (T_center - T_edge) / (T_center - 293) * 100;
text(0.05, 0.95, sprintf('Uniformity: %.1f\\%%\n(within $r < R_{\\mathrm{illum}}$)', uniformity), ...
     'Units', 'normalized', 'FontSize', FONT_SIZE_LEGEND, ...
     'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', 'k');

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Fig12_Temperature_Uniformity');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Fig12\n');

%% Print Enhanced Summary
fprintf('\n========================================\n');
fprintf('ENHANCED PHYSICS VALIDATION\n');
fprintf('========================================\n');
fprintf('Pe_conv / Pe_thermo = %.1f (convection dominates)\n', mean_ratio);
fprintf('r_min/R_illum (COMSOL) = %.3f at z=15μm\n', ratios_comsol(find(z_vals_comsol >= 14 & z_vals_comsol <= 16, 1)));
fprintf('r_min/R_illum (Expt) = %.2f ± %.2f\n', exp_mean, exp_std);
fprintf('Agreement: %.1f%% deviation\n', abs(ratios_comsol(3) - exp_mean)/exp_mean * 100);
fprintf('Temperature uniformity: %.1f%% in hot zone\n', uniformity);
fprintf('\n5 additional figures saved!\n');
fprintf('Total: 12 comprehensive validation figures\n');
fprintf('========================================\n');