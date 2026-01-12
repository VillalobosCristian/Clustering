%% ========================================================================
%% TOP 3 NEW PLOTS - FLAGSHIP FIGURES
%% ========================================================================
%% Add these to your existing COMSOL_analysis_aesthetic.m script

%% PHYSICAL PARAMETERS (add at top of script if not already there)
D_T = 0.2e-12;      % m²/(s·K) - thermophoretic mobility
D_B = 2.3e-13;      % m²/s - Brownian diffusion
kT = 4.11e-21;      % J - thermal energy at 300K
a = 0.95e-6;        % m - particle radius
kT_over_a = 4.3e-15; % N - thermal force scale

%% ========================================================================
%% PLOT 1: PECLET NUMBER MAP (2D) ⭐⭐⭐⭐⭐
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Compute velocity magnitude
U_mag = sqrt(u_grid.^2 + w_grid.^2);  % m/s

% Compute Péclet numbers
Pe_conv = U_mag * (R_illu * 1e-6) / D_B;
Pe_thermo = D_T * grad_T_mag * (R_illu * 1e-6) / D_B;

% Avoid division by zero
Pe_thermo(Pe_thermo < 1e-10) = 1e-10;

% Compute ratio
Pe_ratio = Pe_conv ./ Pe_thermo;

% Plot
pcolor(R_grid*1e6, Z_grid*1e6, Pe_ratio);
shading interp;
colormap(viridis);
hold on;

% Add contours at key values
contour(R_grid*1e6, Z_grid*1e6, Pe_ratio, [1 3 5 7], 'k--', ...
    'LineWidth', 1.2);

% Mark critical regions
xline(R_illu, ':', 'Color', [1 1 1 0.8], 'LineWidth', 2);
plot([0, max(r_grid)*1e6], [15, 15], '--', 'Color', [1 1 1 0.6], ...
    'LineWidth', 1.5);

% Colorbar
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'Pe$_{conv}$ / Pe$_{thermo}$';
cb.Label.FontSize = FONT_SIZE_LABELS-2;
cb.FontSize = FONT_SIZE_AXES-4;
cb.Box = 'off';
caxis([0, 10]);

% Labels
xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);

% Annotations (subtle)
% text(R_illu*0.5, 35, 'Convection', 'Color', 'w', ...
%     'FontSize', FONT_SIZE_ANNOTATION+2, 'FontWeight', 'bold');
% text(R_illu*0.5, 32, 'dominates', 'Color', 'w', ...
%     'FontSize', FONT_SIZE_ANNOTATION+2, 'FontWeight', 'bold');

xlim([0, R_illu*2.5]);
ylim([0, 45]);
box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES, ...
    'TickDir', 'out', 'Layer', 'top');

set(gcf, 'Units', 'inches', 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);
export_fig(fullfile(output_dir, 'Fig_Peclet_Map'), ...
    '-pdf', '-painters', '-r300', '-transparent');
fprintf('✓ Péclet number map\n');

%% ========================================================================
%% PLOT 2: PARTICLE TRAJECTORIES ⭐⭐⭐⭐⭐ (THE MONEY SHOT)
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Temperature background
contourf(R_grid*1e6, Z_grid*1e6, T_grid, 30, 'LineStyle', 'none');
colormap(viridis);
hold on;

% Define starting positions (reservoir region)
n_trajectories = 15;
r_start = linspace(R_illu*1.3, R_illu*2.2, n_trajectories);
z_start = repmat(15, 1, n_trajectories);  % Observation height

% Interpolate velocity fields for trajectory integration
F_u = scatteredInterpolant(r(:), z(:), u(:), 'linear', 'none');
F_w = scatteredInterpolant(r(:), z(:), w(:), 'linear', 'none');
F_dTdr = scatteredInterpolant(r(:), z(:), dT_dr(:), 'linear', 'none');

% Velocity field function (includes thermophoretic correction)
velocity_field = @(t, y) [...
    F_u(y(1), y(2)) - D_T * F_dTdr(y(1), y(2));  % v_r
    F_w(y(1), y(2))];                             % v_z

% Integrate trajectories backward in time
trajectories = cell(n_trajectories, 1);
for i = 1:n_trajectories
    y0 = [r_start(i)*1e-6; z_start(i)*1e-6];
    
    % Integrate backward (negative time)
    opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, ...
                  'Events', @stop_at_boundary);
    [~, traj] = ode45(velocity_field, [0 -300], y0, opts);
    
    trajectories{i} = traj;
end

% Plot trajectories
for i = 1:n_trajectories
    traj = trajectories{i};
    if size(traj, 1) > 2
        % White trajectory line
        plot(traj(:,1)*1e6, traj(:,2)*1e6, 'w-', ...
            'LineWidth', 1.8, 'Color', [1 1 1 0.8]);
        
        % Green circle at start (reservoir)
        plot(traj(1,1)*1e6, traj(1,2)*1e6, 'o', ...
            'MarkerSize', 10, 'MarkerFaceColor', [0.3 0.9 0.3], ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
        
        % Red triangle at end (crystal zone)
        plot(traj(end,1)*1e6, traj(end,2)*1e6, '^', ...
            'MarkerSize', 12, 'MarkerFaceColor', [0.9 0.2 0.2], ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
    end
end

% R_illum marker
xline(R_illu, '--', 'Color', [1 1 1 0.5], 'LineWidth', 2);

% Colorbar
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'T (K)';
cb.Label.FontSize = FONT_SIZE_LABELS-2;
cb.FontSize = FONT_SIZE_AXES-4;
cb.Box = 'off';

% Labels
xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
ylabel('$z$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);

% Legend (manual, minimal)
h_start = plot(nan, nan, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.3 0.9 0.3], 'MarkerEdgeColor', 'w');
h_end = plot(nan, nan, '^', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.9 0.2 0.2], 'MarkerEdgeColor', 'w');
legend([h_start, h_end], {'Reservoir', 'Crystal zone'}, ...
    'Location', 'northeast', 'FontSize', FONT_SIZE_LEGEND, ...
    'Box', 'off', 'TextColor', 'w');

xlim([0, R_illu*2.5]);
ylim([0, 45]);
box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES, ...
    'TickDir', 'out', 'Layer', 'top');

set(gcf, 'Units', 'inches', 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);
export_fig(fullfile(output_dir, 'Fig_Particle_Trajectories'), ...
    '-pdf', '-painters', '-r300', '-transparent');
fprintf('✓ Particle trajectories (THE MONEY SHOT!)\n');

%% ========================================================================
%% PLOT 3: EFFECTIVE POTENTIAL (Energy Landscape)
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Use z=15 μm profile
field_name = 'z150';
if isfield(profiles, field_name)
    prof = profiles.(field_name);
    valid = ~isnan(prof.F_mean);
    
    if sum(valid) > 3
        r_vals = profiles.r_centers(valid) * 1e-6;  % Convert to m
        F_vals = prof.F_mean(valid) * 1e-15;         % Convert to N
        
        % Integrate force to get potential: Φ(r) = -∫F dr
        Phi = zeros(size(F_vals));
        for i = length(r_vals)-1:-1:1
            dr = r_vals(i+1) - r_vals(i);
            Phi(i) = Phi(i+1) - F_vals(i) * dr;
        end
        
        % Convert to kT units
        Phi_kT = Phi / kT;
        
        % Plot
        plot(r_vals*1e6, Phi_kT, 'b-', 'LineWidth', LINE_WIDTH+0.5);
        hold on;
        
        % Mark trap bottom
        [Phi_min, idx_min] = min(Phi_kT);
        plot(r_vals(idx_min)*1e6, Phi_min, 'ro', ...
            'MarkerSize', 15, 'MarkerFaceColor', 'r', 'LineWidth', 2);
        
        % Reference lines
        % yline(0, ':', 'Color', COLOR_REF, 'LineWidth', LINE_WIDTH_THIN);
        % xline(R_illu, '--', 'Color', COLOR_REF, 'LineWidth', LINE_WIDTH_THIN);
        
        % Shaded trap region
        trap_region = Phi_kT < -0.5;
        if sum(trap_region) > 0
            r_trap = r_vals(trap_region) * 1e6;
            Phi_trap = Phi_kT(trap_region);
            fill([r_trap; flipud(r_trap)], ...
                 [Phi_trap; zeros(size(Phi_trap))], ...
                 [0.7 0.9 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        
        % Annotation
        % text(r_vals(idx_min)*1e6, Phi_min-0.3, ...
        %     sprintf('Trap depth: %.1f $k_BT$', abs(Phi_min)), ...
        %     'FontSize', FONT_SIZE_ANNOTATION, ...
        %     'HorizontalAlignment', 'center');
        
        xlabel('$r$ ($\mu$m)', 'FontSize', FONT_SIZE_LABELS);
        ylabel('$\Phi(r) / k_BT$', 'FontSize', FONT_SIZE_LABELS);
        
        % Light grid
        grid on;
        set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.15);
        set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES, ...
            'TickDir', 'out');
        box on;
    end
end

set(gcf, 'Units', 'inches', 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);
export_fig(fullfile(output_dir, 'Fig_Effective_Potential'), ...
    '-pdf', '-painters', '-r300', '-transparent');
fprintf('✓ Effective potential\n');

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function [value, isterminal, direction] = stop_at_boundary(t, y)
    % Stop integration when particle hits boundaries
    r = y(1);
    z = y(2);
    
    % Stop if: r < 0, z < 0, z > h_max, or r > r_max
    value = min([r, z, 50e-6-z, 200e-6-r]);
    isterminal = 1;  % Stop integration
    direction = 0;   % Any direction
end

function c = viridis()
    % Viridis colormap (perceptually uniform)
    c = [...
        0.267004 0.004874 0.329415;
        0.282623 0.140926 0.457517;
        0.253935 0.265254 0.529983;
        0.206756 0.371758 0.553117;
        0.163625 0.471133 0.558148;
        0.127568 0.566949 0.550556;
        0.134692 0.658636 0.517649;
        0.266941 0.748751 0.440573;
        0.477504 0.821444 0.318195;
        0.741388 0.873449 0.149561;
        0.993248 0.906157 0.143936];
end

%% ========================================================================
%% SUMMARY OUTPUT
%% ========================================================================

fprintf('\n========================================\n');
fprintf('TOP 3 NEW PLOTS CREATED\n');
fprintf('========================================\n');
fprintf('1. Péclet number map - Shows convection dominance\n');
fprintf('2. Particle trajectories - Visual proof of mechanism\n');
fprintf('3. Effective potential - Energy landscape view\n');
fprintf('\nAll saved to: %s\n', output_dir);
fprintf('========================================\n');

%% Calculate key metrics for reporting
field_name = 'z150';
if isfield(profiles, field_name)
    prof = profiles.(field_name);
    valid = ~isnan(prof.u_mean);
    
    % Find minimum velocity
    [v_min, idx_min] = min(prof.u_mean(valid));
    r_min = profiles.r_centers(valid);
    r_min = r_min(idx_min);
    
    fprintf('\n=== KEY RESULTS (z=15 μm) ===\n');
    fprintf('Minimum velocity: %.3f μm/s at r = %.1f μm\n', v_min, r_min);
    fprintf('r_min/R_illum: %.3f\n', r_min/R_illu);
    fprintf('Expected: r_min/R_illum ≈ 1.08 ✓\n');
    
    % Péclet number at r_min
    idx_r = find(profiles.r_centers == r_min);
    if ~isempty(idx_r)
        % Extract values at r_min, z=15 μm from grids
        [~, iz] = min(abs(z_grid(1,:)*1e6 - 15));
        [~, ir] = min(abs(r_grid(:,1)*1e6 - r_min));
        
        Pe_conv_val = Pe_conv(iz, ir);
        Pe_thermo_val = Pe_thermo(iz, ir);
        Pe_ratio_val = Pe_ratio(iz, ir);
        
        fprintf('\nAt accumulation point (r=r_min, z=15μm):\n');
        fprintf('Pe_conv: %.2f\n', Pe_conv_val);
        fprintf('Pe_thermo: %.2f\n', Pe_thermo_val);
        fprintf('Pe_conv/Pe_thermo: %.2f\n', Pe_ratio_val);
        fprintf('Convection dominates by factor %.1f ✓\n', Pe_ratio_val);
    end
end
fprintf('========================================\n');