close all
clear all

filename = 'Fields_25um.csv';

% Check if file exists
if ~exist(filename, 'file')
    error('File not found! Please update the filename variable with the correct path.');
end

%% Load and Parse Data
fprintf('==============================================\n');
fprintf('   COMSOL DATA ANALYSIS\n');
fprintf('   Bottom Surface & Force Analysis\n');
fprintf('==============================================\n\n');
fprintf('Loading data...\n');

% Read the header information
fid = fopen(filename, 'r');
header_lines = cell(9, 1);
for i = 1:9
    header_lines{i} = fgetl(fid);
end
fclose(fid);

% Load numerical data (skip first 9 header lines)
data = readmatrix(filename, 'NumHeaderLines', 9);

% Extract columns
r = data(:, 1);        % r coordinate (m)
z = data(:, 2);        % z coordinate (m)
T = data(:, 3);        % Temperature (K)
F_TPr = data(:, 4);    % Thermophoretic force r-component (N)
F_Dr = data(:, 5);     % Drag force r-component (N)
F_TPz = data(:, 6);    % Thermophoretic force z-component (N)
F_Dz = data(:, 7);     % Drag force z-component (N)
F_g = data(:, 8);      % Gravitational force (N)
u = data(:, 9);        % Velocity r-component (m/s)
w = data(:, 10);       % Velocity z-component (m/s)

% Handle complex numbers (take real part)
F_TPr = real(F_TPr);
F_TPz = real(F_TPz);
F_Dr = real(F_Dr);
F_Dz = real(F_Dz);
F_g = real(F_g)*1e-3;

% Calculate total forces
F_total_r = F_TPr + F_Dr;
F_total_z = F_TPz + F_Dz + F_g;
F_total_mag = sqrt(F_total_r.^2 + F_total_z.^2);

fprintf('Data loaded: %d points\n', length(r));
fprintf('File info: %s\n\n', header_lines{2});

%% Create Interpolation Grid
fprintf('Creating interpolation grid...\n');

% Determine grid boundaries
r_min = min(r); r_max = max(r);
z_min = min(z); z_max = max(z);

% Create regular grid for interpolation
grid_resolution = 200;
r_grid = linspace(r_min, r_max, grid_resolution);
z_grid = linspace(z_min, z_max, grid_resolution);
[R_grid, Z_grid] = meshgrid(r_grid, z_grid);

% Interpolate data onto regular grid
T_grid = griddata(r, z, T, R_grid, Z_grid, 'natural');
u_grid = griddata(r, z, u, R_grid, Z_grid, 'natural');
F_TPr_grid = griddata(r, z, F_TPr, R_grid, Z_grid, 'natural');
F_Dr_grid = griddata(r, z, F_Dr, R_grid, Z_grid, 'natural');
F_TPz_grid = griddata(r, z, F_TPz, R_grid, Z_grid, 'natural');
F_Dz_grid = griddata(r, z, F_Dz, R_grid, Z_grid, 'natural');
F_total_r_grid = griddata(r, z, F_total_r, R_grid, Z_grid, 'natural');
F_total_z_grid = griddata(r, z, F_total_z, R_grid, Z_grid, 'natural');

fprintf('Interpolation complete!\n\n');

%% ========================================================================
%% PART 1: BOTTOM SURFACE ANALYSIS (Temperature & Slip)
%% ========================================================================

fprintf('==============================================\n');
fprintf('   BOTTOM SURFACE ANALYSIS\n');
fprintf('==============================================\n\n');

% Find the row closest to the bottom
[~, bottom_idx] = min(abs(Z_grid(:,1) - z_min));

% Extract bottom data
r_bottom = R_grid(bottom_idx, :);
z_bottom = Z_grid(bottom_idx, :);
T_bottom = T_grid(bottom_idx, :);
u_bottom_COMSOL = u_grid(bottom_idx, :);  % Includes slip from COMSOL
F_TPr_bottom = F_TPr_grid(bottom_idx, :);
F_Dr_bottom = F_Dr_grid(bottom_idx, :);
F_TPz_bottom = F_TPz_grid(bottom_idx, :);
F_Dz_bottom = F_Dz_grid(bottom_idx, :);
F_total_r_bottom = F_total_r_grid(bottom_idx, :);
F_total_z_bottom = F_total_z_grid(bottom_idx, :);

fprintf('Bottom boundary location: z = %.3f μm\n', z_bottom(1)*1e6);
fprintf('\nTemperature statistics:\n');
fprintf('  T_max = %.2f K (at center)\n', max(T_bottom));
fprintf('  T_min = %.2f K (at edge)\n', min(T_bottom));
fprintf('  ΔT = %.2f K\n', max(T_bottom) - min(T_bottom));

%% Calculate Temperature Gradient at Bottom
% ∂T/∂r using gradient function
grad_T_r_bottom = gradient(T_bottom, r_grid);

fprintf('\nTemperature gradient ∂T/∂r:\n');
fprintf('  At center (r=0): %.2e K/m = %.2f K/μm\n', ...
    grad_T_r_bottom(1), grad_T_r_bottom(1)*1e-6);
fprintf('  Maximum magnitude: %.2e K/m = %.2f K/μm\n', ...
    max(abs(grad_T_r_bottom)), max(abs(grad_T_r_bottom))*1e-6);
fprintf('  Most negative: %.2e K/m = %.2f K/μm\n', ...
    min(grad_T_r_bottom), min(grad_T_r_bottom)*1e-6);

%% Calculate Theoretical Slip Velocity
fprintf('\n==============================================\n');
fprintf('   THERMO-OSMOTIC SLIP ANALYSIS\n');
fprintf('==============================================\n\n');

% Thermo-osmotic coefficient from COMSOL
chi = -1.28e-10;  % m²/s

fprintf('Parameters:\n');
fprintf('  χ = %.2e m²/s (from COMSOL setup)\n', chi);
fprintf('  Converting: χ = %.2f μm²/(s·K)\n', chi*1e12);
fprintf('  Note: This is equivalent to K^(d) ~ %.1f μm²/(s·K) at T=300K\n\n', chi*300*1e12);

% Theoretical slip velocity: v_slip = (χ/T) * (∂T/∂r)
u_slip_theory = (chi ./ T_bottom) .* grad_T_r_bottom;

fprintf('Slip velocity analysis:\n');
fprintf('  Theory formula: u_slip = (χ/T) × (∂T/∂r)\n\n');

fprintf('Theoretical slip:\n');
fprintf('  At r=0: %.3f μm/s\n', u_slip_theory(1)*1e6);
fprintf('  Maximum magnitude: %.3f μm/s\n', max(abs(u_slip_theory))*1e6);
fprintf('  Mean value: %.3f μm/s\n', mean(u_slip_theory)*1e6);

fprintf('\nCOMSOL slip (u at z=0):\n');
fprintf('  At r=0: %.3f μm/s\n', u_bottom_COMSOL(1)*1e6);
fprintf('  Maximum magnitude: %.3f μm/s\n', max(abs(u_bottom_COMSOL))*1e6);
fprintf('  Mean value: %.3f μm/s\n', mean(u_bottom_COMSOL)*1e6);

% Determine flow direction
if mean(u_slip_theory) > 0
    fprintf('\n>>> SLIP DIRECTION: OUTWARD (positive r) <<<\n');
    fprintf('    Slip flows from HOT → COLD (away from center)\n');
else
    fprintf('\n>>> SLIP DIRECTION: INWARD (negative r) <<<\n');
    fprintf('    Slip flows from COLD → HOT (toward center)\n');
end

% Calculate agreement
difference = abs(u_slip_theory - u_bottom_COMSOL);
relative_error = mean(difference ./ (abs(u_bottom_COMSOL) + 1e-20)) * 100;

fprintf('\nAgreement check:\n');
fprintf('  Mean absolute difference: %.3e m/s\n', mean(difference));
fprintf('  Relative difference: %.1f%%\n', relative_error);
if relative_error < 10
    fprintf('  ✓ Excellent agreement (< 10%%)\n');
elseif relative_error < 50
    fprintf('  ⚠ Moderate agreement (%.1f%%)\n', relative_error);
else
    fprintf('  ✗ Poor agreement (%.1f%%)\n', relative_error);
end

%% ========================================================================
%% PART 2: FORCE ANALYSIS AT BOTTOM
%% ========================================================================

fprintf('\n==============================================\n');
fprintf('   FORCE ON PARTICLES (at bottom)\n');
fprintf('==============================================\n\n');

F_total_mag_bottom = sqrt(F_total_r_bottom.^2 + F_total_z_bottom.^2);

fprintf('Force components at r=0 (center):\n');
fprintf('  F_TP,r:    %+.3f fN\n', F_TPr_bottom(1)*1e15);
fprintf('  F_Drag,r:  %+.3f fN\n', F_Dr_bottom(1)*1e15);
fprintf('  F_total,r: %+.3f fN ', F_total_r_bottom(1)*1e15);
if F_total_r_bottom(1) < 0
    fprintf('(INWARD ← toward center)\n');
else
    fprintf('(OUTWARD → away from center)\n');
end

fprintf('\n  F_TP,z:    %+.3f fN\n', F_TPz_bottom(1)*1e15);
fprintf('  F_Drag,z:  %+.3f fN\n', F_Dz_bottom(1)*1e15);
fprintf('  F_gravity: %+.3f fN\n', F_g(1)*1e15);
fprintf('  F_total,z: %+.3f fN ', F_total_z_bottom(1)*1e15);
if F_total_z_bottom(1) < 0
    fprintf('(DOWNWARD ↓ toward surface)\n');
else
    fprintf('(UPWARD ↑ away from surface)\n');
end

fprintf('\n  |F_total|:  %.3f fN\n\n', F_total_mag_bottom(1)*1e15);

% Overall force direction assessment
fprintf('PARTICLE MOTION PREDICTION:\n');
if mean(F_total_r_bottom) < 0
    fprintf('  ✓ Radial force is INWARD - particles PULLED to center\n');
    fprintf('    This creates particle accumulation at hot spot!\n');
else
    fprintf('  ✗ Radial force is OUTWARD - particles PUSHED away\n');
    fprintf('    This prevents accumulation at center\n');
end

if mean(F_total_z_bottom) > 0
    fprintf('  • Vertical force is UPWARD - particles levitate\n');
    fprintf('    Particles hover above surface\n\n');
else
    fprintf('  • Vertical force is DOWNWARD - particles settle\n');
    fprintf('    Particles stay near bottom surface\n\n');
end

% Find equilibrium positions (where F_r ≈ 0)
F_threshold = 0.1e-15;  % 0.1 fN threshold
equilibrium_indices = find(abs(F_total_r_bottom) < F_threshold);
if ~isempty(equilibrium_indices)
    fprintf('Radial equilibrium positions (|F_r| < 0.1 fN):\n');
    for i = 1:min(3, length(equilibrium_indices))
        idx = equilibrium_indices(i);
        fprintf('  r = %.1f μm, F_r = %+.3f fN\n', ...
            r_bottom(idx)*1e6, F_total_r_bottom(idx)*1e15);
    end
    fprintf('  → Particles accumulate near these radii\n\n');
else
    fprintf('No equilibrium positions found (F_r never crosses zero)\n\n');
end

%% ========================================================================
%% FIGURE: COMPREHENSIVE BOTTOM SURFACE ANALYSIS
%% ========================================================================

fprintf('Creating comprehensive figure...\n\n');

figure('Name', 'Bottom Surface Analysis: Temperature, Slip & Forces', ...
    'Position', [50, 50, 1500, 1100]);
set(gcf, 'Color', 'white');
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Subplot 1: Temperature at Bottom
nexttile
plot(r_bottom*1e6, T_bottom, 'LineWidth', 2.5, 'Color', [0.8 0.2 0.2]);
xlabel('r [μm]', 'FontSize', 12);
ylabel('Temperature [K]', 'FontSize', 12);
title('1. Temperature at Bottom Surface (z ≈ 0)', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);
xlim([0, max(r_bottom)*1e6]);

% Add annotations
text(0.02, 0.95, sprintf('T_{max} = %.2f K', max(T_bottom)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'white', 'EdgeColor', 'k');
text(0.02, 0.85, sprintf('ΔT = %.2f K', max(T_bottom) - min(T_bottom)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'white', 'EdgeColor', 'k');

%% Subplot 2: Radial Temperature Gradient
nexttile
plot(r_bottom*1e6, grad_T_r_bottom*1e-6, 'LineWidth', 2.5, 'Color', [0.2 0.2 0.8]);
hold on;
yline(0, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('r [μm]', 'FontSize', 12);
ylabel('∂T/∂r [K/μm]', 'FontSize', 12);
title('2. Radial Temperature Gradient (Driving Force for Slip)', ...
    'FontSize', 13, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);
xlim([0, max(r_bottom)*1e6]);

% Add annotations
text(0.02, 0.05, sprintf('(∂T/∂r)_{min} = %.2f K/μm', min(grad_T_r_bottom)*1e-6), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'bottom', ...
    'BackgroundColor', 'white', 'EdgeColor', 'k');
text(0.02, 0.15, 'Negative gradient → T decreases outward', ...
    'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'bottom', ...
    'BackgroundColor', [1 1 0.8]);

%% Subplot 3: Slip Velocity Comparison
nexttile
hold on;

% Plot theoretical slip
plot(r_bottom*1e6, u_slip_theory*1e6, 'LineWidth', 2.5, 'Color', [0.8 0.2 0.2], ...
    'DisplayName', 'Theory: v_{slip} = (χ/T)(∂T/∂r)');

% Plot COMSOL slip
plot(r_bottom*1e6, u_bottom_COMSOL*1e6, '--', 'LineWidth', 2.5, 'Color', [0.2 0.6 0.2], ...
    'DisplayName', 'COMSOL: u(r, z=0)');

yline(0, ':k', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('r [μm]', 'FontSize', 12);
ylabel('Slip Velocity [μm/s]', 'FontSize', 12);
title('3. Thermo-Osmotic Slip Velocity', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 11);
xlim([0, max(r_bottom)*1e6]);

% Add annotations
if mean(u_slip_theory) > 0
    direction_text = 'SLIP DIRECTION: OUTWARD → (Hot to Cold)';
    text_color = [0.8 0.2 0.2];
else
    direction_text = 'SLIP DIRECTION: INWARD ← (Cold to Hot)';
    text_color = [0.2 0.2 0.8];
end
text(0.98, 0.95, direction_text, ...
    'Units', 'normalized', 'FontSize', 11, 'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'right', 'BackgroundColor', 'white', ...
    'EdgeColor', 'k', 'FontWeight', 'bold', 'Color', text_color);
text(0.98, 0.85, sprintf('χ = %.2f μm²/(s·K)', chi*1e12), ...
    'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'right', 'BackgroundColor', 'white', 'EdgeColor', 'k');
text(0.98, 0.75, sprintf('Agreement: %.1f%%', 100-relative_error), ...
    'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'right', 'BackgroundColor', 'white', 'EdgeColor', 'k');

% %% Subplot 4: Total Force on Particles (Single Axis)
% % nexttile
% close all
% hold on;
% 
% % Plot radial force (negative = toward center = good for accumulation!)
% plot(r_bottom*1e6, F_total_r_bottom*1e15, 'LineWidth', 3, 'Color', [0.8 0.2 0.2], ...
%     'DisplayName', 'F_{total,r} (radial)');
% 
% % Plot vertical force
% plot(r_bottom*1e6, F_total_z_bottom*1e15, 'LineWidth', 3, 'Color', [0.2 0.2 0.8], ...
%     'DisplayName', 'F_{total,z} (vertical)');
% 
% % Zero line
% yline(0, '--k', 'LineWidth', 2, 'HandleVisibility', 'off');
% 
% xlabel('r [μm]', 'FontSize', 12);
% ylabel('Force [fN]', 'FontSize', 12);
% title('4. Total Force on Particles at Bottom (TP + Drag + Gravity)', ...
%     'FontSize', 13, 'FontWeight', 'bold');
% legend('Location', 'best', 'FontSize', 11);
% grid on;
% set(gca, 'FontSize', 11);
% xlim([0, max(r_bottom)*1e6]);
% 
% % Add critical information box
% if mean(F_total_r_bottom) < 0
%     radial_status = '✓ F_r < 0: INWARD ← (ATTRACTS to center)';
%     radial_color = [0.2 0.8 0.2];  % Green for good
%     box_color = [0.8 1 0.8];  % Light green background
% else
%     radial_status = '✗ F_r > 0: OUTWARD → (REPELS from center)';
%     radial_color = [0.8 0.2 0.2];  % Red for bad
%     box_color = [1 0.8 0.8];  % Light red background
% end
% 
% if mean(F_total_z_bottom) > 0
%     vertical_status = 'F_z > 0: UPWARD ↑ (levitates)';
% else
%     vertical_status = 'F_z < 0: DOWNWARD ↓ (sediments)';
% end
% 
% % Create information box
% annotation('textbox', [0.15, 0.12, 0.25, 0.08], ...
%     'String', {radial_status, vertical_status}, ...
%     'FontSize', 11, 'FontWeight', 'bold', ...
%     'Color', radial_color, ...
%     'BackgroundColor', box_color, ...
%     'EdgeColor', 'k', 'LineWidth', 2, ...
%     'FitBoxToText', 'on');
% 
% % Add interpretation note
% text(0.98, 0.95, 'NEGATIVE F_r = toward center = ACCUMULATION', ...
%     'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top', ...
%     'HorizontalAlignment', 'right', 'BackgroundColor', [1 1 0.8], ...
%     'EdgeColor', 'k', 'FontWeight', 'bold');
% 
