close all
clear all

filename = 'Fields_75um.csv';

% Check if file exists
if ~exist(filename, 'file')
    error('File not found! Please update the filename variable with the correct path.');
end

%% Load and Parse Data
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
u = data(:, 9);        % Velocity r-component (m/s)
w = data(:, 10);       % Velocity z-component (m/s)

fprintf('Data loaded: %d points\n\n', length(r));

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

fprintf('Interpolation complete!\n\n');

%% Extract Bottom Boundary Data (z ≈ z_min)
fprintf('==============================================\n');
fprintf('   BOTTOM SURFACE ANALYSIS\n');
fprintf('==============================================\n\n');

% Find the row closest to the bottom
[~, bottom_idx] = min(abs(Z_grid(:,1) - z_min));

% Extract bottom data
r_bottom = R_grid(bottom_idx, :);
z_bottom = Z_grid(bottom_idx, :);
T_bottom = T_grid(bottom_idx, :);
u_bottom_COMSOL = u_grid(bottom_idx, :);  % This includes slip from COMSOL

fprintf('Bottom boundary at z = %.2f μm\n', z_bottom(1)*1e6);
fprintf('Temperature range: %.2f - %.2f K\n', min(T_bottom), max(T_bottom));
fprintf('ΔT = %.2f K\n\n', max(T_bottom) - min(T_bottom));

%% Calculate Temperature Gradient at Bottom
% ∂T/∂r using gradient function
grad_T_r_bottom = gradient(T_bottom, r_grid);

fprintf('Temperature gradient ∂T/∂r:\n');
fprintf('  Maximum (most positive): %.2f K/μm\n', max(grad_T_r_bottom)*1e-6);
fprintf('  Minimum (most negative): %.2f K/μm\n', min(grad_T_r_bottom)*1e-6);
fprintf('  At center (r=0): %.2f K/μm\n\n', grad_T_r_bottom(1)*1e-6);

%% Calculate Theoretical Slip Velocity
fprintf('==============================================\n');
fprintf('   THERMO-OSMOTIC SLIP CALCULATION\n');
fprintf('==============================================\n\n');

% Thermo-osmotic coefficient from COMSOL
chi = 1.28e-10;  % m²/s

fprintf('Using χ = %.2e m²/s (from COMSOL)\n', chi);
fprintf('Converting: χ = %.2f μm²/(s·K)\n\n', chi*1e12);

% Theoretical slip velocity: v_slip = (χ/T) * (∂T/∂r)
u_slip_theory = (chi ./ T_bottom) .* grad_T_r_bottom;

fprintf('Theoretical slip velocity:\n');
fprintf('  Maximum: %.3f μm/s\n', max(abs(u_slip_theory))*1e6);
fprintf('  At r=0: %.3f μm/s\n', u_slip_theory(1)*1e6);
fprintf('  Direction: ');
if mean(u_slip_theory) > 0
    fprintf('OUTWARD (positive r)\n\n');
else
    fprintf('INWARD (negative r)\n\n');
end

fprintf('  Maximum: %.3f μm/s\n', max(abs(u_bottom_COMSOL))*1e6);
fprintf('  At r=0: %.3f μm/s\n', u_bottom_COMSOL(1)*1e6);
fprintf('  Direction: ');
if mean(u_bottom_COMSOL) > 0
    fprintf('OUTWARD (positive r)\n\n');
else
    fprintf('INWARD (negative r)\n\n');
end

% Calculate difference/agreement
difference = abs(u_slip_theory - u_bottom_COMSOL);
relative_error = mean(difference ./ (abs(u_bottom_COMSOL) + 1e-15)) * 100;

fprintf('Comparison:\n');
fprintf('  Mean absolute difference: %.3e m/s\n', mean(difference));
fprintf('  Relative difference: %.1f%%\n\n', relative_error);


%% Subplot 1: Temperature at Bottom
% nexttile
plot(r_bottom*1e6, T_bottom, 'LineWidth', 2.5, 'Color', [0.8 0.2 0.2]);
xlabel('r [μm]', 'FontSize', 12);
ylabel('Temperature [K]', 'FontSize', 12);
% title('Temperature at Bottom Surface (z ≈ 0)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);
xlim([0, max(r_bottom)*1e6]);

% Add annotations
text(0.05, 0.95, sprintf('T_{max} = %.2f K', max(T_bottom)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'white', 'EdgeColor', 'k');
text(0.05, 0.85, sprintf('ΔT = %.2f K', max(T_bottom) - min(T_bottom)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'white', 'EdgeColor', 'k');

%% Subplot 2: Radial Temperature Gradient
% nexttile
figure
plot(r_bottom*1e6, grad_T_r_bottom*1e-6, 'LineWidth', 2.5, 'Color', [0.2 0.2 0.8]);
hold on;
yline(0, '--k', 'LineWidth', 1.5);
xlabel('r [μm]', 'FontSize', 12);
ylabel('∂T/∂r [K/μm]', 'FontSize', 12);
% title('Radial Temperature Gradient at Bottom (Driving Force)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);
xlim([0, max(r_bottom)*1e6]);

