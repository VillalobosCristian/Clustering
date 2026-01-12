close all
filename = 'Fields_75um.csv';

% Check if file exists
if ~exist(filename, 'file')
    error('File not found! Please update the filename variable with the correct path.');
end

%% Load and Parse Data
% fprintf('Loading data...\n');

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

% Handle complex numbers (take real part for forces)
F_TPr = real(F_TPr);
F_TPz = real(F_TPz);
F_Dr = real(F_Dr);
F_Dz = real(F_Dz);
F_g = real(F_g);

% Calculate magnitudes
F_TP_mag = sqrt(F_TPr.^2 + F_TPz.^2);  % Thermophoretic force magnitude
F_D_mag = sqrt(F_Dr.^2 + F_Dz.^2);     % Drag force magnitude
vel_mag = sqrt(u.^2 + w.^2);           % Velocity magnitude

% Calculate TOTAL FORCE (sum of all forces)
F_total_r = F_TPr + F_Dr;               % Total r-component (no r-component of gravity)
F_total_z = F_TPz + F_Dz + F_g;         % Total z-component (including gravity)
F_total_mag = sqrt(F_total_r.^2 + F_total_z.^2);  % Total force magnitude


%% Create Interpolation Grid

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
w_grid = griddata(r, z, w, R_grid, Z_grid, 'natural');
F_TPr_grid = griddata(r, z, F_TPr, R_grid, Z_grid, 'natural');
F_TPz_grid = griddata(r, z, F_TPz, R_grid, Z_grid, 'natural');
F_Dr_grid = griddata(r, z, F_Dr, R_grid, Z_grid, 'natural');
F_Dz_grid = griddata(r, z, F_Dz, R_grid, Z_grid, 'natural');
F_g_grid = griddata(r, z, F_g, R_grid, Z_grid, 'natural');
vel_mag_grid = griddata(r, z, vel_mag, R_grid, Z_grid, 'natural');
F_TP_mag_grid = griddata(r, z, F_TP_mag, R_grid, Z_grid, 'natural');
F_D_mag_grid = griddata(r, z, F_D_mag, R_grid, Z_grid, 'natural');

% Interpolate TOTAL FORCE
F_total_r_grid = griddata(r, z, F_total_r, R_grid, Z_grid, 'natural');
F_total_z_grid = griddata(r, z, F_total_z, R_grid, Z_grid, 'natural');
F_total_mag_grid = griddata(r, z, F_total_mag, R_grid, Z_grid, 'natural');


%% Figure 1: Temperature Field with Velocity Vectors

figure('Name', 'Temperature & Velocity Field', 'Position', [50, 50, 900, 750]);
set(gcf, 'Color', 'white');

skip_vec = 10;
contourf(R_grid*1e6, Z_grid*1e6, T_grid, 50, 'LineStyle', 'none');
colormap(gca, jet);  % Use jet colormap (or try: hot, parula, turbo)
cb = colorbar;
ylabel(cb, 'Temperature [K]', 'FontSize', 11);
hold on

% Add velocity vectors
quiver(R_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       Z_grid(1:skip_vec:end, 1:skip_vec:end)*1e6, ...
       u_grid(1:skip_vec:end, 1:skip_vec:end), ...
       w_grid(1:skip_vec:end, 1:skip_vec:end), 2, 'w', 'LineWidth', 1.2);

xlabel('r [\mum]', 'FontSize', 12);
ylabel('z [\mum]', 'FontSize', 12);
%title('Temperature Field with Velocity Vectors', 'FontSize', 14, 'FontWeight', 'bold');
axis equal tight;
grid on;
set(gca, 'FontSize', 11);
fig_dir='pwd'
fileName='v_t_clfield'
savefigures(fig_dir, fileName);

%% Figure 2: Individual Force Fields

figure('Name', 'Force Fields Analysis', 'Position', [100, 100, 1600, 700]);
set(gcf, 'Color', 'white');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

skip_force = 10;

% Subplot 1: Thermophoretic Force
nexttile
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_TPr_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_TPz_grid(1:skip_force:end, 1:skip_force:end)*1e15, 2, 'r', 'LineWidth', 1.2);
xlabel('r [\mum]', 'FontSize', 12);
ylabel('z [\mum]', 'FontSize', 12);
%title('Thermophoretic Force Vectors [fN]', 'FontSize', 13, 'FontWeight', 'bold');
axis equal tight;
grid on;
set(gca, 'FontSize', 11);

% Subplot 2: Drag Force
nexttile
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_Dr_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_Dz_grid(1:skip_force:end, 1:skip_force:end)*1e15, 2, 'b', 'LineWidth', 1.2);
xlabel('r [\mum]', 'FontSize', 12);
ylabel('z [\mum]', 'FontSize', 12);
%title('Drag Force Vectors [fN]', 'FontSize', 13, 'FontWeight', 'bold');
axis equal tight;
grid on;
set(gca, 'FontSize', 11);
fig_dir='pwd'
fileName='vfield'
savefigures(fig_dir, fileName);


%% Figure 3: Total Force Field

figure('Name', 'Total Force Analysis', 'Position', [150, 150, 1600, 700]);
set(gcf, 'Color', 'white');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Total Force Vectors
nexttile
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_total_r_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_total_z_grid(1:skip_force:end, 1:skip_force:end)*1e15, 2, 'k', 'LineWidth', 1.2);
xlabel('r [\mum]', 'FontSize', 12);
ylabel('z [\mum]', 'FontSize', 12);
%title('Total Force Vectors (TP + Drag + Gravity) [fN]', 'FontSize', 13, 'FontWeight', 'bold');
axis equal tight;
grid on;
set(gca, 'FontSize', 11);

% Subplot 2: Total Force Magnitude Contour
nexttile
contourf(R_grid*1e6, Z_grid*1e6, T_grid, 50, 'LineStyle', 'none');
colormap(gca, parula);
cb = colorbar;
ylabel(cb, 'Temperature field', 'FontSize', 11);
hold on

% Overlay force vectors
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_total_r_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_total_z_grid(1:skip_force:end, 1:skip_force:end)*1e15, 2, 'w', 'LineWidth', 1);

xlabel('r [\mum]', 'FontSize', 12);
ylabel('z [\mum]', 'FontSize', 12);
axis equal tight;
grid on;
set(gca, 'FontSize', 11);
fig_dir='pwd'
fileName='forces'
savefigures(fig_dir, fileName);
