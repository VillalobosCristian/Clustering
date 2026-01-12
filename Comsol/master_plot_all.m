close all
clear all

filename = 'Fields_75um.csv';

if ~exist(filename, 'file')
    error('File not found! Please update the filename variable.');
end

%% Load Data
fprintf('Loading data...\n');

fid = fopen(filename, 'r');
for i = 1:9
    fgetl(fid);
end
fclose(fid);

data = readmatrix(filename, 'NumHeaderLines', 9);

r = data(:, 1);
z = data(:, 2);
T = data(:, 3);
F_TPr = real(data(:, 4));
F_Dr = real(data(:, 5));
F_TPz = real(data(:, 6));
F_Dz = real(data(:, 7));
F_g = real(data(:, 8));

F_total_r = F_TPr + F_Dr;
F_total_z = F_TPz + F_Dz + F_g*1e-2;

fprintf('Data loaded: %d points\n', length(r));

%% Interpolate
fprintf('Interpolating...\n');

r_min = min(r); r_max = max(r);
z_min = min(z); z_max = max(z);

grid_resolution = 200;
r_grid = linspace(r_min, r_max, grid_resolution);
z_grid = linspace(z_min, z_max, grid_resolution);
[R_grid, Z_grid] = meshgrid(r_grid, z_grid);

T_grid = griddata(r, z, T, R_grid, Z_grid, 'natural');
F_total_r_grid = griddata(r, z, F_total_r, R_grid, Z_grid, 'natural');
F_total_z_grid = griddata(r, z, F_total_z, R_grid, Z_grid, 'natural');

fprintf('Interpolation complete!\n\n');

%% Figure: Color-Coded Force Vectors by Vertical Component
fprintf('Creating color-coded force vectors...\n');

figure('Name', 'Total Force Analysis', 'Position', [150, 150, 1600, 700]);
set(gcf, 'Color', 'white');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

skip_force = 8;

% Prepare data for plotting
R_plot = R_grid(1:skip_force:end, 1:skip_force:end);
Z_plot = Z_grid(1:skip_force:end, 1:skip_force:end);
F_r_plot = F_total_r_grid(1:skip_force:end, 1:skip_force:end);
F_z_plot = F_total_z_grid(1:skip_force:end, 1:skip_force:end);

%% Subplot 1: Color-coded vectors only
nexttile
hold on;

% Plot each arrow colored by F_z
for i = 1:size(R_plot, 1)
    for j = 1:size(R_plot, 2)
        if ~isnan(F_z_plot(i,j))
            % Color based on vertical component
            if F_z_plot(i,j) < 0
                % Downward (toward surface) - RED/ORANGE
                arrow_color = [0.9 0.3 0.2];  % Red
                line_width = 2;
            else
                % Upward (away from surface) - BLUE
                arrow_color = 'k';  % Blue
                line_width = 1.5;
            end
            
            % Plot arrow
            quiver(R_plot(i,j)*1e6, Z_plot(i,j)*1e6, ...
                   F_r_plot(i,j)*1e15, F_z_plot(i,j)*1e15, ...
                   2, 'Color', arrow_color, 'LineWidth', 1, 'MaxHeadSize', 0.5);
        end
    end
end

xlabel('r [\mum]', 'FontSize', 12);
ylabel('z [\mum]', 'FontSize', 12);
axis equal tight;
% grid on;
box on
set(gca, 'FontSize', 11);

% Add legend
% text(0.02, 0.95, '{\bf Red:} Downward $\downarrow$ (F_z < 0)', ...
%     'Units', 'normalized', 'FontSize', 11, 'Interpreter', 'latex', ...
%     'BackgroundColor', [0.9 0.3 0.2], 'ForeGroundColor', 'w', ...
%     'EdgeColor', 'k', 'LineWidth', 1.5, 'VerticalAlignment', 'top');
% text(0.02, 0.88, '{\bf Blue:} Upward $\uparrow$ (F_z > 0)', ...
%     'Units', 'normalized', 'FontSize', 11, 'Interpreter', 'latex', ...
%     'BackgroundColor', [0.2 0.4 0.9], 'ForeGroundColor', 'w', ...
%     'EdgeColor', 'k', 'LineWidth', 1.5, 'VerticalAlignment', 'top');

%% Subplot 2: Color-coded vectors on temperature background
nexttile
hold on;

% Temperature background
contourf(R_grid*1e6, Z_grid*1e6, T_grid, 50, 'LineStyle', 'none');
colormap(gca, turbo);

% Overlay colored force vectors
for i = 1:size(R_plot, 1)
    for j = 1:size(R_plot, 2)
        if ~isnan(F_z_plot(i,j))
            % Color based on vertical component
            if F_z_plot(i,j) < 0
                % Downward (toward surface) - RED/ORANGE (stands out on blue background)
                arrow_color = [1 0.3 0.2];  % Bright red
                line_width = 2.5;
            else
                % Upward (away from surface) - WHITE (contrasts with blue)
                arrow_color = [1 1 1];  % White
                line_width = 1.5;
            end
            
            % Plot arrow
            quiver(R_plot(i,j)*1e6, Z_plot(i,j)*1e6, ...
                   F_r_plot(i,j)*1e15, F_z_plot(i,j)*1e15, ...
                   2, 'Color', arrow_color, 'LineWidth', 1, 'MaxHeadSize', 0.5);
        end
    end
end

xlabel('r [\mum]', 'FontSize', 12);
ylabel('z [\mum]', 'FontSize', 12);
axis equal tight;
% 
box on
set(gca, 'FontSize', 11);

% Colorbar for temperature
cb = colorbar;
ylabel(cb, 'Temperature field', 'FontSize', 11);



%% Optional: Save figure
% Uncomment to save
% fig_dir = pwd;
% fileName = 'force_vectors_vertical_color';
% saveas(gcf, fullfile(fig_dir, [fileName '.png']));
% saveas(gcf, fullfile(fig_dir, [fileName '.fig']));
% fprintf('Figure saved to: %s\n', fig_dir);
%% Interpolate
fprintf('Interpolating...\n');

r_min = min(r); r_max = max(r);
z_min = min(z); z_max = max(z);

grid_resolution = 200;
r_grid = linspace(r_min, r_max, grid_resolution);
z_grid = linspace(z_min, z_max, grid_resolution);
[R_grid, Z_grid] = meshgrid(r_grid, z_grid);

T_grid = griddata(r, z, T, R_grid, Z_grid, 'natural');
F_total_r_grid = griddata(r, z, F_total_r, R_grid, Z_grid, 'natural');
F_total_z_grid = griddata(r, z, F_total_z, R_grid, Z_grid, 'natural');

skip_force = 8;

% Prepare data for plotting
R_plot = R_grid(1:skip_force:end, 1:skip_force:end);
Z_plot = Z_grid(1:skip_force:end, 1:skip_force:end);
F_r_plot = F_total_r_grid(1:skip_force:end, 1:skip_force:end);
F_z_plot = F_total_z_grid(1:skip_force:end, 1:skip_force:end);
%% Subplot 2: Temperature + Gradient magnitude surface + Force vectors
% nexttile

figure
hold on;

% --- Temperature background ---
contourf(R_grid*1e6, Z_grid*1e6, T_grid, 50, 'LineStyle', 'none');
colormap(gca, turbo);

% --- Compute temperature gradients ---
[dT_dr, dT_dz] = gradient(T_grid, r_grid(2)-r_grid(1), z_grid(2)-z_grid(1));
grad_T_magnitude = sqrt(dT_dr.^2 + dT_dz.^2);

% --- Surface plot of gradient magnitude (color = |gradT|) ---
s = surf(R_grid*1e6, Z_grid*1e6, zeros(size(grad_T_magnitude)), ... % flat in z
         grad_T_magnitude, ...                                        % color = gradT
         'EdgeColor', 'none', 'FaceAlpha', 0.7);                      % semi-transparent
colormap(gca, hot);  % colormap for gradient

% --- Overlay colored force vectors ---
for i = 1:size(R_plot, 1)
    for j = 1:size(R_plot, 2)
        if ~isnan(F_z_plot(i,j))
            if F_z_plot(i,j) < 0
                arrow_color = [1 0.3 0.2];  % Red/orange downward
            else
                arrow_color = [1 1 1];      % White upward
            end
            quiver(R_plot(i,j)*1e6, Z_plot(i,j)*1e6, ...
                   F_r_plot(i,j)*1e15, F_z_plot(i,j)*1e15, ...
                   2, 'Color', arrow_color, 'LineWidth', 1, 'MaxHeadSize', 0.5);
        end
    end
end

% --- Axis labels & formatting ---
xlabel('r [\mum]', 'FontSize', 12);
ylabel('z [\mum]', 'FontSize', 12);
axis equal tight;
box on
set(gca, 'FontSize', 11);

% --- Colorbar for gradient magnitude ---
cb = colorbar;
ylabel(cb, '|\nabla T| [K/m]', 'FontSize', 11);

title('Temperature + |∇T| Surface + Force Vectors', 'FontSize', 12);
