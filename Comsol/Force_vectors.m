
fprintf('Creating force vector visualization...\n');

%% Figure: Individual Force Vector Fields
figure('Name', 'Force Vector Fields Analysis', 'Position', [50, 50, 1800, 1000]);
set(gcf, 'Color', 'white');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

skip_force = 14;  % Plot every 8th vector for clarity

%% Subplot 1: Thermophoretic Force
nexttile
contourf(R_grid*1e6, Z_grid*1e6, sqrt(F_TPr_grid.^2 + F_TPz_grid.^2)*1e15, 20, 'LineStyle', 'none');
colormap(gca, hot);
hold on;
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_TPr_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_TPz_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       2, 'k', 'LineWidth', 1.5);
cb = colorbar;
ylabel(cb, '|F_{TP}| [fN]', 'FontSize', 10);
xlabel('r [μm]', 'FontSize', 11);
ylabel('z [μm]', 'FontSize', 11);
title('Thermophoretic Force', 'FontSize', 12, 'FontWeight', 'bold');
axis equal tight;
grid on;
set(gca, 'FontSize', 10);

%% Subplot 2: Drag Force
nexttile
contourf(R_grid*1e6, Z_grid*1e6, sqrt(F_Dr_grid.^2 + F_Dz_grid.^2)*1e15, 20, 'LineStyle', 'none');
colormap(gca, cool);
hold on;
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_Dr_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_Dz_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       2, 'k', 'LineWidth', 1.5);
cb = colorbar;
ylabel(cb, '|F_D| [fN]', 'FontSize', 10);
xlabel('r [μm]', 'FontSize', 11);
ylabel('z [μm]', 'FontSize', 11);
title('Drag Force', 'FontSize', 12, 'FontWeight', 'bold');
axis equal tight;
grid on;
set(gca, 'FontSize', 10);

%% Subplot 3: Gravitational Force (constant, pointing down)
nexttile
% Create a uniform field showing gravity
F_g_grid_plot = ones(size(R_grid)) * mean(F_g);
contourf(R_grid*1e6, Z_grid*1e6, abs(F_g_grid_plot)*1e15, 1, 'LineStyle', 'none');
colormap(gca, gray);
hold on;
% Plot downward arrows
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       zeros(size(R_grid(1:skip_force:end, 1:skip_force:end))), ...
       F_g_grid_plot(1:skip_force:end, 1:skip_force:end)*1e15, ...
       2, 'r', 'LineWidth', 1.5);
cb = colorbar;
ylabel(cb, '|F_g| [fN]', 'FontSize', 10);
xlabel('r [μm]', 'FontSize', 11);
ylabel('z [μm]', 'FontSize', 11);
title('Gravitational Force (uniform)', 'FontSize', 12, 'FontWeight', 'bold');
axis equal tight;
grid on;
set(gca, 'FontSize', 10);

%% Subplot 4: Total Force (TP + Drag + Gravity)
nexttile
F_total_mag_grid = sqrt(F_total_r_grid.^2 + F_total_z_grid.^2);
contourf(R_grid*1e6, Z_grid*1e6, F_total_mag_grid*1e15, 20, 'LineStyle', 'none');
colormap(gca, parula);
hold on;
quiver(R_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       Z_grid(1:skip_force:end, 1:skip_force:end)*1e6, ...
       F_total_r_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       F_total_z_grid(1:skip_force:end, 1:skip_force:end)*1e15, ...
       2, 'w', 'LineWidth', 1.5);
cb = colorbar;
ylabel(cb, '|F_{total}| [fN]', 'FontSize', 10);
xlabel('r [μm]', 'FontSize', 11);
ylabel('z [μm]', 'FontSize', 11);
title('Total Force (TP + Drag + Gravity)', 'FontSize', 12, 'FontWeight', 'bold');
axis equal tight;
grid on;
set(gca, 'FontSize', 10);

%% Subplot 5: Force Balance - Radial Component
nexttile
contourf(R_grid*1e6, Z_grid*1e6, F_total_r_grid*1e15, 30, 'LineStyle', 'none');
colormap(gca, redblue(256));
hold on;
% Add zero contour line
contour(R_grid*1e6, Z_grid*1e6, F_total_r_grid*1e15, [0 0], 'k', 'LineWidth', 3);
cb = colorbar;
ylabel(cb, 'F_{r,total} [fN]', 'FontSize', 10);
xlabel('r [μm]', 'FontSize', 11);
ylabel('z [μm]', 'FontSize', 11);
title('Radial Force Component (Black: F_r=0)', 'FontSize', 12, 'FontWeight', 'bold');
axis equal tight;
grid on;
set(gca, 'FontSize', 10);
% Add interpretation
text(0.5, 0.95, 'Blue = INWARD (attracts to center)', ...
    'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'k');

%% Subplot 6: Force Balance - Vertical Component
nexttile
contourf(R_grid*1e6, Z_grid*1e6, F_total_z_grid*1e15, 30, 'LineStyle', 'none');
colormap(gca, redblue(256));
hold on;
% Add zero contour line
contour(R_grid*1e6, Z_grid*1e6, F_total_z_grid*1e15, [0 0], 'k', 'LineWidth', 3);
cb = colorbar;
ylabel(cb, 'F_{z,total} [fN]', 'FontSize', 10);
xlabel('r [μm]', 'FontSize', 11);
ylabel('z [μm]', 'FontSize', 11);
title('Vertical Force Component (Black: F_z=0)', 'FontSize', 12, 'FontWeight', 'bold');
axis equal tight;
grid on;
set(gca, 'FontSize', 10);
% Add interpretation
text(0.5, 0.95, 'Red = UPWARD | Blue = DOWNWARD', ...
    'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'k');

%% Print force statistics
fprintf('\n==============================================\n');
fprintf('   FORCE VECTOR STATISTICS\n');
fprintf('==============================================\n\n');

fprintf('Thermophoretic Force:\n');
fprintf('  Max magnitude: %.3f fN\n', max(sqrt(F_TPr_grid(:).^2 + F_TPz_grid(:).^2))*1e15);
fprintf('  Mean magnitude: %.3f fN\n', mean(sqrt(F_TPr_grid(:).^2 + F_TPz_grid(:).^2))*1e15);
fprintf('  Radial component range: %.3f to %.3f fN\n', ...
    min(F_TPr_grid(:))*1e15, max(F_TPr_grid(:))*1e15);

fprintf('\nDrag Force:\n');
fprintf('  Max magnitude: %.3f fN\n', max(sqrt(F_Dr_grid(:).^2 + F_Dz_grid(:).^2))*1e15);
fprintf('  Mean magnitude: %.3f fN\n', mean(sqrt(F_Dr_grid(:).^2 + F_Dz_grid(:).^2))*1e15);
fprintf('  Radial component range: %.3f to %.3f fN\n', ...
    min(F_Dr_grid(:))*1e15, max(F_Dr_grid(:))*1e15);

fprintf('\nGravitational Force:\n');
fprintf('  Constant: %.3f fN (always downward)\n', mean(F_g)*1e15);

fprintf('\nTotal Force:\n');
fprintf('  Max magnitude: %.3f fN\n', max(sqrt(F_total_r_grid(:).^2 + F_total_z_grid(:).^2))*1e15);
fprintf('  Mean magnitude: %.3f fN\n', mean(sqrt(F_total_r_grid(:).^2 + F_total_z_grid(:).^2))*1e15);
fprintf('  Radial component range: %.3f to %.3f fN\n', ...
    min(F_total_r_grid(:))*1e15, max(F_total_r_grid(:))*1e15);
fprintf('  Vertical component range: %.3f to %.3f fN\n', ...
    min(F_total_z_grid(:))*1e15, max(F_total_z_grid(:))*1e15);

fprintf('\n==============================================\n');
fprintf('Force visualization complete!\n');
fprintf('==============================================\n\n');

%% Custom colormap function (red-white-blue)
function cmap = redblue(m)
    if nargin < 1
        m = 256;
    end
    % Red-White-Blue diverging colormap
    r = [ones(1,m/2), linspace(1,0,m/2)]';
    g = [linspace(0,1,m/2), linspace(1,0,m/2)]';
    b = [linspace(0,1,m/2), ones(1,m/2)]';
    cmap = [r, g, b];
end