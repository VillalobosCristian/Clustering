clear all; clc; close all;

%% TEMPERATURE GRADIENT PLOT - Normalized by R_illum
%% Two figures: v_r and ∂T/∂r vs r/R_illum

%% FIGURE PARAMETERS
FIGURE_WIDTH = 5;
FIGURE_HEIGHT = 5;
FONT_SIZE_LABELS = 28;
FONT_SIZE_AXES = 28;
FONT_SIZE_LEGEND = 14;
LINE_WIDTH = 2.5;
AXES_WIDTH = 1.5;

%% SETUP
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');

%% DEFINE OBJECTIVES
objectives = {'60x', '40x', '20x'};
R_illums = [25, 35, 75];  % μm
colors_obj = [0.8 0.2 0.2;   % Red for 60x
              0.2 0.6 0.2;   % Green for 40x
              0.2 0.2 0.8];  % Blue for 20x

n_obj = length(objectives);


%% LOAD ALL DATASETS
fprintf('Select COMSOL CSV files in order:\n');
for i = 1:n_obj
    fprintf('  %d. %s objective (R_illum = %d μm)\n', i, objectives{i}, R_illums(i));
end
fprintf('\n');

all_data = cell(n_obj, 1);

for i_obj = 1:n_obj
    [file_name, path_name] = uigetfile('*.csv', ...
        sprintf('Select %s objective file', objectives{i_obj}));
    
    if isequal(file_name, 0)
        error('File selection cancelled');
    end
    
    filename = fullfile(path_name, file_name);
    fprintf('Loading %s: %s\n', objectives{i_obj}, file_name);
    
    % Load data
    data = readmatrix(filename, 'NumHeaderLines', 9);
    
    % Extract fields
    all_data{i_obj}.r = real(data(:, 1));
    all_data{i_obj}.z = real(data(:, 2));
    all_data{i_obj}.T = real(data(:, 3));
    all_data{i_obj}.u = real(data(:, 9));
    
    all_data{i_obj}.R_illum = R_illums(i_obj);
    all_data{i_obj}.objective = objectives{i_obj};
    
    fprintf('  Points: %d\n', length(all_data{i_obj}.r));
    fprintf('  ΔT = %.2f K\n\n', max(all_data{i_obj}.T) - min(all_data{i_obj}.T));
end

%% OUTPUT DIRECTORY
if ~exist('path_name', 'var')
    path_name = pwd;
end
output_dir = fullfile(path_name, 'figures_gradient');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% PROCESS EACH DATASET - EXTRACT PROFILES AT z=15 μm
% fprintf('Processing temperature gradient at z=15μm...\n\n');

z_target = 15e-6;  % 15 μm
z_tolerance = 2e-6;  % ±2 μm
n_bins = 40;

profiles = cell(n_obj, 1);

for i_obj = 1:n_obj
    data = all_data{i_obj};
    R_illum = R_illums(i_obj);
    
    % Extract layer at z=15 μm
    idx_layer = abs(data.z - z_target) < z_tolerance;
    
    if sum(idx_layer) < 10
        warning('Insufficient points at z=15μm for %s', objectives{i_obj});
        continue;
    end
    
    % Extract data
    r_layer = data.r(idx_layer);
    T_layer = data.T(idx_layer);
    u_layer = data.u(idx_layer);  % Radial velocity
    
    % Create bins
    r_max = max(r_layer);
    r_edges = linspace(0, r_max, n_bins+1);
    r_centers = (r_edges(1:end-1) + r_edges(2:end))/2;
    
    % Bin temperature and velocity
    T_mean = nan(n_bins, 1);
    u_mean = nan(n_bins, 1);
    
    for i_bin = 1:n_bins
        idx_bin = r_layer >= r_edges(i_bin) & r_layer < r_edges(i_bin+1);
        if sum(idx_bin) >= 2
            T_mean(i_bin) = mean(T_layer(idx_bin));
            u_mean(i_bin) = mean(u_layer(idx_bin));
        end
    end
    
    % Smooth
    T_mean = movmean(T_mean, 3, 'omitnan');
    u_mean = movmean(u_mean, 3, 'omitnan');
    
    % Calculate temperature gradient ∂T/∂r
    r_centers_col = r_centers(:);
    T_mean_col = T_mean(:);
    
    dr = diff(r_centers_col);
    dT = diff(T_mean_col);
    gradT_r = dT ./ dr;
    r_grad = (r_centers_col(1:end-1) + r_centers_col(2:end))/2;
    
    % Smooth gradient
    gradT_r = movmean(gradT_r, 3, 'omitnan');
    
    % Add r=0 point (gradient = 0 by symmetry at center)
    r_grad = [0; r_grad(:)];  % Prepend r=0
    gradT_r = [0; gradT_r(:)];  % Gradient is 0 at center
    
    % Store results
    profiles{i_obj}.r_centers = r_centers * 1e6;  % μm
    profiles{i_obj}.u_mean = u_mean * 1e6;  % μm/s
    profiles{i_obj}.r_normalized = (r_centers * 1e6) / R_illum;  % r/R_illum
    profiles{i_obj}.r_grad = r_grad * 1e6;  % μm
    profiles{i_obj}.gradT_r = gradT_r * 1e-6;  % K/μm
    profiles{i_obj}.r_grad_normalized = (r_grad * 1e6) / R_illum;  % r/R_illum
    profiles{i_obj}.R_illum = R_illum;
    
    fprintf('%s: Processed (%d points)\n', objectives{i_obj}, sum(~isnan(gradT_r)));
end

fprintf('\n');

%% ========================================================================
%% FIGURE 1: Radial Velocity v_r vs r/R_illum
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');
hold on;

for i_obj = 1:n_obj
    if isempty(profiles{i_obj}), continue; end
    
    % Get velocity data
    r_norm = profiles{i_obj}.r_normalized(:);
    v_r = profiles{i_obj}.u_mean(:);
    
    % Valid data
    valid = ~isnan(v_r) & ~isnan(r_norm);
    
    if sum(valid) < 2
        warning('Insufficient valid velocity points for %s', objectives{i_obj});
        continue;
    end
    
    % Plot velocity
    plot(r_norm(valid), v_r(valid), '-', ...
         'Color', colors_obj(i_obj,:), 'LineWidth', LINE_WIDTH, ...
         'DisplayName', sprintf('%s ($R_{\\mathrm{illum}} = %d$ $\\mu$m)', ...
         objectives{i_obj}, profiles{i_obj}.R_illum));
end

% Reference lines
yline(0, 'k--', 'LineWidth', 1.5);
xline(1.0, 'r--', 'LineWidth', 2, 'DisplayName', '$r = R_{\\mathrm{illum}}$');

% Shade heating zone (r < R_illum)
% y_range = ylim;
% fill([0, 1, 1, 0], [y_range(1), y_range(1), y_range(2), y_range(2)], ...
%      'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

xlabel('$r / R_{\mathrm{illum}}$', 'FontSize', FONT_SIZE_LABELS);
ylabel('$v_r$ ($\mu$m/s)', 'FontSize', FONT_SIZE_LABELS);
% title('Radial Velocity at $z = 15$ $\mu$m', 'FontSize', FONT_SIZE_LABELS-2);
% legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND);
grid off; box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);
xlim([0, 2.5]);

% Add annotation
% text(0.05, 0.95, 'Heating zone', ...
%      'Units', 'normalized', 'FontSize', FONT_SIZE_LEGEND, ...
%      'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', 'k');

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Velocity_vs_R_illum');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Velocity_vs_R_illum.pdf\n');

%% ========================================================================
%% FIGURE 2: Temperature Gradient ∂T/∂r vs r/R_illum
%% ========================================================================

figure('Position', [100, 100, 900, 700], 'Color', 'w');
hold on;

for i_obj = 1:n_obj
    if isempty(profiles{i_obj}), continue; end
    
    % Get gradient data
    r_norm = profiles{i_obj}.r_grad_normalized(:);
    grad_T = profiles{i_obj}.gradT_r(:);
    
    % Valid data
    valid = ~isnan(grad_T) & ~isnan(r_norm);
    
    if sum(valid) < 2
        warning('Insufficient valid gradient points for %s', objectives{i_obj});
        continue;
    end
    
    % Plot gradient
    plot(r_norm(valid), grad_T(valid), '-', ...
         'Color', colors_obj(i_obj,:), 'LineWidth', LINE_WIDTH, ...
         'DisplayName', sprintf('%s ($R_{\\mathrm{illum}} = %d$ $\\mu$m)', ...
         objectives{i_obj}, profiles{i_obj}.R_illum));
end

% Reference lines
yline(0, 'k--', 'LineWidth', 1.5);
xline(1.0, 'r--', 'LineWidth', 2, 'DisplayName', '$r = R_{\\mathrm{illum}}$');

% Shade heating zone (r < R_illum)
% y_range = ylim;
% fill([0, 1, 1, 0], [y_range(1), y_range(1), y_range(2), y_range(2)], ...
%      'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

xlabel('$r / R_{\mathrm{illum}}$', 'FontSize', FONT_SIZE_LABELS);
ylabel('$\partial T/\partial r$ (K/$\mu$m)', 'FontSize', FONT_SIZE_LABELS);
% title('Temperature Gradient at $z = 15$ $\mu$m', 'FontSize', FONT_SIZE_LABELS-2);
% legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND);
grid off; box on;
set(gca, 'LineWidth', AXES_WIDTH, 'FontSize', FONT_SIZE_AXES);
xlim([0, 2.5]);

% Add annotation
% text(0.05, 0.95, 'Heating zone', ...
%      'Units', 'normalized', 'FontSize', FONT_SIZE_LEGEND, ...
%      'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', 'k');

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, FIGURE_WIDTH, FIGURE_HEIGHT]);

fileName = fullfile(output_dir, 'Temperature_Gradient_vs_R_illum');
export_fig(fileName, '-pdf', '-painters', '-r300', '-transparent');
fprintf('Saved: Temperature_Gradient_vs_R_illum.pdf\n');

%% ========================================================================
%% FINAL MESSAGE
