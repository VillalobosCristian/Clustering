%% ========================================================================
%  Simple COMSOL Data Plotter (MATLAB Version)
%  - Velocity field (2D)
%  - Forces on particle
%  - Radial velocity profile <v_r>
%  Cristian
% ========================================================================

clear; clc; close all;

%% FIGURE PARAMS
set(0,'defaulttextInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

FIGURE_DPI = 300;

%% ========================================================================
%  COMPLEX PARSER
% ========================================================================
function val = parse_complex(s)
    % Convert to char if it's a string or cell
    if iscell(s)
        s = s{1};
    elseif isstring(s)
        s = char(s);
    end

    if isempty(s)
        val = NaN;
        return;
    end

    s = strtrim(s);

    % If complex number
    if contains(s,'i')
        s = erase(s,'i');

        % Find all + or - not preceded by e or E
        idx = regexp(s,'(?<![eE])[+-]');
        if ~isempty(idx)
            idx = idx(1);   % take first
            s = s(1:idx-1); % extract before
        end
    end

    % Convert to number
    val = str2double(s);
end

%% ========================================================================
%  LOAD COMSOL CSV FILE
% ========================================================================
function df = load_data(filename)
    fprintf("Loading %s...\n", filename);

    opts = detectImportOptions(filename);
    opts.DataLines = [10 Inf];     % skip 9 header lines
    opts = setvartype(opts, repmat("char",1,10));

    T = readtable(filename, opts);

    % Parse complex values column-by-column
    for k = 1:width(T)
        T.(k) = arrayfun(@parse_complex, T.(k));
    end

    T.Properties.VariableNames = ...
        {'r','z','T','F_TPr','F_Dr','F_TPz','F_Dz','F_g','u','w'};

    % Unit conversions
    T.r_um = T.r * 1e6;
    T.z_um = T.z * 1e6;
    T.u_um = T.u * 1e6;
    T.w_um = T.w * 1e6;

    T.F_r_fN = (T.F_TPr + T.F_Dr) * 1e15;
    T.F_z_fN = (T.F_TPz + T.F_Dz + T.F_g) * 1e15;

    T.Delta_T = T.T - 293.15;

    df = T;

    fprintf("  %d points loaded\n", height(df))
end

%% ========================================================================
%  RADIAL PROFILE <v_r>
% ========================================================================
function [r_centers, v_r] = get_radial_profile(df, zmax)

    edges = 0:2:80;
    r_centers = (edges(1:end-1) + edges(2:end))/2;

    mask = df.z_um < zmax;
    sub = df(mask,:);

    v_r = nan(size(r_centers));

    for i = 1:length(r_centers)
        mask_bin = sub.r_um >= edges(i) & sub.r_um < edges(i+1);
        if any(mask_bin)
            v_r(i) = mean(sub.u_um(mask_bin));
        end
    end
end

%% ========================================================================
% MAIN
% ========================================================================
filename = "Fields_35um.csv";  % change
R_illum = 35;

df = load_data(filename);


%% ========================================================================
% FIGURE 1 — VELOCITY FIELD
% ========================================================================
figure('Position',[50 50 1000 600]);

r_max = 80;  z_max = 25;
r_grid = linspace(0, r_max, 60);
z_grid = linspace(0.5, z_max, 40);
[R,Z] = meshgrid(r_grid, z_grid);

U = griddata(df.r_um, df.z_um, df.u_um, R, Z);
W = griddata(df.r_um, df.z_um, df.w_um, R, Z);
speed = sqrt(U.^2 + W.^2);

pcolor(R, Z, speed); shading interp; colormap viridis;
colorbar; title('|v| ($\mu$m/s)');
hold on;

streamslice(R, Z, U, W);
plot([R_illum R_illum], [0 z_max], 'r--','LineWidth',2)

xlabel('r ($\mu$m)');
ylabel('z ($\mu$m)');
title('Velocity Field');

saveas(gcf,'plot_velocity_field.png');


%% ========================================================================
% FIGURE 2 — FORCES
% ========================================================================
figure('Position',[50 50 1200 500]);

% Panel (a) — F_r
subplot(1,2,1)
F_r = griddata(df.r_um, df.z_um, df.F_r_fN, R, Z);
vmax = prctile(abs(df.F_r_fN),95);
pcolor(R,Z,F_r); shading interp; colormap thermal; caxis([-vmax vmax]);
colorbar
hold on
plot([R_illum R_illum],[0 z_max],'k--','LineWidth',2)
xlabel('r ($\mu$m)');
ylabel('z ($\mu$m)');
title('Radial Force $F_r$')

% Panel (b) — F_z
subplot(1,2,2)
F_z = griddata(df.r_um, df.z_um, df.F_z_fN, R, Z);
vmaxz = prctile(abs(df.F_z_fN),95);
pcolor(R,Z,F_z); shading interp; colormap thermal; caxis([-vmaxz vmaxz]);
colorbar
hold on
plot([R_illum R_illum],[0 z_max],'k--','LineWidth',2)
xlabel('r ($\mu$m)');
ylabel('z ($\mu$m)');
title('Vertical Force $F_z$')

saveas(gcf,'plot_forces.png');


%% ========================================================================
% FIGURE 3 — RADIAL VELOCITY PROFILE
% ========================================================================
figure('Position',[50 50 700 500]);

[r_prof, vr_prof] = get_radial_profile(df, 5);

plot(r_prof, vr_prof, 'o-','LineWidth',2); hold on
plot([R_illum R_illum], ylim, 'r--','LineWidth',2)

yline(0,'k-')
xlabel('r ($\mu$m)')
ylabel('$\langle v_r \rangle$ ($\mu$m/s)')
title('Radial Velocity Profile (z < 5 $\mu$m)')
grid on

saveas(gcf,'plot_vr_profile.png');


%% ========================================================================
% FIGURE 4 — MULTIPLE IRIS SIZES
% ========================================================================
figure('Position',[50 50 900 600]);

files = {
    '25 \mum', "Fields_25um.csv", 25, [0.8 0.2 0.2];
    '35 \mum', "Fields_35um.csv", 35, [0.2 0.6 0.2];
    '75 \mum', "Fields_75um.csv", 75, [0.2 0.3 0.8];
};

hold on
for i = 1:size(files,1)
    name  = files{i,1};
    fpath = files{i,2};
    Rval  = files{i,3};
    color = files{i,4};

    dfi = load_data(fpath);
    [r_i, vr_i] = get_radial_profile(dfi, 5);

    plot(r_i, vr_i, 'LineWidth',2, 'Color',color, 'DisplayName', name);
    xline(Rval,'--','Color',color,'LineWidth',1.2)
end

yline(0,'k-')
xlabel('r ($\mu$m)')
ylabel('$\langle v_r \rangle$ ($\mu$m/s)')
title('Radial Velocity Profile — All Iris Sizes')
legend show
grid on

saveas(gcf,'plot_vr_comparison.png');

fprintf("\nDone! All plots saved.\n");
