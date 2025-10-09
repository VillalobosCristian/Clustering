clear all; clc; close all;

%% ------------------------- Load CSV data -------------------------------
[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), return; end

data = readtable(fullfile(path_name, file_name));
data = data(4:end,:); % Remove TrackMate headers

% Convert columns to numeric if needed
if iscell(data.TRACK_ID),    data.TRACK_ID    = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X),  data.POSITION_X  = cellfun(@str2double, data.POSITION_X); end
if iscell(data.POSITION_Y),  data.POSITION_Y  = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T),  data.POSITION_T  = cellfun(@str2double, data.POSITION_T); end

% Sort for easier processing
sorted_data   = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

%% ------------------------- Extract trajectories ------------------------
min_track_length = 10;
X = {}; Y = {}; T = {}; particle_ids = [];
count = 0;

for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    if height(track_data) >= min_track_length
        count = count + 1;
        X{count} = track_data.POSITION_X;
        Y{count} = track_data.POSITION_Y;
        T{count} = track_data.POSITION_T;
        particle_ids(count) = unique_tracks(i);
    end
end

nTracks = numel(X);
if nTracks == 0
    error('No tracks with at least %d points were found.', min_track_length);
end

%% ------------------------- Time grid for analysis ----------------------
% Collect all time points present in kept tracks
all_times = [];
for i = 1:nTracks, all_times = [all_times; T{i}(:)]; end
unique_times = unique(all_times);

time_step = 1;  % can subsample frames here to speed up
analysis_times = unique_times(1:time_step:end);
nT = numel(analysis_times);

%% ------------------------- Parameters ----------------------------------
calib = 1.0;                 % px -> µm
distance_threshold = 5;      % µm
persistence_frames = 5;      % min frames to be "stable"
tolerance_frames   = 5;      % allowed consecutive gaps (absent or disconnected)

%% ------------------------- Plot limits 
all_x = []; all_y = [];
for i = 1:nTracks
    all_x = [all_x; calib * X{i}(:)];
    all_y = [all_y; calib * Y{i}(:)];
end
x_padding = 0.05 * (max(all_x) - min(all_x));
y_padding = 0.05 * (max(all_y) - min(all_y));
x_limits = [min(all_x) - x_padding, max(all_x) + x_padding];
y_limits = [min(all_y) - y_padding, max(all_y) + y_padding];

%% ------------------------- Frame buckets (no per-frame find) -----------
frame_pts = cell(nT,1);   % Nx2 positions (µm)
frame_ids = cell(nT,1);   % corresponding global particle indices (1..nTracks)

for i = 1:nTracks
    if isempty(T{i}), continue; end
    [tf, loc] = ismember(T{i}, analysis_times); % exact match
    if any(tf)
        xi = calib * X{i}(tf);
        yi = calib * Y{i}(tf);
        li = loc(tf);
        for k = 1:numel(li)
            t = li(k);
            frame_pts{t} = [frame_pts{t}; xi(k), yi(k)];
            frame_ids{t} = [frame_ids{t}; i];
        end
    end
end

%% ------------------------- Build connection history --------------------
fprintf('Building particle connection history...\n');
particle_history = false(nTracks, nT); % logical: connected / not

haveRangesearch = exist('rangesearch', 'file') == 2;
r  = distance_threshold;
r2 = r*r;

for t = 1:nT
    P = frame_pts{t};
    ids = frame_ids{t};
    Np = size(P,1);
    if Np < 2
        continue
    end

    if haveRangesearch
        nbrs = rangesearch(P, P, r);
        isConn = cellfun(@(c) numel(c) > 1, nbrs); % any neighbor inside r
    else
        % Squared distance comparison (avoid sqrt)
        D2 = pdist2(P, P, 'squaredeuclidean');
        D2(1:Np+1:end) = inf; % ignore self
        isConn = any(D2 <= r2, 2);
    end

    particle_history(ids, t) = isConn;
end

%% ------------------------- Video writer --------------------------------
video_filename = [file_name(1:end-4), '_delaunay_persistence_animation.mp4'];
video_path = fullfile(path_name, video_filename);
v = VideoWriter(video_path, 'MPEG-4');
v.FrameRate = 10;
v.Quality   = 95;
open(v);

%% ------------------------- Figure & handles (reuse!) -------------------
fig = figure('Position', [100, 100, 1400, 600], 'Color', 'w', 'Renderer', 'opengl');

tiledlayout(fig, 1, 2, 'Padding','compact','TileSpacing','compact');

% Colors
disconnected_color   = [0.7, 0.7, 0.7];
transient_color      = [0.2, 0.4, 0.8];
stable_color         = [0.9, 0.2, 0.2];
delaunay_edge_color  = [0.8, 0.8, 0.8];
valid_triangle_color = [0.9, 0.2, 0.2];

% Left: Full Delaunay
ax1 = nexttile; hold(ax1, 'on'); axis(ax1, 'equal'); grid(ax1, 'on');
ax1.GridColor = [0.9 0.9 0.9]; ax1.GridAlpha = 0.6; ax1.Color = 'w';
xlim(ax1, x_limits); ylim(ax1, y_limits);
title(ax1, 'Full Delaunay Triangulation', 'FontSize', 14, 'FontWeight','normal', 'Color', [0.2 0.2 0.2]);
xlabel(ax1, '$x \,(\mu m)$', 'Interpreter','latex'); ylabel(ax1, '$y \,(\mu m)$', 'Interpreter','latex');

% Right: Distance-Filtered + Persistence
ax2 = nexttile; hold(ax2, 'on'); axis(ax2, 'equal'); grid(ax2, 'on');
ax2.GridColor = [0.9 0.9 0.9]; ax2.GridAlpha = 0.6; ax2.Color = 'w';
xlim(ax2, x_limits); ylim(ax2, y_limits);
title(ax2, 'Distance-Filtered + Persistence Analysis', 'FontSize', 14, 'FontWeight','normal', 'Color', [0.2 0.2 0.2]);
xlabel(ax2, '$x \,(\mu m)$', 'Interpreter','latex'); ylabel(ax2, '$y \,(\mu m)$', 'Interpreter','latex');

% Pre-create artists
hScatter1 = scatter(ax1, NaN, NaN, 60, [0 0 0], 'filled', ...
                    'MarkerEdgeColor','w','MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.8,'LineWidth',1);
hTri1     = patch(ax1,'Faces',zeros(0,3),'Vertices',zeros(0,2), ...
                  'FaceColor','none','EdgeColor',delaunay_edge_color,'LineWidth',0.8);

hScatter2 = scatter(ax2, NaN, NaN, 60, [0 0 0], 'filled', ...
                    'MarkerEdgeColor','w','MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.8,'LineWidth',1);
hTri2     = patch(ax2,'Faces',zeros(0,3),'Vertices',zeros(0,2), ...
                  'FaceColor','none','EdgeColor',valid_triangle_color,'LineWidth',1.5);

% Text labels (update only strings)
txt1 = text(ax1, 0.02, 0.98, '', 'Units','normalized','FontSize',11, ...
            'BackgroundColor','w','Margin',5,'EdgeColor',[0.8 0.8 0.8], 'VerticalAlignment','top');
txt2 = text(ax2, 0.02, 0.98, '', 'Units','normalized','FontSize',11, ...
            'BackgroundColor','w','Margin',5,'EdgeColor',[0.8 0.8 0.8], 'VerticalAlignment','top');
txt3 = text(ax2, 0.02, 0.88, '', 'Units','normalized','FontSize',10, ...
            'BackgroundColor','w','Margin',5,'EdgeColor',[0.8 0.8 0.8], 'VerticalAlignment','top');
txt4 = text(ax2, 0.02, 0.78, '', 'Units','normalized','FontSize',9, 'Color',[0.5 0.5 0.5], ...
            'BackgroundColor','w','Margin',5,'EdgeColor',[0.8 0.8 0.8], 'VerticalAlignment','top');

% Legend (once)
hStable    = scatter(ax2, NaN,NaN,60,stable_color,     'filled','MarkerEdgeColor','w');
hTransient = scatter(ax2, NaN,NaN,60,transient_color,  'filled','MarkerEdgeColor','w');
hDisc      = scatter(ax2, NaN,NaN,60,disconnected_color,'filled','MarkerEdgeColor','w');
legend(ax2, [hStable, hTransient, hDisc], {'Stable Cluster','Transient','Disconnected'}, 'Location','northeast','FontSize',9);

%% ------------------------- Persistence counters ------------------------
persist_cnt = zeros(nTracks,1,'uint16'); % consecutive connected frames (considering tolerance)
gap_cnt     = zeros(nTracks,1,'uint16'); % consecutive gaps (absent or disconnected)
tol  = uint16(tolerance_frames);
need = uint16(persistence_frames);

%% ------------------------- Animation loop ------------------------------
fprintf('Creating video: %s\n', video_filename);
fprintf('Processing %d time points...\n', nT);

for t_idx = 1:nT
    % Positions & ids present in this frame
    P   = frame_pts{t_idx};
    ids = frame_ids{t_idx};
    N   = size(P,1);

    % --- Integer-safe counter updates ---
    % 1) Everyone accumulates a gap by default
    gap_cnt = gap_cnt + uint16(1);

    if ~isempty(ids)
        % connectivity for present tracks
        isConn = particle_history(ids, t_idx);   % logical

        % 2) Reset gap counter for connected-present tracks
        if any(isConn)
            gap_cnt(ids(isConn)) = uint16(0);
        end
        % (present+disconnected keep the +1 we added above)

        % 3) If too many consecutive gaps for any track, reset persistence
        resetMask = gap_cnt > tol;
        if any(resetMask)
            persist_cnt(resetMask) = uint16(0);
        end

        % 4) Increase persistence only for connected-present tracks
        if any(isConn)
            persist_cnt(ids(isConn)) = persist_cnt(ids(isConn)) + uint16(1);
        end
    end

    % --- Classification for present particles ---
    isStable    = false(N,1);
    isTransient = false(N,1);
    isDisc      = false(N,1);
    if ~isempty(ids)
        isConn = particle_history(ids, t_idx);
        isStable    = isConn & (persist_cnt(ids) >= need);
        isTransient = isConn & ~isStable;
        isDisc      = ~isConn;
    end

    % Colors for present particles
    C = repmat(disconnected_color, N, 1);
    if any(isTransient), C(isTransient,:) = repmat(transient_color, sum(isTransient), 1); end
    if any(isStable),    C(isStable,:)    = repmat(stable_color,    sum(isStable),    1); end

    % --- Left panel: Full Delaunay triangulation ---
    if N >= 3
        tri = delaunay(P(:,1), P(:,2));
    else
        tri = zeros(0,3);
    end
    set(hTri1,     'Faces', tri, 'Vertices', P);
    set(hScatter1, 'XData', P(:,1), 'YData', P(:,2), 'CData', C);
    set(txt1,      'String', sprintf('%d triangles', size(tri,1)));

    % --- Right panel: distance-filtered triangles (vectorized) ---
    if ~isempty(tri)
        i12 = tri(:,[1 2]); i23 = tri(:,[2 3]); i31 = tri(:,[3 1]);
        e12 = sum((P(i12(:,1),:) - P(i12(:,2),:)).^2,2);
        e23 = sum((P(i23(:,1),:) - P(i23(:,2),:)).^2,2);
        e31 = sum((P(i31(:,1),:) - P(i31(:,2),:)).^2,2);
        maxe2 = max([e12 e23 e31], [], 2);
        keep = maxe2 <= r2;
        tri_valid = tri(keep,:);
    else
        tri_valid = zeros(0,3);
    end

    set(hTri2,     'Faces', tri_valid, 'Vertices', P);
    set(hScatter2, 'XData', P(:,1), 'YData', P(:,2), 'CData', C);

    % Stats text
    stable_count       = sum(isStable);
    transient_count    = sum(isTransient);
    disconnected_count = sum(isDisc);

    set(txt2, 'String', sprintf('%d valid triangles', size(tri_valid,1)));
    set(txt3, 'String', sprintf('Stable: %d | Transient: %d', stable_count, transient_count));
    set(txt4, 'String', sprintf('Threshold: %.1f \\mum | Persist: %d frames', distance_threshold, persistence_frames));

    % Window title
    fig.Name = sprintf('Time: %g | Stable: %d | Transient: %d | Disconnected: %d', ...
        analysis_times(t_idx), stable_count, transient_count, disconnected_count);

    % Draw & write frame
    drawnow limitrate nocallbacks;
    frame = getframe(fig); % for speed use: getframe(ax2)
    writeVideo(v, frame);

    if mod(t_idx, 10) == 0 || t_idx == nT
        fprintf('Progress: %d/%d frames (%.1f%%)\n', t_idx, nT, 100*t_idx/nT);
    end
end

% Close video and figure
close(v);
close(fig);

% Summary
fprintf('\nVideo saved successfully: %s\n', video_path);
fprintf('Video duration: %.2f seconds at %d fps\n', nT / v.FrameRate, v.FrameRate);
fprintf('\n=== PERSISTENCE ANALYSIS COMPLETE ===\n');
fprintf('Persistence threshold: %d frames\n', persistence_frames);
fprintf('Gap tolerance: %d frames\n', tolerance_frames);
fprintf('Distance threshold: %.1f μm\n', distance_threshold);
