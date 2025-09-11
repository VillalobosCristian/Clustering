clear all; clc; close all;

%% Load CSV data
[file_name, path_name] = uigetfile('*.csv', 'Select CSV file');
if isequal(file_name, 0), return; end

data = readtable(fullfile(path_name, file_name));
data = data(4:end,:); % Remove TrackMate headers

% Convert to numeric if needed
if iscell(data.TRACK_ID), data.TRACK_ID = cellfun(@str2double, data.TRACK_ID); end
if iscell(data.POSITION_X), data.POSITION_X = cellfun(@str2double, data.POSITION_X); end
if iscell(data.POSITION_Y), data.POSITION_Y = cellfun(@str2double, data.POSITION_Y); end
if iscell(data.POSITION_T), data.POSITION_T = cellfun(@str2double, data.POSITION_T); end

sorted_data = sortrows(data, {'TRACK_ID', 'POSITION_T'});
unique_tracks = unique(sorted_data.TRACK_ID);

%% Extract trajectories
min_track_length = 100;
X = {}; Y = {}; T = {};
count = 0;

for i = 1:length(unique_tracks)
    track_data = sorted_data(sorted_data.TRACK_ID == unique_tracks(i), :);
    if height(track_data) >= min_track_length
        count = count + 1;
        X{count} = track_data.POSITION_X;
        Y{count} = track_data.POSITION_Y;
        T{count} = track_data.POSITION_T;
    end
end

targetTime = 1000;  
calib = 1.0;      

positions = [];
for i = 1:length(X)
    timeIdx = find(T{i} == targetTime, 1);
    if ~isempty(timeIdx)
        positions = [positions; X{i}(timeIdx), Y{i}(timeIdx)];
    end
end
positions = positions * calib;

N = size(positions,1);
k = 6;  % six nearest neighbors

% find the k+1 nearest neighbors (first neighbor is self, so we'll skip it)
[idxs, dists] = knnsearch(positions, positions, 'K', k+1);

psi6 = zeros(N,1);

for i = 1:N
    neigh = idxs(i,2:end);   % skip the first column = self
    dx = positions(neigh,1) - positions(i,1);
    dy = positions(neigh,2) - positions(i,2);
    thetas = atan2(dy, dx);   % bond angles from i to each neighbor
    
    psi6(i) = abs( mean( exp(6i * thetas) ) );
end

% --- now plot colored by |psi6| ---
figure; hold on;
colormap(jet);
scatter(positions(:,1), positions(:,2), 50, psi6, 'filled');
colorbar;
caxis([0 1]);
axis equal; box on;
xlabel('x (\mum)');
ylabel('y (\mum)');
title('Local Hexatic Order');
%%
