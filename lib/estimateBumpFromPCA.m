function td = estimateBumpFromPCA(trial_data)
% ESTIMATEBUMPFROMPCA Estimates the direction of a bump from S1 PC space

[~,td] = getTDidx(trial_data,'result','R');

% for all directions
%   compute distribution of distances to clusters in baseline condition
ubumpDirs = unique([td.bumpDir]);

for iDir = 1:length(ubumpDirs)
    [~,td_bl] = getTDidx(td,'epoch','BL','bumpDir',ubumpDirs(iDir));
    
    for iTrial = 1:length(td_bl)
        