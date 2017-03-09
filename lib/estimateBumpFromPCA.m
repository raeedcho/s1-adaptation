function td = estimateBumpFromPCA(trial_data)
% ESTIMATEBUMPFROMPCA Estimates the direction of a bump from S1 PC space

ndim = 10;

[~,td] = getTDidx(trial_data,'result','R');

% S1 pca
td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
td = truncateAndBin(td,{'idx_bumpTime',-10},{'idx_bumpTime',15});
[td,pca_info] = getPCA(td,struct('signals','S1_spikes'));

% for all directions
%   compute distribution of distances to clusters in baseline condition
[~,td_bl] = getTDidx(td,'epoch','BL');
ubumpDirs = sort(unique([td.bumpDir]));
bl_clust = cell(8,1);

% compose clusters
for iDir = 1:length(ubumpDirs)
    [~,td_temp] = getTDidx(td_bl,'bumpDir',ubumpDirs(iDir));
    bl_clust{iDir} = cat(1,td_temp.S1_pca(20,1:ndim));
end
