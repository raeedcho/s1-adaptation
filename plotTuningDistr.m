function plotTuningDistr(trial_data)
% plot distribution of preferred directions in baseline

[~,td] = getTDidx(trial_data,'epoch','BL','result','R');

td = truncateAndBin(td,4,{'idx_bumpTime',0},{'idx_bumpTime',12});

% get tuning curves from unbinned data
sort_idx = find(td(1).S1_unit_guide(:,2));
[~,glm_info] = getModel(td,struct('model_type','glm','model_name','linmodel',...
                                  'in_signals',{'vel'},'out_signals',{{'S1_spikes',sort_idx}}));
PD = atan2(glm_info.b(3,:),glm_info.b(2,:));
offset = glm_info.b(1,:);
moddepth = sqrt(sum(glm_info.b(2:3,:).^2,1));

figure
h = polar(repmat(PD,2,1),[0;1]*ones(size(PD)),'ko-');
set(h,'linewidth',2,'markersize',10)
