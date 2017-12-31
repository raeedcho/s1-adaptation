function bumpcurl_analysis(trial_data)
% BUMPCURL_ANALYSIS set of analyses to examine bumpcurl data

%% Look at adaptation behavior
clearvars -except trial_data
% get only rewards
[~,td] = getTDidx(trial_data,'result','R');
td = removeBadTrials(td,struct('remove_nan_idx',false));
td = getSpeed(td);
td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime'));
% td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','which_method','thresh','s_thresh',20));

metric = getLearningMetrics(td,struct('which_metric','angle','fit_bl_ref_curve',true,'vel_or_pos','pos','time_window',{{'idx_movement_on',0;'idx_peak_speed',0}}));
% metric = getLearningMetrics(td,struct('which_metric','corr'));
% metric = getLearningMetrics(td,struct('which_metric','curvature'));

bl_idx = getTDidx(td,'epoch','BL');
ad_idx = getTDidx(td,'epoch','AD');
wo_idx = getTDidx(td,'epoch','WO');

% plot metrics
figure
plot(metric)
hold on
plot(repmat(ad_idx(1),2,1),[-1 1],'k--','linewidth',3)
plot(repmat(wo_idx(1),2,1),[-1 1],'k--','linewidth',3)
plot([0 wo_idx(end)],[0 0],'k-','linewidth',2)

%% Visualize baseline data
clearvars -except trial_data
[~,td] = getTDidx(trial_data,'result','R','epoch','BL');
td = removeBadTrials(td,struct('remove_nan_idx',false));
td = getSpeed(td);
td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime'));

visData(td,struct('trials_to_plot',getTDidx(td,'bumpDir',90)));

%% Look at tuning of neurons in baseline and washout
% upshot -
%   - Many neurons are tuned to bumps
%   - Tuning doesn't change between baseline and washout
clearvars -except trial_data
% get only rewards
[~,td] = getTDidx(trial_data,'result','R');
td = removeBadTrials(td,struct('remove_nan_idx',false));

% scale spikes by binsize
for i = 1:length(td)
    td(i).S1_spikes = td(i).S1_spikes/td(i).bin_size;
end
% split into baseline and washout
[~,td_bl] = getTDidx(td,'epoch','BL');
[~,td_wo] = getTDidx(td,'epoch','WO');

% trim to just bumps
td_bl = trimTD(td_bl,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_wo = trimTD(td_wo,{'idx_bumpTime',0},{'idx_bumpTime',15});

% opensim_idx = find(contains(td(1).opensim_names,'_vel'));
opensim_idx = find(contains(td(1).opensim_names,'_muscVel'));
signal_name = 'opensim';
num_boots = 1000;
opensimPDs{1} = getTDPDs(td_bl,struct('out_signals',{{signal_name,opensim_idx}},'out_signal_names',{td(1).opensim_names(opensim_idx)},...
                                          'distribution','normal','move_corr','vel','num_boots',num_boots));
opensimPDs{2} = getTDPDs(td_wo,struct('out_signals',{{signal_name,opensim_idx}},'out_signal_names',{td(1).opensim_names(opensim_idx)},...
                                          'distribution','normal','move_corr','vel','num_boots',num_boots));
spikesPDs{1} = getTDPDs(td_bl,struct('out_signals',{{'S1_spikes'}},'out_signal_names',{td(1).S1_unit_guide},...
                                          'distribution','poisson','move_corr','vel','num_boots',num_boots));
spikesPDs{2} = getTDPDs(td_wo,struct('out_signals',{{'S1_spikes'}},'out_signal_names',{td(1).S1_unit_guide},...
                                          'distribution','poisson','move_corr','vel','num_boots',num_boots));

% get tuning curves
num_bins = 8;
[opensimCurves{1},bins] = getTuningCurves(td_bl,struct('out_signals',{{signal_name,opensim_idx}},'out_signal_names',{td(1).opensim_names(opensim_idx)},'num_bins',num_bins));
opensimCurves{2} = getTuningCurves(td_wo,struct('out_signals',{{signal_name,opensim_idx}},'out_signal_names',{td(1).opensim_names(opensim_idx)},'num_bins',num_bins));
spikesCurves{1} = getTuningCurves(td_bl,struct('out_signals',{{'S1_spikes'}},'out_signal_names',{td(1).S1_unit_guide},'num_bins',num_bins));
spikesCurves{2} = getTuningCurves(td_wo,struct('out_signals',{{'S1_spikes'}},'out_signal_names',{td(1).S1_unit_guide},'num_bins',num_bins));

% plot muscle baseline tuning against washout tuning
isTuned_params = struct('move_corr','vel','CIthresh',pi/4);
isTuned = checkIsTuned(opensimPDs{1},isTuned_params)...
            & checkIsTuned(opensimPDs{2},isTuned_params);
h1 = figure;
subplot(1,2,2)
errorbar(opensimPDs{1}.velPD(isTuned),opensimPDs{2}.velPD(isTuned),...
        minusPi2Pi(opensimPDs{2}.velPD(isTuned)-opensimPDs{2}.velPDCI(isTuned,1)),...
        minusPi2Pi(opensimPDs{2}.velPDCI(isTuned,1)-opensimPDs{2}.velPD(isTuned)),...
        minusPi2Pi(opensimPDs{1}.velPD(isTuned)-opensimPDs{1}.velPDCI(isTuned,1)),...
        minusPi2Pi(opensimPDs{1}.velPDCI(isTuned,1)-opensimPDs{1}.velPD(isTuned)),...
        'ro','linewidth',2)
hold on
plot([-pi pi],[-pi pi],'--k','linewidth',2)
xlabel('Baseline PD')
ylabel('Washout PD')
title('Muscle PDs')
set(get(gca,'ylabel'),'rotation',0,'horizontalalignment','right')
set(gca,'box','off','tickdir','out','xlim',[-pi pi],'ylim',[-pi pi],'xtick',[-pi pi],'ytick',[-pi pi],'xticklabel',{'-\pi','\pi'},'yticklabel',{'-\pi','\pi'})
axis equal

% plot all tuning
figure
compareTuning(opensimCurves,opensimPDs,bins,find(isTuned))


% plot neuron baseline tuning against washout tuning (same figure)
figure(h1)
isTuned_params = struct('move_corr','vel','CIthresh',pi/3);
isTuned = checkIsTuned(spikesPDs{1},isTuned_params)...
            & checkIsTuned(spikesPDs{2},isTuned_params);
subplot(1,2,1)
errorbar(spikesPDs{1}.velPD(isTuned),spikesPDs{2}.velPD(isTuned),...
        (spikesPDs{2}.velPD(isTuned)-spikesPDs{2}.velPDCI(isTuned,1)),...
        (spikesPDs{2}.velPDCI(isTuned,1)-spikesPDs{2}.velPD(isTuned)),...
        (spikesPDs{1}.velPD(isTuned)-spikesPDs{1}.velPDCI(isTuned,1)),...
        (spikesPDs{1}.velPDCI(isTuned,1)-spikesPDs{1}.velPD(isTuned)),...
        'mo','linewidth',2)
hold on
plot([-pi pi],[-pi pi],'--k','linewidth',2)
xlabel('Baseline PD')
ylabel('Washout PD')
title('Neural PDs')
set(get(gca,'ylabel'),'rotation',0,'horizontalalignment','right')
set(gca,'box','off','tickdir','out','xlim',[-pi pi],'ylim',[-pi pi],'xtick',[-pi pi],'ytick',[-pi pi],'xticklabel',{'-\pi','\pi'},'yticklabel',{'-\pi','\pi'})
axis equal

% plot all tuning
figure
compareTuning(spikesCurves,spikesPDs,bins,find(isTuned))

%% Fit GLM to just kinematics during bumps in baseline and evaluate R^2 in washout
% get only rewards
[~,td] = getTDidx(trial_data,'result','R');
td = removeBadTrials(td,struct('remove_nan_idx',false));

% split into baseline and washout
[~,td_bl] = getTDidx(td,'epoch','BL');
[~,td_wo] = getTDidx(td,'epoch','WO');

% trim to just bumps
td_bl = trimTD(td_bl,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_bl = binTD(td_bl,5);
td_wo = trimTD(td_wo,{'idx_bumpTime',0},{'idx_bumpTime',15});
td_wo = binTD(td_wo,5);

% fit glm to baseline
[td_bl,bl_model] = getModel(td_bl,struct('model_type','glm',...
    'model_name','S1_bumps_handle','in_signals',{{'pos';'vel';'force'}},...
    'out_signals',{'S1_spikes'}));
[td_wo] = getModel(td_bl,bl_model);

% evaluate washout period
eval_params = bl_model;
eval_params.eval_metric = 'pr2';
eval_params.trial_idx = 1:5:length(td_wo);
wo_eval = evalModel(td_wo,eval_params);

% plot resulting evals
figure
trial_num = eval_params.trial_idx(1:end-1);
for i=1:size(wo_eval,2)
    clf;
    single_eval = squeeze(wo_eval(:,i,:));
    patch([trial_num';flipud(trial_num')],[single_eval(:,1);flipud(single_eval(:,2))],[1,0,0]);
    waitforbuttonpress
end


%% Fit GLM to kinematics at end of adaptation and evaluate R^2 going backwards
