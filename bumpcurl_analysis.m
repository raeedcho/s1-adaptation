function bumpcurl_analysis(trial_data)
% BUMPCURL_ANALYSIS set of analyses to examine bumpcurl data

%% Look at adaptation behavior
    % Han_20171206
    %   - Used angle metric
    %   - Fit baseline reference with 3 sine functions
    %   - clear washout, less clear adaptation
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

    % scale spikes by binsize
    for i = 1:length(td)
        td(i).S1_spikes = td(i).S1_spikes/td(i).bin_size;
    end

    [~,td_bl] = getTDidx(td,'epoch','BL');
    td_bl = trimTD(td_bl,{'idx_bumpTime',0},{'idx_bumpTime',15});
    
    % visData(td_bl,struct('trials_to_plot',getTDidx(td_bl,'bumpDir',90)));

    % get PEH for all neurons and all directions
    td_bl = trialAverage(td_bl,'bumpDir');

    % plot PEH for each direction on super plot for each neuron
    % num_neurons = size(td_bl(1).S1_spikes,2);
    num_neurons = 10;
    num_dirs = length(td_bl);
    % get max FRs
    spikes = cat(1,td_bl.S1_spikes);
    maxFRs = max(spikes,[],1);
    for fignum = 1:4
        figure(fignum)
        for neuron_ctr = 1:num_neurons
            true_neur_num = num_neurons*(fignum-1)+neuron_ctr;
            % do some sort of figure
            for dir_ctr = 1:num_dirs
                % get plotnumber
                plotnum = (1+num_dirs)*(neuron_ctr-1)+dir_ctr;
                subplot(num_neurons,num_dirs+1,plotnum)
                plot(td_bl(dir_ctr).S1_spikes(:,true_neur_num),'linewidth',2)
                set(gca,'ylim',[0 maxFRs(true_neur_num)])
            end
        end
    end

    td_bl = binTD(td_bl,15);
    dirs = repmat(cat(1,td_bl.bumpDir),2,1);
    spikes = repmat(cat(1,td_bl.S1_spikes),2,1);
    for fignum = 1:4
        figure(fignum)
        for neuron_ctr = 1:num_neurons
            true_neur_num = num_neurons*(fignum-1)+neuron_ctr;
            plotnum = (1+num_dirs)*(neuron_ctr);
            subplot(num_neurons,num_dirs+1,plotnum)
            h = polar(dirs*pi/180,spikes(:,true_neur_num));
            set(h,'linewidth',2);
        end
    end

%% Look at tuning of neurons in baseline and washout
    % Han_20171206
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

%% Decoder analysis - train decoder for x and y velocity on neurons in baseline
    % Han_20171206
    %   - notes go here
    clearvars -except trial_data
    % get only rewards
    [~,td] = getTDidx(trial_data,'result','R');
    td = removeBadTrials(td,struct('remove_nan_idx',false));
    td = getSpeed(td);
    td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime'));

    % get baseline and washout
    [~,td_bl] = getTDidx(td,'epoch','BL');
    [~,td_wo] = getTDidx(td,'epoch','WO');

    % trim to bumps + 150 ms
    td_bl = trimTD(td_bl,{'idx_bumpTime',0},{'idx_bumpTime',15});
    % td_bl = trimTD(td_bl,{'idx_movement_on',0},{'idx_movement_on',30});
    td_bl = binTD(td_bl,15);
    td_wo = trimTD(td_wo,{'idx_bumpTime',0},{'idx_bumpTime',15});
    % td_wo = trimTD(td_wo,{'idx_movement_on',0},{'idx_movement_on',30});
    td_wo = binTD(td_wo,15);

    % crossvalidation sets to look at prediction capabilities
    [train_idx,test_idx] = crossvalind('HoldOut',length(td_bl),0.2);

    % Get spikes and velocities into matrices
    bl_vel = cat(1,td_bl.vel);
    wo_vel = cat(1,td_wo.vel);
    bl_dir = atan2(bl_vel(:,2),bl_vel(:,1));
    wo_dir = atan2(wo_vel(:,2),wo_vel(:,1));

    % fit models
    [td_bl,bl_linmodel] = getModel(td_bl,struct('model_type','linmodel',...
        'model_name','vel_from_spikes','out_signals',{{'vel'}},...
        'in_signals',{'S1_spikes'},'train_idx',find(train_idx)));
    [td_bl,bl_nnmodel] = getModel(td_bl,struct('model_type','nn',...
        'model_name','vel_from_spikes','out_signals',{{'vel'}},...
        'in_signals',{'S1_spikes'},'train_idx',find(train_idx)));

    % predict basline test set
    bl_lmvel_pred = cat(1,td_bl.linmodel_vel_from_spikes);
    bl_nnvel_pred = cat(1,td_bl.nn_vel_from_spikes);
    bl_dir_lmpred = atan2(bl_lmvel_pred(:,2),bl_lmvel_pred(:,1));
    bl_dir_nnpred = atan2(bl_nnvel_pred(:,2),bl_nnvel_pred(:,1));

    % plot predictions vs trial
    figure
    plot((1:sum(test_idx))',minusPi2Pi(bl_dir_lmpred(test_idx)-bl_dir(test_idx)),'o','linewidth',2);
    hold on
    plot([1 sum(test_idx)],[0 0],'--k','linewidth',2);

    % predict washout
    td_wo = getModel(td_wo,bl_linmodel);
    td_wo = getModel(td_wo,bl_nnmodel);
    wo_lmvel_pred = cat(1,td_wo.linmodel_vel_from_spikes);
    wo_nnvel_pred = cat(1,td_wo.nn_vel_from_spikes);
    wo_dir_lmpred = atan2(wo_lmvel_pred(:,2),wo_lmvel_pred(:,1));
    wo_dir_nnpred = atan2(wo_nnvel_pred(:,2),wo_nnvel_pred(:,1));

    % plot predictions vs trial
    figure
    plot((1:length(wo_dir))',minusPi2Pi(wo_dir_lmpred-wo_dir),'ro','linewidth',2);
    hold on
    plot([1 length(wo_dir)],[0 0],'--k','linewidth',2);

%% Fit GLM to just kinematics during bumps in baseline and evaluate R^2 in washout
    % Han_20171206
    %   - washout pr2 is pretty low...
    %   - units don't fire very much
    clearvars -except trial_data
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
    % Han_20171206
    clearvars -except trial_data
    % get only rewards
    [~,td] = getTDidx(trial_data,'result','R');
    td = removeBadTrials(td,struct('remove_nan_idx',false));
    
    % split into baseline and washout
    [~,td_ad] = getTDidx(td,'epoch','AD');
    
    % trim to just bumps
    td_ad = trimTD(td_ad,{'idx_goCueTime',0},{'idx_endTime',0});
    td_ad = binTD(td_ad,5);
    
    % fit glm to baseline
    [td_ad,ad_model] = getModel(td_ad,struct('model_type','glm',...
        'model_name','S1_handle','in_signals',{{'pos';'vel';'force'}},...
        'out_signals',{'S1_spikes'},'train_idx',(length(td_ad)-150):length(td_ad)));
    
    % evaluate washout period
    eval_params = ad_model;
    eval_params.eval_metric = 'pr2';
    eval_params.trial_idx = [1 length(td_ad)];
    full_eval = evalModel(td_ad,eval_params);
    eval_params.trial_idx = [length(td_ad)-150 length(td_ad)];
    part_eval = evalModel(td_ad,eval_params);
    eval_params.trial_idx = 1:30:length(td_ad);
    ad_eval = evalModel(td_ad,eval_params);
    
    % plot resulting evals
    figure
    trial_num = eval_params.trial_idx(1:end-1);
    for i=1:size(ad_eval,2)
        clf;
        single_eval = squeeze(ad_eval(:,i,:));
        patch([trial_num';flipud(trial_num')],[single_eval(:,1);flipud(single_eval(:,2))],[1,0,0]);
        hold on
        plot(trial_num,zeros(size(trial_num)),'--k','linewidth',2)
        waitforbuttonpress
    end

%% look at pR2 vs num spikes
