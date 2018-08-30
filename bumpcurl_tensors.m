%% test out tensors on bumpcurl data
    % First make tensor out of trial data
        [~,td] = getTDidx(trial_data,'result','R');
        td = removeBadTrials(td,struct('remove_nan_idx',false));
        
        td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
        td = smoothSignals(td,struct('signals','markers'));
        td = getDifferential(td,struct('signals','markers','alias','marker_vel'));
        % add firing rates rather than spike counts
        td = addFiringRates(td,struct('array','S1'));
        td = smoothSignals(td,struct('signals',{{'S1_spikes'}},'calc_rate',true,'kernel_SD',0.05));
        td = sqrtTransform(td,'S1_spikes');

        % Remove unsorted channels
        keepers = (td(1).S1_unit_guide(:,2)~=0);
        for trial = 1:length(td)
            td(trial).S1_unit_guide = td(trial).S1_unit_guide(keepers,:);
            td(trial).S1_spikes = td(trial).S1_spikes(:,keepers);
        end
        
        % trim to just movement
        % num_bins_before = floor(still_bins/2);
        num_bins_before = 30;
        num_bins_after = 60;
        
        td= trimTD(td,struct('idx_start',{{'idx_movement_on',-num_bins_before}},'idx_end',{{'idx_movement_on',num_bins_after}},'zero_pad',true));
        % clean nans out...?
        nanners = isnan(cat(1,td.target_direction));
        if any(nanners)
            error('there are nanners');
            td_act = td_act(~nanners);
        end

        % colormap
        cm_viridis = viridis(200);

    % make tensor
        num_neurons = length(td(1).S1_unit_guide);
        num_markers = size(td(1).markers,2);
        num_timepoints = num_bins_before+num_bins_after+1;
        num_trials = length(td);
        neural_data = tensor(cat(3,td.S1_spikes)); % dimensions are timepoints x neurons x trials
        behave_data = tensor(cat(2,cat(3,td.marker_vel),cat(3,td.markers))); % dimensions are timepoints x markers x trials

    % decompose
        num_factors = 10;
        M_behave = cp_als(behave_data,num_factors);
        M_neural = cp_als(neural_data,num_factors);

    % Look at temporal factors
        neural_temporal = M_neural.U{1};
        timevec = ((1:num_timepoints)-num_bins_before-1)*td(1).bin_size;
        figure
        for i = 1:num_factors
            subplot(5,2,i)
            plot(timevec,neural_temporal(:,i),'-k','linewidth',3)
            hold on
            plot(timevec([1 end]),[0 0],'-k','linewidth',2)
            plot([0 0],ylim,'--k','linewidth',2)
        end

    % Look at trial factors
        % color by direction
        num_colors = 4;
        dir_colors = linspecer(num_colors);
        trial_dirs = cat(1,td.target_direction);
        dir_idx = mod(round(trial_dirs/(2*pi/num_colors)),num_colors)+1;
        trial_colors = dir_colors(dir_idx,:);

        neural_trial = M_neural.U{3};
        trialvec = 1:num_trials;
        figure
        for i = 1:num_factors
            subplot(5,2,i)
            scatter(trialvec,neural_trial(:,i),[],trial_colors,'filled')
            hold on
            plot(trialvec([1 end]),[0 0],'-k','linewidth',2)
        end

    % Look at neural/marker factors
        neural_markers = M_neural.U{2};
        neuralvec = 1:num_neurons;
        figure
        for i = 1:num_factors
            subplot(5,2,i)
            bar(neuralvec,neural_markers(:,i))
        end
        

