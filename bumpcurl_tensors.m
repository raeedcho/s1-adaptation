%% test out tensors on bumpcurl data
    % First make tensor out of trial data
        [~,td] = getTDidx(trial_data,'result','R');
        td = removeBadTrials(td,struct('remove_nan_idx',false));
        
        td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
        td = smoothSignals(td,struct('signals','markers'));
        td = getDifferential(td,struct('signals','markers','alias','marker_vel'));
        % add firing rates rather than spike counts
        % td = smoothSignals(td,struct('signals',{{'S1_spikes'}},'calc_rate',true,'kernel_SD',0.05));
        % td = sqrtTransform(td,'S1_spikes');

        % Remove unsorted channels
        keepers = (td(1).S1_unit_guide(:,2)~=0);
        for trial = 1:length(td)
            td(trial).S1_unit_guide = td(trial).S1_unit_guide(keepers,:);
            td(trial).S1_spikes = td(trial).S1_spikes(:,keepers);
        end
        
        % trim to just movement
        % num_bins_before = floor(still_bins/2);
        num_bins_before = 30;
        num_bins_after = 100;
        
        td= trimTD(td,struct('idx_start',{{'idx_movement_on',-num_bins_before}},'idx_end',{{'idx_movement_on',num_bins_after}},'zero_pad',true));
        % clean nans out...?
        nanners = isnan(cat(1,td.target_direction));
        if any(nanners)
            error('there are nanners');
            td_act = td_act(~nanners);
        end

        % colormap
        cm_viridis = viridis(200);

    % make tensors and decompose
        num_neurons = length(td(1).S1_unit_guide);
        num_markers = size(td(1).markers,2);
        num_timepoints = num_bins_before+num_bins_after+1;
        num_trials = length(td);
        num_factors = 10;

        if isfield(td,'emg')
            emg_data = cat(3,td.emg); % dimensions are timepoints x emg x trials
            emg_data = emg_data-repmat(mean(mean(emg_data,3),1),num_timepoints,1,num_trials);
            emg_tensor = tensor(emg_data);
            M_emg = cp_als(emg_tensor,num_factors,'maxiters',100,'printitn',10);
        end
        
        behave_data = cat(2,cat(3,td.marker_vel),cat(3,td.markers)); % dimensions are timepoints x markers x trials
        behave_data = behave_data-repmat(mean(mean(behave_data,3),1),num_timepoints,1,num_trials);
        behave_tensor = tensor(behave_data);
        M_behave = cp_als(behave_tensor,num_factors,'maxiters',200,'printitn',10);

        neural_data = cat(3,td.S1_spikes); % dimensions are timepoints x neurons x trials
        neural_tensor = tensor(neural_data);
        tic;
        M_neural = cp_apr(neural_tensor,num_factors,'maxiters',100,'printitn',10);
        toc

    % color by direction
        num_colors = 4;
        dir_colors = linspecer(num_colors);
        trial_dirs = cat(1,td.target_direction);
        dir_idx = mod(round(trial_dirs/(2*pi/num_colors)),num_colors)+1;
        trial_colors = dir_colors(dir_idx,:);

    % plot tensor decomposition
        plotTensorDecomp(M_neural,struct('trial_colors',trial_colors,'bin_size',td(1).bin_size,'temporal_zero',num_bins_before+1))
