function trial_data_cell = load_curl_data(filenames)
% loads all curl field data from folder, conditions it, and puts it into a
% cell array of TrialData structures

    tgt_block_assign = [0 pi/2 pi 3*pi/2];
    learning_block_ranges = 0:0.2:1;

    trial_data_cell = cell(1,length(filenames));
    for filenum = 1:length(filenames)
        %% load and preprocess data
        td = load([filenames{filenum}]);
    
        % rename trial_data for ease
        td = td.trial_data;

        % because Duncan's S1 is already smoothed for some reason
        if isfield(td,'S1_spikes_bins')
            td = aliasTDfields(td,struct('alias_list',{{...
                'S1_spikes','S1_smooth'}}));
        end
        % make Kyle's and Matt's data consistent with my conventions
        td = aliasTDfields(td,struct('overwrite_fields',true,'alias_list',{{...
            'date','date_time';...
            'trial_type','epoch';...
            'S1_spikes_bins','S1_spikes';...
            'target_direction','tgtDir';...
            'trial_id','trialID';...
            'idx_trial_start','idx_startTime';...
            'idx_trial_end','idx_endTime';...
            'idx_go_cue','idx_goCueTime';...
            }}));

        
        % get speed and ds
        td = getNorm(td,struct('signals','vel','field_extra','_norm'));
        td = getDifferential(td,struct('signals','vel_norm','alias','dvel_norm'));
        
        % remove unsorted neurons
        arrays_in_td = getTDfields(td,'arrays');
        for arraynum = 1:length(arrays_in_td)
            arrayname = arrays_in_td{arraynum};
            for trialnum = 1:length(td)
                unit_ids = td(trialnum).([arrayname '_unit_guide']);
                unsorted_units = (unit_ids(:,2)==0);
                new_unit_guide = unit_ids(~unsorted_units,:);

                td(trialnum).(sprintf('%s_unit_guide',arrayname)) = new_unit_guide;
                
                spikes = td(trialnum).(sprintf('%s_spikes',arrayname));
                spikes(:,unsorted_units) = [];
                td(trialnum).(sprintf('%s_spikes',arrayname)) = spikes;
            end
        end

        % check to make sure unit guides are consistent
        if length(td)>1
            bad_units = checkUnitGuides({td.([arrayname '_unit_guide'])});
            assert(isempty(bad_units), 'Unit guides have different units in them!')
        end
    
        % prep trial data by getting only rewards
        % split into trials
        td = splitTD(...
            td,...
            struct(...
                'split_idx_name','idx_startTime',...
                'linked_fields',{{...
                    'trialID',...
                    'result',...
                    'tgtDir',...
                    }},...
                'start_name','idx_startTime',...
                'end_name','idx_endTime'));
        reward_idx = [td.result]=='R';
        incomplete_idx = [td.result]=='I';
        % td = td(reward_idx | incomplete_idx);
        td = td(reward_idx);
        td = reorderTDfields(td);
        
        % clean nans out...?
        nanners = isnan(cat(1,td.tgtDir));
        td = td(~nanners);
        if any(nanners)
            fprintf('Removed %d trials because of missing target direction\n',sum(nanners))
        end
        
        % find the relevant movmement onsets
        if ~isfield(td,'idx_movement_on')
            % kyle's params
            td = getMoveOnsetAndPeak(td,struct(...
                'start_idx','idx_goCueTime',...
                'end_idx','idx_endTime',...
                'method','thresh',...
                'min_ds',2.5));
        end
        
        % get rid of trials that for some reason have a movement onset after peak speed
        bad_trial = vertcat(td.idx_peak_speed)<vertcat(td.idx_movement_on);
        td(bad_trial) = [];
        if any(bad_trial)
            fprintf('Removed %d trials because of peak speed happening before movement onset\n',sum(bad_trial))
        end
        
        % remove low firing neurons
        td = removeBadNeurons(td,struct(...
            'min_fr',1,...
            'fr_window',{{'idx_movement_on',0;'idx_movement_on',11}},...
            'calc_fr',true));
        
        % Switch to radians instead of degrees
        tgt_dirs = cat(2,td.tgtDir);
        unique_tgt_dirs = sort(unique(tgt_dirs));
        if any(unique_tgt_dirs>4)
            % we have degrees instead of rad
            for trialnum = 1:length(td)
                td(trialnum).tgtDir = td(trialnum).tgtDir*pi/180;
            end
        end
        
        % assign targets to blocks for cardinal directions
        for trialnum = 1:length(td)
            block_dist = angleDiff(td(trialnum).tgtDir,tgt_block_assign);
            [~,blocknum] = min(abs(block_dist));
            td(trialnum).target_block = tgt_block_assign(blocknum);
        end
        
        % assign learning blocks
        learning_blocks = cell(1,length(learning_block_ranges)-1);
        for blocknum = 1:length(learning_blocks)
            learning_blocks{blocknum} = getTDidx(td,'epoch','AD','range',learning_block_ranges([blocknum blocknum+1]));
            [td(learning_blocks{blocknum}).learning_block] = deal(blocknum);
        end

        % add fake learning block indices to baseline and washout to match learning blocks in adaptation
        num_trials_in_learning_block = mode(cellfun(@length,learning_blocks));
        bl_idx = getTDidx(td,'epoch','BL');
        wo_idx = getTDidx(td,'epoch','WO');
        %reference to last baseline idx
        bl_blocknum = num2cell(floor((bl_idx-(bl_idx(end)+1))/num_trials_in_learning_block)+1);
        [td(bl_idx).learning_block] = deal(bl_blocknum{:});
        wo_blocknum = num2cell(floor((wo_idx-(wo_idx(1)))/num_trials_in_learning_block)+length(learning_blocks)+1);
        [td(wo_idx).learning_block] = deal(wo_blocknum{:});
        
        % add learning metric to trial_data
        bin_size = td(1).bin_size;
        td_temp=td;
%         td_temp = smoothSignals(td,struct('signals','vel','width',0.1));
        % kyle's params
        metric = getLearningMetrics(td_temp,struct(...
            'which_metric','angle',...
            'use_bl_ref',true,...
            'fit_bl_ref_curve',true,...
            'vel_or_pos','vel',...
            'target_dir_fieldname','tgtDir',...
            'time_window',{{'idx_movement_on',0;'idx_movement_on',ceil(30e-3/bin_size)}}));
        
        metric_cell = num2cell(metric);
        [td.learning_metric] = deal(metric_cell{:});
        
        % bin signals if needed
        if td(1).bin_size == 0.005
            td = binTD(td,2);
        end
    
        trial_data_cell{filenum} = td;
    end