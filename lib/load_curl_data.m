function trial_data_cell = load_curl_data(filenames)
% loads all curl field data from folder, conditions it, and puts it into a
% cell array of TrialData structures

    trial_data_cell = cell(1,length(filenames));
    for filenum = 1:length(filenames)
        %% load and preprocess data
        td = load([filenames{filenum}]);
    
        % rename trial_data for ease
        td = td.trial_data;

        % make Kyle's and Matt's data consistent with my conventions
        td = aliasTDfields(td,struct('alias_list',{{...
            'date','date_time';...
            'trial_type','epoch';...
            'S1_spikes_bins','binned_S1';...
            'target_direction','tgtDir';...
            'trial_id','trialID';...
            'idx_trial_start','idx_startTime';...
            'idx_trial_end','idx_endTime';...
            'idx_go_cue','idx_goCueTime';...
            }}));

        % process marker data (if it's there)
        if isfield(td,'markers')
            % find times when markers are NaN and replace with zeros temporarily
            for trialnum = 1:length(td)
                markernans = isnan(td(trialnum).markers);
                td(trialnum).markers(markernans) = 0;
                td(trialnum) = smoothSignals(td(trialnum),struct('signals','markers'));
                td(trialnum).markers(markernans) = NaN;
                clear markernans
            end
    
            % get marker velocity
            td = getDifferential(td,struct('signals','markers','alias','marker_vel'));
        end
        
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
        fprintf('Removed %d trials because of missing target direction\n',sum(nanners))
        % biggers = ~isnan(cat(1,td.bumpDir)) & abs(cat(1,td.bumpDir))>360;
        % td = td(~biggers);
        % fprintf('Removed %d trials because bump direction makes no sense\n',sum(biggers))
    
        % remove trials where markers aren't present
        if isfield(td,'markers')
            bad_trial = false(length(td),1);
            for trialnum = 1:length(td)
                if any(any(isnan(td(trialnum).markers)))
                    bad_trial(trialnum) = true;
                end
            end
            td(bad_trial) = [];
            fprintf('Removed %d trials because of missing markers\n',sum(bad_trial))
        end
        
        % remove trials where muscles aren't present
        if isfield(td,'muscle_len')
            bad_trial = false(length(td),1);
            for trialnum = 1:length(td)
                if any(any(isnan(td(trialnum).muscle_len) | isnan(td(trialnum).muscle_vel)))
                    bad_trial(trialnum) = true;
                end
            end
            td(bad_trial) = [];
            fprintf('Removed %d trials because of missing muscles\n',sum(bad_trial))
        end
        
        % find the relevant movmement onsets
        if ~isfield(td,'idx_movement_on')
            td = getMoveOnsetAndPeak(td,struct(...
                'start_idx','idx_goCueTime',...
                'start_idx_offset',15,...
                'peak_idx_offset',20,...
                'end_idx','idx_endTime',...
                'method','peak',...
                'peak_divisor',10,...
                'min_ds',1));
        end

        % remove low firing neurons
        td = removeBadNeurons(td,struct(...
            'min_fr',1,...
            'fr_window',{{'idx_movement_on',0;'idx_movement_on',11}},...
            'calc_fr',true));
    
        trial_data_cell{filenum} = td;
    end