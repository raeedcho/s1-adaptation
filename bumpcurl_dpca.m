%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bumpcurl_analysis - 
%%      script to run bumpcurl task analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up meta info
    if ispc
        dataroot = 'G:\raeed\project-data\limblab\s1-adapt';
    else
        dataroot = '/data/raeed/project-data/limblab/s1-adapt';
    end
    
    file_info = dir(fullfile(dataroot,'td-library','*CO*.mat'));
    filenames = horzcat({file_info.name})';

%% Loop through trial data files to clean them up
    trial_data_cell = cell(1,length(filenames));
    for filenum = 1:length(filenames)
        %% load and preprocess data
        td = load(fullfile(dataroot,'td-library',[filenames{filenum}]));
    
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

%% Plot trial info (hand speed and example rasters)
    for filenum = 1:length(trial_data_cell)
        %% load and preprocess data
        td = trial_data_cell{filenum};

        % trim from go cue to end time (skip bump)
        td = trimTD(td,{'idx_goCueTime',0},{'idx_endTime',0});
    
        bl_idx = getTDidx(td,'epoch','BL');
        ad_idx = getTDidx(td,'epoch','AD');
        wo_idx = getTDidx(td,'epoch','WO');
        
        metric = getLearningMetrics(td,struct(...
            'which_metric','angle',...
            'use_bl_ref',true,...
            'fit_bl_ref_curve',false,...
            'vel_or_pos','pos',...
            'target_dir_fieldname','tgtDir',...
            'time_window',{{'idx_movement_on',0;'idx_movement_on',40}}));

        % plot metrics
        smoothing_window_size = 5;
        figure('defaultaxesfontsize',18)
        plot(conv(metric,ones(1,smoothing_window_size)/smoothing_window_size,'same'),'k-')
        hold on
        plot(repmat(ad_idx(1),2,1),[-1 1],'k--','linewidth',3)
        plot(repmat(wo_idx(1),2,1),[-1 1],'k--','linewidth',3)
        plot([0 wo_idx(end)],[0 0],'k-','linewidth',2)

        % plot trials
        figure('defaultaxesfontsize',18)
        num_trials_to_plot = 60;
        % baseline
        subplot(3,3,4)
        for trialnum = getTDidx(td,'epoch','BL','rand',num_trials_to_plot)
            plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
            hold on
        end
        axis equal
        title 'Baseline'

        % adaptation
        subplot(3,3,2)
        for trialnum = getTDidx(td,'epoch','AD','rand',num_trials_to_plot,'range',[0 0.3])
            plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
            hold on
        end
        axis equal
        title 'Adaptation'
        subplot(3,3,5)
        for trialnum = getTDidx(td,'epoch','AD','rand',num_trials_to_plot,'range',[0.3 0.8])
            plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
            hold on
        end
        axis equal
        subplot(3,3,8)
        for trialnum = getTDidx(td,'epoch','AD','rand',num_trials_to_plot,'range',[0.8 1])
            plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
            hold on
        end
        axis equal

        % washout
        subplot(3,3,3)
        for trialnum = getTDidx(td,'epoch','WO','rand',num_trials_to_plot,'range',[0 0.3])
            plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
            hold on
        end
        axis equal
        title 'Adaptation'
        subplot(3,3,6)
        for trialnum = getTDidx(td,'epoch','WO','rand',num_trials_to_plot,'range',[0.3 0.8])
            plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
            hold on
        end
        axis equal
        subplot(3,3,9)
        for trialnum = getTDidx(td,'epoch','WO','rand',num_trials_to_plot,'range',[0.8 1])
            plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
            hold on
        end
        axis equal
    end

%% try dPCA on data
    for filenum = 1:length(trial_data_cell)
        %% load and preprocess data
        td = trial_data_cell{filenum};

        % trim from go cue to end time (skip bump)
        spikes_in_td = getTDfields(td,'spikes');
        td = smoothSignals(td,struct('signals',{spikes_in_td}));
        if td(1).bin_size == 0.005
            td = binTD(td,2);
        end
        td = trimTD(td,{'idx_movement_on',-10},{'idx_movement_on',50});
    
        % assign target direction blocks
        tgt_dirs = cat(2,td.tgtDir);
        unique_tgt_dirs = unique(tgt_dirs);
        % assume 16 targets and reassign to one of 4
        switch length(unique_tgt_dirs)
            case 16
                tgt_block_assign = reshape(repmat([1 2 3 4],4,1),1,16);
            case 8
                tgt_block_assign = reshape(repmat([1 2 3 4],2,1),1,8);
            case 4
                tgt_block_assign = [1 2 3 4];
        end
        for dirnum = 1:length(unique_tgt_dirs)
            trial_idx = getTDidx(td,'tgtDir',unique_tgt_dirs(dirnum));
            [td(trial_idx).tgt_dir_block] = deal(tgt_block_assign(dirnum));
        end

        bl_idx = getTDidx(td,'epoch','BL');
        ad_idx = getTDidx(td,'epoch','AD');
        wo_idx = getTDidx(td,'epoch','WO');

        td_dpca = td;

        % get dPCA conditions
        learning_blocks = {...
            getTDidx(td_dpca,'epoch','BL','range',[0.5 1]),...
            getTDidx(td_dpca,'epoch','AD','range',[0 0.2]),...
            getTDidx(td_dpca,'epoch','AD','range',[0.2 0.5]),...
            getTDidx(td_dpca,'epoch','AD','range',[0.5 1]),...
            getTDidx(td_dpca,'epoch','WO','range',[0 0.3]),...
            };

        % get actual dPCA
        marg_colors = [150 150 150; 23 100 171; 187 20 25; 114 97 171]/256; % grey, blue, red, purple
        for arraynum = 1:length(spikes_in_td)
            % td_dpca = runDPCA(td_dpca,'tgt_dir_block',learning_blocks,struct(...
            %     'signals',spikes_in_td(arraynum),...
            %     'marg_names',{{'time','direction','learning'}},...
            %     'combined_params',{{{1},{2,[1 2]},{3,[1 3],[2 3],[1 2 3]}}},...
            %     'marg_colors',marg_colors(1:3,:),...
            %     'do_plot',true,'num_dims',10));
            td_dpca = runDPCA(td_dpca,'tgt_dir_block',learning_blocks,struct(...
                'signals',spikes_in_td(arraynum),...
                'do_plot',true,'num_dims',10));
            saveas(gcf,sprintf('%s_%s_%s_dPCA.png',td(1).monkey,strrep(td(1).date_time,'/','-'),spikes_in_td{arraynum}))
        end
    end

