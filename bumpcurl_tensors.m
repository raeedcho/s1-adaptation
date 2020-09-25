%% Script to run tensor decomposition on curl field data

%% Set up meta info
    if ispc
        % dataroot = 'G:\raeed\project-data';
        dataroot = 'C:\Users\Raeed\data\project-data';
    else
        dataroot = '/data/raeed/project-data';
    end
    
    file_info = dir(fullfile(dataroot,'limblab','s1-adapt','td-library','*CO*.mat'));
    filenames = horzcat({file_info.name})';

%% Loop through trial data files to clean them up
    trial_data_cell = load_curl_data(fullfile(dataroot,'limblab','s1-adapt','td-library',filenames));

%% test out tensors on bumpcurl data
    trim_start = -0.3;
    trim_end = 0.5;
    num_factors = 5;
    for filenum = 1:length(trial_data_cell)
        td = trial_data_cell{filenum};
        
        % trim from go cue to end time (skip bump)
        spikes_in_td = getTDfields(td,'spikes');
        td = smoothSignals(td,struct('signals',{spikes_in_td},'width',0.1));
        if td(1).bin_size == 0.005
            td = binTD(td,2);
        end
        td = trimTD(td,struct(...
            'idx_start',{{'idx_movement_on',trim_start/td(1).bin_size}},...
            'idx_end',{{'idx_movement_on',trim_end/td(1).bin_size}},...
            'remove_short',true));
        
        timevec = trim_start:td(1).bin_size:(length(td(1).vel)*td(1).bin_size+trim_start-0.001);
        for trialnum = 1:length(td)
            td(trialnum).timevec = timevec';
        end
    
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
            [td(trial_idx).target_block] = deal(tgt_block_assign(dirnum));
        end
        
        cm_viridis = viridis(200);
        
        % subselect baseline and adaptation trials
        ad_offset = 100;
        td_idx = cat(2,...
            getTDidx(td,'epoch','BL','range',[1 ad_offset]),...
            getTDidx(td,'epoch','AD','range',[1 200]));
        td_select = td(td_idx);

        for arraynum = 1:length(spikes_in_td)
            % make tensors and decompose
%             num_neurons = length(td(1).S1_unit_guide);
%             num_markers = size(td(1).markers,2);
%             num_timepoints = num_bins_before+num_bins_after+1;
%             num_trials = length(td);

%             if isfield(td,'emg')
%                 emg_data = cat(3,td.emg); % dimensions are timepoints x emg x trials
%                 emg_data = emg_data-repmat(mean(mean(emg_data,3),1),num_timepoints,1,num_trials);
%                 emg_tensor = tensor(emg_data);
%                 M_emg = cp_als(emg_tensor,num_factors,'maxiters',100,'printitn',10);
%             end
% 
%             if isfield(td,'markers')
%                 behave_data = cat(2,cat(3,td.marker_vel),cat(3,td.markers)); % dimensions are timepoints x markers x trials
%                 behave_data = behave_data-repmat(mean(mean(behave_data,3),1),num_timepoints,1,num_trials);
%                 behave_tensor = tensor(behave_data);
%                 M_behave = cp_als(behave_tensor,num_factors,'maxiters',200,'printitn',10);
%             end

            target_blocks = unique([td.target_block]);
            % want neural data to be timepoints x neurons x trials x target block
            neural_data = zeros(...
                size(td(1).(spikes_in_td{arraynum}),1),...
                size(td(1).(spikes_in_td{arraynum}),2),...
                55,...
                length(target_blocks));
                
            for blocknum = 1:length(target_blocks)
                [~,td_bl_dir] = getTDidx(td,'target_block',target_blocks(blocknum),'epoch','BL','range',[1 25]);
                [~,td_ad_dir] = getTDidx(td,'target_block',target_blocks(blocknum),'epoch','AD','range',[1 30]);
                neural_data(:,:,1:25,blocknum) = cat(3,td_bl_dir.(spikes_in_td{arraynum}));
                neural_data(:,:,26:end,blocknum) = cat(3,td_ad_dir.(spikes_in_td{arraynum}));
            end
            neural_data = cat(3,td_select.(spikes_in_td{arraynum})); % dimensions are timepoints x neurons x trials
            neural_tensor = tensor(neural_data);
            tic;
%             M_neural = cp_apr(neural_tensor,num_factors,'maxiters',100,'printitn',10);
            M_neural_als = cp_als(neural_tensor,num_factors,'maxiters',200,'printitn',10);
            toc

            % color by direction
            dir_colors = linspecer(length(unique([td_select.target_block])));
            trial_colors = dir_colors(cat(1,td_select.target_block),:);

            % plot tensor decomposition
%             plotTensorDecomp(M_neural_als,struct('trial_colors','k','timevec',td_select(1).timevec))
            M = M_neural_als;
            figure('position',[50 50 800 1200])
            ncol = length(M.U);
            nrow = size(M.U{1},2);

            % Look at signal factors
            signal_factors = M.U{2};
            markervec = 1:size(signal_factors,1);
            for i = 1:size(signal_factors,2)
                subplot(nrow,ncol,(i-1)*ncol+1)
                bar(markervec,signal_factors(:,i))
                set(gca,'box','off','tickdir','out')
            end

            % Plot temporal factors
            temporal_factors = M.U{1};
            if isempty(timevec)
                timevec = 1:size(temporal_factors,1);
            end
            for i = 1:size(temporal_factors,2)
                subplot(nrow,ncol,(i-1)*ncol+2)
                plot(timevec,temporal_factors(:,i),'-k','linewidth',3)
                hold on
                plot(timevec([1 end]),[0 0],'-k','linewidth',2)
                plot([0 0],ylim,'--k','linewidth',2)
                set(gca,'box','off','tickdir','out')
            end

            % Look at trial factors
            trial_factors = M.U{3};
            trialvec = (1:size(trial_factors,1))-ad_offset;
            for i = 1:size(trial_factors,2)
                subplot(nrow,ncol,(i-1)*ncol+3)
                plot(trialvec([1 end]),[0 0],'-k','linewidth',2)
                hold on
                scatter(trialvec,trial_factors(:,i),10,trial_colors,'filled')
                plot(trialvec,smooth_data(trial_factors(:,i),1,15),'k')
                plot([0 0],0.15*[-1 1],'-r')
                set(gca,'box','off','tickdir','out')
            end
            
%             % target block factors
%             target_factors = M.U{4};
%             targetvec = 1:size(target_factors);
%             for i = 1:size(target_factors,2)
%                 subplot(nrow,ncol,(i-1)*ncol+4)
%                 plot(targetvec([1 end]),[0 0],'-k','linewidth',2)
%                 hold on
%                 scatter(targetvec,target_factors(:,i),[],dir_colors,'filled')
%                 set(gca,'box','off','tickdir','out')
%             end
            
            suptitle(sprintf('%s %s',td(1).monkey,strrep(spikes_in_td{arraynum},'_spikes','')))
        end
    end
