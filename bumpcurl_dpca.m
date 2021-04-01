function bumpcurl_dpca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bumpcurl_dpca - 
%%      function to run bumpcurl dpca analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up meta info
    if ispc
        % dataroot = 'G:\raeed\project-data';
        dataroot = 'C:\Users\Raeed\data\project-data';
    else
        % dataroot = '/data/raeed/project-data';
        dataroot = '/home/raeed/data/project-data';
    end
    
    file_info = dir(fullfile(dataroot,'limblab','s1-adapt','td-library','*CO*.mat'));
    filenames = horzcat({file_info.name})';

%% Loop through trial data files to clean them up
    trial_data_cell = load_curl_data(fullfile(dataroot,'limblab','s1-adapt','td-library',filenames));

%% Plot trial info (hand speed and example rasters)
    plot_learning_behavior(trial_data_cell);
    
%% Calculate variance due to learning using bootstrapping over trials
    num_learning_groups = 5;
    do_shuffle_control = false;
    num_control_shuffles = 1;
    [margvar_cell,margvar_bl_cell,learning_metric_cell] = deal(cell(length(trial_data_cell),1));
    [margvar_learning_shuffle_cell,margvar_target_shuffle_cell] = deal(cell(length(trial_data_cell),num_control_shuffles));
    for filenum = 1:length(filenames)
        td = trial_data_cell{filenum};
        [~,td_bl] = getTDidx(td,'epoch','BL');
        [~,td_ad] = getTDidx(td,'epoch','AD');
        
        num_trials = min(length(td_bl),length(td_ad));
        
        % get baseline values for metrics
        % bl_trial_frac = ad_trial_frac*length(td_ad)/length(td_bl);
        margvar_bl_cell{filenum} = get_dpca_var(td_bl,struct(...
            'num_boots',1,...
            'plot_marg',false,...
            'plot_comp_dpcs',false,...
            'learning_block_ranges',linspace(0,num_trials/length(td_bl),num_learning_groups+1)));
        
        margvar_cell{filenum} = get_dpca_var(td_ad,struct(...
            'num_boots',1,...
            'plot_marg',false,...
            'plot_comp_dpcs',false,...
            'comp_to_plot','learning',...
            'learning_block_ranges',linspace(0,num_trials/length(td_ad),num_learning_groups+1)));
        % [margvar_cell{filenum},learning_metric_cell{filenum}] = get_dpca_var(td_ad,struct(...
        %     'num_boots',1,...
        %     'plot_marg',false,...
        %     'plot_comp_dpcs',false,...
        %     'learning_block_ranges',0:ad_trial_frac:1));

        % get trial-shuffled values for metrics
        if do_shuffle_control
            for shufflenum = 1:num_control_shuffles
                try
                    margvar_learning_shuffle_cell{filenum,shufflenum} = get_dpca_var(td_ad,struct(...
                        'num_boots',1,...
                        'plot_marg',false,...
                        'plot_comp_dpcs',false,...
                        'comp_to_plot','learning',...
                        'shuffle_conds','learning_block',...
                        'learning_block_ranges',0:ad_trial_frac:1));
                catch ME
                    continue
                end
    %             margvar_target_shuffle_cell{filenum,shufflenum} = get_dpca_var(td_ad,struct(...
    %                 'num_boots',1,...
    %                 'plot_marg',false,...
    %                 'plot_comp_dpcs',false,...
    %                 'comp_to_plot','target',...
    %                 'shuffle_conds','target_block',...
    %                 'learning_block_ranges',0:ad_trial_frac:1));

                shuffleID = table(shufflenum,'VariableNames',{'shuffleID'});
                shuffleID.Properties.VariableDescriptions = {'meta'};
                margvar_learning_shuffle_cell{filenum,shufflenum} = horzcat(...
                    margvar_learning_shuffle_cell{filenum,shufflenum},...
                    repmat(shuffleID,height(margvar_learning_shuffle_cell{filenum,shufflenum}),1));
    %             margvar_target_shuffle_cell{filenum,shufflenum} = horzcat(...
    %                 margvar_target_shuffle_cell{filenum,shufflenum},...
    %                 repmat(shuffleID,height(margvar_target_shuffle_cell{filenum,shufflenum}),1));

                if mod(shufflenum,floor(num_control_shuffles*0.1))==0
                    fprintf('\tFinished shuffle control %d of %d\n',shufflenum,num_control_shuffles)
                end
            end
        end

        fprintf('Finished file %d of %d\n',filenum,length(filenames))

    end
    margvar_table = vertcat(margvar_cell{:});
    margvar_bl_table = vertcat(margvar_bl_cell{:});
    margvar_learning_shuffle_table = vertcat(margvar_learning_shuffle_cell{:});
%     margvar_target_shuffle_table = vertcat(margvar_target_shuffle_cell{:});
    % learning_metric_table = vertcat(learning_metric_cell{:});
    
%% compare learning metric with dPC movement
    % arrays = {'PMd_spikes','M1_spikes','S1_spikes'};
    % for arraynum = 1:length(arrays)
    %     [~,array_metric_table] = getNTidx(learning_metric_table,'array',arrays{arraynum});
    %     sessions = unique(array_metric_table.date_time);
    %     session_colors = linspecer(length(sessions));
    %     figure('defaultaxesfontsize',18)
    %     for sessionnum = 1:length(sessions)
    %         [~,session_table] = getNTidx(array_metric_table,'date_time',sessions{sessionnum});
    %         hold on
    %         scatter(session_table.rel_learning_metric,session_table.learning_dPC_dist_norm,[],session_colors(sessionnum,:),'filled')
    %     end
    %     set(gca,'box','off','tickdir','out','ylim',[0 0.75],'xlim',[-0.35 0.35])
    %     title(arrays{arraynum})
    %     xlabel('Learning metric diff')
    %     ylabel('Norm learning dPC diff')
    % end

%% compare dpca learning total margvar between baseline and adaptation (run marginalization in baseline)
    
%     arrays = {'S1_spikes','PMd_spikes','M1_spikes'};
    arrays = {'M1_spikes','S1_spikes'};
    arraynames = strrep(arrays,'_spikes','');
    monkeys = {'Chewie','Mihili','MrT','Duncan','Han'};
    monkey_colors = linspecer(length(monkeys));
    marker_size = 15;
    offset_size = 0.2;
    margs_to_plot = {'time','target','learning'};
    figure(...
        'defaultaxesfontsize',10,...
        'units','inches',...
        'position',[10 10 2*length(margs_to_plot) 3],...
        'papersize',[2*length(margs_to_plot) 3])
    for arraynum = 1:length(arrays)
        [~,array_metric_table] = getNTidx(margvar_table,'array',arrays{arraynum});
        sessions = unique(array_metric_table.date_time);
        for sessionnum = 1:length(sessions)
            [~,session_table] = getNTidx(array_metric_table,'date_time',sessions{sessionnum});
            if do_shuffle_control
                [~,learning_shuffle_table] = getNTidx(margvar_learning_shuffle_table,...
                    'array',arrays{arraynum},...
                    'date_time',sessions{sessionnum});
            end
            bl_tab_idx = strncmp(margvar_bl_table.date_time,sessions{sessionnum},10) & ...
                strcmpi(margvar_bl_table.array,arrays{arraynum});
            bl_table = margvar_bl_table(bl_tab_idx,:);
            assert(height(session_table)==1) % assume all sessions are on different days
            point_color = monkey_colors(strcmpi(session_table.monkey,monkeys),:);

            for margnum = 1:length(margs_to_plot)
                marg_idx = strcmpi(margs_to_plot{margnum},session_table.marg_names);

                subplot(1,length(margs_to_plot),margnum)
                hold on
%                 if ad_val>prctile(learning_shuffle_val,0.95)
%                     scatter(...
%                         arraynum+0.1,...
%                         ad_val,...
%                         marker_size,point_color,'filled')
%                 else
%                     scatter(...
%                         arraynum+0.1,...
%                         ad_val,...
%                         marker_size,point_color)
%                 end

                ad_val = session_table.marg_var_pct(marg_idx);
                scatter(...
                    arraynum+offset_size,...
                    ad_val,...
                    marker_size,point_color,'filled')
                if do_shuffle_control
                    learning_shuffle_val =(learning_shuffle_table.marg_var_pct(:,marg_idx));
                    scatter(...
                        repmat(arraynum-offset_size,num_control_shuffles,1),...
                        learning_shuffle_val,...
                        marker_size,point_color)
                    scatter(...
                        arraynum-offset_size,...
                        mean(learning_shuffle_val),...
                        marker_size,point_color)
                    plot(...
                        arraynum+[1 -1]*offset_size,...
                        horzcat(ad_val, mean(learning_shuffle_val)),...
                        '-','color',point_color)
                else
                    bl_val = bl_table.marg_var_pct(marg_idx);
                    scatter(...
                        arraynum-offset_size,...
                        bl_val,...
                        marker_size,point_color)
                    plot(...
                        arraynum+[1 -1]*offset_size,...
                        horzcat(ad_val, bl_val),...
                        '-','color',point_color)
                end
            end
        end
    end
    for margnum = 1:length(margs_to_plot)
        subplot(1,length(margs_to_plot),margnum)
        
        if margnum == 1
            ylabel('% of total neural variance')
            ytick_labels = 0:20:100;
        else
            ytick_labels = 0:20:100;
        end
        set(gca,...
            'box','off',...
            'tickdir','out',...
            'ylim',[0 100],'xlim',[0 length(arrays)+1],...
            'xtick',1:length(arrays),'xticklabel',arraynames,...
            'ytick',0:20:100,'yticklabel',ytick_labels,...
            'ticklabelinterpreter','none')
        
        title(margs_to_plot{margnum})
    end
    if do_shuffle_control
        suptitle('Condition-related variance in adaptation vs learning-block-shuffle control')
    else
        suptitle('Comparison between variance in baseline and adaptation')
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_learning_behavior(trial_data_cell)
    for filenum = 1:length(trial_data_cell)
        % load and preprocess data
        td = trial_data_cell{filenum};

        % trim from go cue to end time (skip bump)
        td = trimTD(td,{'idx_goCueTime',0},{'idx_endTime',0});
    
        % bl_idx = getTDidx(td,'epoch','BL');
        ad_idx = getTDidx(td,'epoch','AD');
        wo_idx = getTDidx(td,'epoch','WO');
        
        metric = cat(1,td.learning_metric);

        % plot metrics
        smoothing_window_size = 12;
        figure('defaultaxesfontsize',18)
        plot(conv(metric,ones(1,smoothing_window_size)/smoothing_window_size,'same'),'k-')
        hold on
        plot(repmat(ad_idx(1),2,1),[-1 1],'k--','linewidth',3)
        plot(repmat(wo_idx(1),2,1),[-1 1],'k--','linewidth',3)
        plot([0 wo_idx(end)],[0 0],'k-','linewidth',2)
        title(sprintf('%s %s',td(1).monkey,td(1).date_time))

        % plot trials
        % figure('defaultaxesfontsize',18)
        % num_trials_to_plot = 60;
        % % baseline
        % subplot(3,3,4)
        % for trialnum = getTDidx(td,'epoch','BL','rand',num_trials_to_plot)
        %     plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
        %     hold on
        % end
        % axis equal
        % title 'Baseline'

        % % adaptation
        % subplot(3,3,2)
        % for trialnum = getTDidx(td,'epoch','AD','rand',num_trials_to_plot,'range',[0 0.3])
        %     plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
        %     hold on
        % end
        % axis equal
        % title 'Adaptation'
        % subplot(3,3,5)
        % for trialnum = getTDidx(td,'epoch','AD','rand',num_trials_to_plot,'range',[0.3 0.8])
        %     plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
        %     hold on
        % end
        % axis equal
        % subplot(3,3,8)
        % for trialnum = getTDidx(td,'epoch','AD','rand',num_trials_to_plot,'range',[0.8 1])
        %     plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
        %     hold on
        % end
        % axis equal

        % % washout
        % subplot(3,3,3)
        % for trialnum = getTDidx(td,'epoch','WO','rand',num_trials_to_plot,'range',[0 0.3])
        %     plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
        %     hold on
        % end
        % axis equal
        % title 'Washout'
        % subplot(3,3,6)
        % for trialnum = getTDidx(td,'epoch','WO','rand',num_trials_to_plot,'range',[0.3 0.8])
        %     plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
        %     hold on
        % end
        % axis equal
        % subplot(3,3,9)
        % for trialnum = getTDidx(td,'epoch','WO','rand',num_trials_to_plot,'range',[0.8 1])
        %     plot(td(trialnum).pos(:,1),td(trialnum).pos(:,2),'k','linewidth',0.5)
        %     hold on
        % end
        % axis equal
    end
end
