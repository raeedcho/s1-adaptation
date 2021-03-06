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
    % plot_learning_behavior(trial_data_cell);
    
%% Calculate variance due to learning using bootstrapping over trials
    [margvar_cell,margvar_bl_cell,learning_metric_cell] = deal(cell(length(trial_data_cell)));
    for filenum = [1:4 6:11]
        td = trial_data_cell{filenum};
        [~,td_ad] = getTDidx(td,'epoch','AD');
        ad_trial_frac = 0.2;
        [margvar_cell{filenum},learning_metric_cell{filenum}] = get_dpca_var(td_ad,struct(...
            'num_boots',1,...
            'plot_marg',false,...
            'plot_learning_dpcs',false,...
            'learning_block_ranges',0:ad_trial_frac:1));

        % get baseline values for metrics
        [~,td_bl] = getTDidx(td,'epoch','BL');
        % bl_trial_frac = ad_trial_frac*length(td_ad)/length(td_bl);
        bl_trial_frac = 0.2;
        [margvar_bl_cell{filenum},~] = get_dpca_var(td_bl,struct(...
            'num_boots',1,...
            'plot_marg',false,...
            'plot_learning_dpcs',false,...
            'learning_block_ranges',0:bl_trial_frac:1));
    end
    margvar_table = vertcat(margvar_cell{:});
    margvar_bl_table = vertcat(margvar_bl_cell{:});
    learning_metric_table = vertcat(learning_metric_cell{:});
    
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
    figure('defaultaxesfontsize',18)
    arrays = {'S1_spikes','PMd_spikes','M1_spikes'};
    monkeys = {'Chewie','Mihili','MrT','Duncan','Han'};
    monkey_colors = linspecer(length(monkeys));
    for arraynum = 1:length(arrays)
        [~,array_metric_table] = getNTidx(margvar_table,'array',arrays{arraynum});
        sessions = unique(array_metric_table.date_time);
        for sessionnum = 1:length(sessions)
            [~,session_table] = getNTidx(array_metric_table,'date_time',sessions{sessionnum});
            bl_tab_idx = strncmp(margvar_bl_table.date_time,sessions{sessionnum},10) & ...
                strcmpi(margvar_bl_table.array,arrays{arraynum});
            bl_table = margvar_bl_table(bl_tab_idx,:);
            assert(height(session_table)==1 && height(bl_table)==1) % assume all sessions are on different days
            point_color = monkey_colors(strcmpi(session_table.monkey,monkeys),:);

            % % figure out which td file we're on
            % for filenum = 1:length(filenames)
            %     if strcmpi(trial_data_cell{filenum}(1).monkey,session_table.monkey)
            %         if strncmp(trial_data_cell{filenum}(1).date_time,session_table.date_time,10)
            %             td_select = trial_data_cell{filenum};
            %             td_select = smoothSignals(td_select,struct('signals',arrays{arraynum},'width',0.1));
            %             if td_select(1).bin_size == 0.005
            %                 td_select = binTD(td_select,2);
            %             end
            %             td_select = trimTD(td_select,struct(...
            %                 'idx_start',{{'idx_movement_on',-0.1/td_select(1).bin_size}},...
            %                 'idx_end',{{'idx_movement_on',0.5/td_select(1).bin_size}},...
            %                 'remove_short',true));
            %             [~,td_bl] = getTDidx(td_select,'epoch','BL');
            %             [~,td_ad] = getTDidx(td_select,'epoch','AD');
            %             break
            %         end
            %     end
            % end

            % % calculate marg var for adaptation
            % ad_expl_var = dpca_explainedVariance(...
            %     session_table.fr_tensor{1},...
            %     session_table.decoder{1}, session_table.encoder{1},...
            %     'combinedParams', session_table.combined_params);
            % bl_expl_var = dpca_explainedVariance(...
            %     bl_table.fr_tensor{1},...
            %     session_table.decoder{1}, session_table.encoder{1},...
            %     'combinedParams', session_table.combined_params);
            % bl_expl_var = dpca_explainedVariance(...
            %     bl_table.fr_tensor{1},...
            %     bl_table.decoder{1}, bl_table.encoder{1},...
            %     'combinedParams', bl_table.combined_params);

            for margnum = 1:length(session_table.marg_names)
                % ad_val = sum(ad_expl_var.margVar(margnum,session_table.which_marg==margnum),2);
                % bl_val = sum(bl_expl_var.margVar(margnum,session_table.which_marg==margnum),2);

                ad_val = 100*session_table.marg_var(margnum)/sum(session_table.marg_var);
                bl_val = 100*bl_table.marg_var(margnum)/sum(bl_table.marg_var);

                % spike_data_ad = getSig(td_ad,arrays{arraynum});
                % spike_data_ad = spike_data_ad-mean(spike_data_ad);
                % spike_data_bl = getSig(td_bl,arrays{arraynum});
                % spike_data_bl = spike_data_bl-mean(spike_data_bl);
                % marg_proj_ad = spike_data_ad * session_table.decoder{1}(:,session_table.which_marg==margnum);
                % marg_proj_bl = spike_data_bl * session_table.decoder{1}(:,session_table.which_marg==margnum);
                % spike_var_ad = sum(sum(spike_data_ad.^2));
                % spike_var_bl = sum(sum(spike_data_bl.^2));
                % marg_pred_ad = marg_proj_ad * transpose(session_table.encoder{1}(:,session_table.which_marg==margnum));
                % marg_pred_bl = marg_proj_bl * transpose(session_table.encoder{1}(:,session_table.which_marg==margnum));
                % ad_val = 100*(1-sum(sum((spike_data_ad-marg_pred_ad).^2))/spike_var_ad);
                % bl_val = 100*(1-sum(sum((spike_data_bl-marg_pred_bl).^2))/spike_var_bl);
                % ad_val = 100*sum(sum(marg_pred_ad.^2))/spike_var_ad;
                % bl_val = 100*sum(sum(marg_pred_bl.^2))/spike_var_bl;

                subplot(1,length(session_table.marg_names),margnum)
                hold on
                scatter(...
                    arraynum+0.1,...
                    ad_val,...
                    [],point_color,'filled')
                scatter(...
                    arraynum-0.1,...
                    bl_val,...
                    [],point_color)
                plot(...
                    arraynum+[0.1 -0.1],...
                    horzcat(ad_val, bl_val),...
                    '-','color',point_color)
            end
        end
    end
    for margnum = 1:length(margvar_table(1,:).marg_names)
        subplot(1,length(margvar_table(1,:).marg_names),margnum)
        set(gca,...
            'box','off',...
            'tickdir','out',...
            'ylim',[0 100],'xlim',[0 length(arrays)+1],...
            'xtick',1:length(arrays),'xticklabel',arrays,...
            'ticklabelinterpreter','none')
        xlabel('Brain region')
        ylabel('% of total neural variance')
        title(margvar_table(1,:).marg_names{margnum})
    end
    suptitle('Comparison between variance in (projected) baseline and adaptation')

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
        smoothing_window_size = 5;
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

