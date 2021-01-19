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
    [~,learning_metric_table] = get_dpca_var(trial_data_cell([1 2 3 4 6 7]));
    
%% compare learning metric with dPC movement
    arrays = {'PMd_spikes','M1_spikes','S1_spikes'};
    for arraynum = 1:length(arrays)
        [~,array_metric_table] = getNTidx(learning_metric_table,'array',arrays{arraynum});
        sessions = unique(array_metric_table.date_time);
        session_colors = linspecer(length(sessions));
        figure('defaultaxesfontsize',18)
        for sessionnum = 1:length(sessions)
            [~,session_table] = getNTidx(array_metric_table,'date_time',sessions{sessionnum});
            hold on
            scatter(session_table.rel_learning_metric,session_table.learning_dPC_dist_norm,[],session_colors(sessionnum,:),'filled')
        end
        set(gca,'box','off','tickdir','out','ylim',[0 0.75],'xlim',[-0.35 0.35])
        title(arrays{arraynum})
        xlabel('Learning metric diff')
        ylabel('Norm learning dPC diff')
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

