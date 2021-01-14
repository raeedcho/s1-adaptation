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
    dpca_results = get_dpca_var(trial_data_cell([1 2 3 4 6 7]));
    
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

function margvar_table = get_dpca_var(trial_data_cell,params)
    % helper function to get table of dPCA results

    % set up params
    margvar_cell = cell(length(trial_data_cell));
    num_boots = 1;
    num_dims = 16;
    filetic = tic;
    do_dpca_plot = (num_boots==1);
    trim_start = -0.1;
    trim_end = 0.5;
    comp_to_plot = 'learning';
    learning_block_ranges = 0:0.1:1;
    colorfunc = @viridis;
    if nargin>1
        assignParams(who,params)
    end

    for filenum = 1:length(trial_data_cell)
        % load and preprocess data
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

        learning_blocks = cell(1,length(learning_block_ranges)-1);
        for blocknum = 1:length(learning_blocks)
            learning_blocks{blocknum} = getTDidx(td,'epoch','AD','range',learning_block_ranges([blocknum blocknum+1]));
            [td(learning_blocks{blocknum}).learning_block] = deal(blocknum);
        end
        
        % subselect trials for dpca
        all_block_inds = horzcat(learning_blocks{:});
        td_dpca = td(all_block_inds);

        % get actual bootstrapped dPCA
        margvar_array = cell(length(spikes_in_td),num_boots);
        for arraynum = 1:length(spikes_in_td)
            for bootnum = 1:num_boots
                % inds to use in bootstrap
                if num_boots==1
                    td_boot = td_dpca;
                else
                    td_boot = randsample(td_dpca,length(td_dpca),true);
                end
                
                % get dPCA results for this bootstrap iteration
                [td_boot, dpca_info] = runDPCA(td_boot,'target_block','learning_block',struct(...
                    'signals',spikes_in_td(arraynum),...
                    'marg_names',{{'time','target','learning','target_learning'}},...
                    'do_plot',false,'num_dims',num_dims,'out_sig_prefix',strcat(spikes_in_td{arraynum},'_','dpca')));

                margvar = dpca_info.expl_var.totalMarginalizedVar/dpca_info.expl_var.totalVar;
                margvar_array{arraynum,bootnum} = table(...
                    {td_boot(1).monkey},{td_boot(1).date_time},spikes_in_td(arraynum),bootnum,dpca_info.marg_names,margvar,...
                    'VariableNames',{'monkey','date_time','array','bootID','marg_names','margvar'});
                margvar_array{arraynum,bootnum}.Properties.VariableDescriptions = {'meta','meta','meta','meta','meta','linear'};
                
                if do_dpca_plot && isfield(td_boot,sprintf('%s_dpca_%s',spikes_in_td{arraynum},comp_to_plot))
                    num_blocks = length(unique([td_boot.(sprintf('%s_block',comp_to_plot))]));
                    block_colors = colorfunc(num_blocks);
                    figure
                    
                    % plot one time point per trial in dPC space
                    td_avg = binTD(td_boot,'average');
                    data = get_vars(td_avg,{sprintf('%s_dpca_%s',spikes_in_td{arraynum},comp_to_plot),1:2});
                    scatter(data(:,1),data(:,2),[],colorfunc(length(data)),'filled')
                    hold on
                    
                    for blocknum = unique([td_boot.(sprintf('%s_block',comp_to_plot))])
                        [~,td_block] = getTDidx(td_boot,sprintf('%s_block',comp_to_plot),blocknum);
                        td_block = binTD(td_block,'average');
                        
                        data = get_vars(...
                            td_block,...
                            {sprintf('%s_dpca_%s',spikes_in_td{arraynum},comp_to_plot),1:2});
                        
                        plotErrorEllipse(mean(data,1),cov(data),0.95,'color',block_colors(blocknum,:),'linewidth',2)
                    end
                    xlabel(sprintf('%s dim 1',comp_to_plot))
                    ylabel(sprintf('%s dim 2',comp_to_plot))
                    title(sprintf('%s %s projection into %s dims',td_boot(1).monkey,strrep(spikes_in_td{arraynum},'_spikes',''),comp_to_plot))
                    set(gca,'box','off','tickdir','out')
                end
                fprintf('Filenum %d, arraynum %d: Finished bootstrap iteration %d of %d at time %f\n',filenum,arraynum,bootnum,num_boots,toc(filetic))
            end
        end
        margvar_cell{filenum} = vertcat(margvar_array{:});
    end
    margvar_table = vertcat(margvar_cell{:});
end
