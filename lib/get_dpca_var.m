function [margvar_table,learning_metric_table] = get_dpca_var(td,params)
    % helper function to get table of dPCA results (for bumpcurl_dpca)

    % set up params
    num_boots = 1;
    num_dims = 16;
    filetic = tic;
    shuffle_conds = {}; % use this to shuffle trials in one contition (e.g. shuffle_conds={'learning_block'} will shuffle learning block labels). Can be more than one
    plot_comp_dpcs = false;
    plot_marg = false;
    trim_start = -0.1;
    trim_end = 0.5;
    comp_to_plot = 'learning';
    learning_block_ranges = 0:0.1:1;
    colorfunc = @viridis;
    if nargin>1
        assignParams(who,params)
    end

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

    if isfield(td,'learning_block')
        % in case it was already assigned for some reason, we want to reassign it here
        td = rmfield(td,'learning_block');
    end
    learning_blocks = cell(1,length(learning_block_ranges)-1);
    for blocknum = 1:length(learning_blocks)
        learning_blocks{blocknum} = getTDidx(td,'range',learning_block_ranges([blocknum blocknum+1]));
        [td(learning_blocks{blocknum}).learning_block] = deal(blocknum);
    end

    % subselect trials for dpca
    all_block_inds = horzcat(learning_blocks{:});
    td_dpca = td(all_block_inds);

    % shuffle conditions if we want to
    if ~isempty(shuffle_conds)
        td_dpca = shuffle_td_labels(td_dpca,shuffle_conds);
    end

    % get actual bootstrapped dPCA
    margvar_array = cell(length(spikes_in_td),num_boots);
    learning_metric_array = cell(length(spikes_in_td),num_boots,length(learning_blocks));
    for arraynum = 1:length(spikes_in_td)
        for bootnum = 1:num_boots
            % inds to use in bootstrap
            if num_boots==1
                td_boot = td_dpca;
            else
                td_boot = randsample(td_dpca,length(td_dpca),true);
            end
            
            % get total neural covariance to start with
            data = getSig(td_boot,spikes_in_td(arraynum));
            neural_cov = cov(data);

            % check to make sure there are enough neurons
            if size(data,2)<=num_dims
                break;
            end
            
            % get dPCA results for this bootstrap iteration
            [td_boot, dpca_info] = runDPCA(td_boot,'target_block','learning_block',struct(...
                'signals',spikes_in_td(arraynum),...
                'combined_params',{{{1}, {2,[1 2]}, {3,[1 3],[2 3],[1 2 3]}}},...
                'marg_names',{{'time','target','learning'}},...
                'marg_colors',[150 150 150; 23 100 171; 187 20 25]/256,... % grey, blue, red
                'do_plot',plot_marg,'num_dims',num_dims,'out_sig_prefix',strcat(spikes_in_td{arraynum},'_','dpca')));
            % [td_boot, dpca_info] = runDPCA(td_boot,'target_block','learning_block',struct(...
            %     'signals',spikes_in_td(arraynum),...
            %     'marg_names',{{'time','target','learning','target_learning'}},...
            %     'do_plot',plot_marg,'num_dims',num_dims,'out_sig_prefix',strcat(spikes_in_td{arraynum},'_','dpca')));

            if plot_marg
                suptitle(sprintf('%s %s %s',td_boot(1).monkey,td_boot(1).date_time,spikes_in_td{arraynum}))
                % set(gcf,...
                %     'Name',sprintf('%s %s %s',td_boot(1).monkey,td_boot(1).date_time,spikes_in_td{arraynum}),...
                %     'NumberTitle','off')
            end

            % set up margvar table
            margvar_array{arraynum,bootnum} = table(...
                {td_boot(1).monkey},...
                {td_boot(1).date_time},...
                spikes_in_td(arraynum),...
                bootnum,...
                dpca_info.marg_names,...
                dpca_info.expl_var.totalMarginalizedVar,...
                100*dpca_info.expl_var.totalMarginalizedVar/sum(dpca_info.expl_var.totalMarginalizedVar),...
                dpca_info.which_marg,...
                dpca_info.combined_params,...
                {dpca_info.W},...
                {dpca_info.V},...
                {dpca_info.fr_tensor},...
                'VariableNames',{'monkey','date_time','array','bootID','marg_names','marg_var','marg_var_pct','which_marg','combined_params','decoder','encoder','fr_tensor'});
            margvar_array{arraynum,bootnum}.Properties.VariableDescriptions = {'meta','meta','meta','meta','meta','linear','linear','meta','meta','meta','meta','meta'};
            
            % set up learning metric table
            if nargout>1 && isfield(td_boot,sprintf('%s_dpca_learning',spikes_in_td{arraynum}))
                [~,td_block] = getTDidx(td_boot,'learning_block',length(learning_blocks));
                td_block = binTD(td_block,'average');
                
                learning_metric_end = mean(cat(1,td_block.learning_metric));
                data = get_vars(...
                        td_block,...
                        {sprintf('%s_dpca_%s',spikes_in_td{arraynum},comp_to_plot),1:2});
                learning_dPC_end = mean(data);
                
                for blocknum = unique([td_boot.learning_block])
                    [~,td_block] = getTDidx(td_boot,'learning_block',blocknum);
                    td_block = binTD(td_block,'average');
                    block_learning_metric = mean(cat(1,td_block.learning_metric));
                    
                    data = get_vars(...
                        td_block,...
                        {sprintf('%s_dpca_%s',spikes_in_td{arraynum},comp_to_plot),1:2});
                    block_learning_dPC = mean(data);
                    
                    rel_learning_dPC = block_learning_dPC-learning_dPC_end;
                    
                    learning_metric_array{arraynum,bootnum,blocknum} = table(...
                        {td_block(1).monkey},...
                        {td_block(1).date_time},...
                        spikes_in_td(arraynum),...
                        bootnum,...
                        blocknum,...
                        block_learning_metric-learning_metric_end,...
                        rel_learning_dPC,...
                        sqrt(sum(rel_learning_dPC.^2,2)),...
                        sqrt(sum(rel_learning_dPC.^2,2)/trace(neural_cov)),...
                        'VariableNames',{'monkey','date_time','array','bootID','learning_block','rel_learning_metric','rel_learning_dPC','learning_dPC_dist','learning_dPC_dist_norm'});
                    learning_metric_array{arraynum,bootnum,blocknum}.Properties.VariableDescriptions = {'meta','meta','meta','meta','meta','linear','linear','linear','linear'};
                end
            end
            
            % do a plot of the dPC space
            if plot_comp_dpcs && isfield(td_boot,sprintf('%s_dpca_%s',spikes_in_td{arraynum},comp_to_plot))
                block_labels = unique([td_boot.(sprintf('%s_block',comp_to_plot))]);
                num_blocks = length(block_labels);
                block_colors = colorfunc(num_blocks);
                figure
                
                % plot one time point per trial in dPC space
                td_avg = binTD(td_boot,'average');
                data = get_vars(td_avg,{sprintf('%s_dpca_%s',spikes_in_td{arraynum},comp_to_plot),1:2});
                scatter(data(:,1),data(:,2),[],colorfunc(length(data)),'filled')
                hold on
                
                for blocknum = 1:num_blocks
                    [~,td_block] = getTDidx(td_boot,sprintf('%s_block',comp_to_plot),block_labels(blocknum));
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

            if num_boots>1
                fprintf('Arraynum %d: Finished bootstrap iteration %d of %d at time %f\n',arraynum,bootnum,num_boots,toc(filetic))
            end
        end
    end
    margvar_table = vertcat(margvar_array{:});
    learning_metric_table = vertcat(learning_metric_array{:});
end

