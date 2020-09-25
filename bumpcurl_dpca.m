%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bumpcurl_analysis - 
%%      script to run bumpcurl task analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Plot trial info (hand speed and example rasters)
    for filenum = 1:length(trial_data_cell)
        % load and preprocess data
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
        title 'Washout'
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
    
%% Calculate variance due to learning using bootstrapping over trials
    margvar_cell = cell(length(trial_data_cell));
    num_boots = 1;
    num_dims = 15;
    filetic = tic;
    do_dpca_plot = (num_boots==1);
    trim_start = -0.1;
    trim_end = 0.5;
    comp_to_plot = 'target';
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

        % get learning block
        learning_blocks = {...
            getTDidx(td,'epoch','AD','range',[0 0.33]),...
            getTDidx(td,'epoch','AD','range',[0.33 0.67]),...
            getTDidx(td,'epoch','AD','range',[0.67 1]),...
            };
        
        for blocknum = 1:length(learning_blocks)
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
                
                % [~, dpca_info] = runDPCA(td_boot,'tgt_dir_block','learning_block',struct(...
                %     'signals',spikes_in_td(arraynum),...
                %     'marg_names',{{'time','direction','learning'}},...
                %     'combined_params',{{{1},{2,[1 2]},{3,[1 3],[2 3],[1 2 3]}}},...
                %     'marg_colors',marg_colors(1:3,:),...
                %     'do_plot',do_dpca_plot,'num_dims',10));
                % [~, dpca_info] = runDPCA(td_boot,'tgt_dir_block','learning_block',struct(...
                %     'signals',spikes_in_td(arraynum),...
                %     'marg_names',{{'learning','target/learning','other'}},...
                %     'combined_params',{{{3,[1 3]},{[2 3],[1 2 3]},{1,2,[1 2]}}},...
                %     'do_plot',do_dpca_plot,'num_dims',10));
                % [~, dpca_info] = runDPCA(td_boot,'tgt_dir_block','learning_block',struct(...
                %     'signals',spikes_in_td(arraynum),...
                %     'marg_names',{{'learning','other'}},...
                %     'combined_params',{{{3,[1 3],[2 3],[1 2 3]},{1,2,[1 2]}}},...
                %     'do_plot',do_dpca_plot,'num_dims',10));
                [td_boot, dpca_info] = runDPCA(td_boot,'target_block','learning_block',struct(...
                    'signals',spikes_in_td(arraynum),...
                    'marg_names',{{'time','target','learning','target_learning'}},...
                    'do_plot',false,'num_dims',num_dims,'out_sig_prefix',strcat(spikes_in_td{arraynum},'_','dpca')));
                
                if do_dpca_plot && isfield(td_boot,sprintf('%s_dpca_%s',spikes_in_td{arraynum},comp_to_plot))
                    num_blocks = length(unique([td_boot.(sprintf('%s_block',comp_to_plot))]));
                    block_colors = linspecer(num_blocks);
                    figure
                    for blocknum = unique([td_boot.(sprintf('%s_block',comp_to_plot))])
                        % plot(td_boot(trialnum).(sprintf('%s_dpca_learning',spikes_in_td{arraynum}))(:,1),'color',learning_colors(td_boot(trialnum).learning_block,:))
                        
%                         plot_traces(td_boot,struct(...
%                             'signals',{{'timevec',1;sprintf('%s_dpca_learning',spikes_in_td{arraynum}),1}},...
%                             'plot_dim',2,...
%                             'linestyle','-',...
%                             'color',learning_colors(blocknum,:),...
%                             'saturation',1,...
%                             'alpha',0.6,...
%                             'trials_to_use',getTDidx(td_boot,'learning_block',blocknum),...
%                             'trials_to_plot',getTDidx(td_boot,'learning_block',blocknum,'rand',20)))
                        plot_traces(td_boot,struct(...
                            'signals',{{sprintf('%s_dpca_%s',spikes_in_td{arraynum},comp_to_plot),1:2}},...
                            'plot_dim',2,...
                            'linestyle','-',...
                            'color',block_colors(blocknum,:),...
                            'saturation',1,...
                            'alpha',0.2,...
                            'trials_to_use',getTDidx(td_boot,sprintf('%s_block',comp_to_plot),blocknum),...
                            'trials_to_plot',getTDidx(td_boot,sprintf('%s_block',comp_to_plot),blocknum,'rand',20)))
                        
                        hold on
                    end
%                     xlabel('Time from movement onset (s)')
%                     ylabel('Normalized FR?')
                    xlabel(sprintf('%s dim 1',comp_to_plot))
                    ylabel(sprintf('%s dim 2',comp_to_plot))
                    title(sprintf('%s %s projection into %s dims',td_boot(1).monkey,strrep(spikes_in_td{arraynum},'_spikes',''),comp_to_plot))
                    set(gca,'box','off','tickdir','out')
                end

                margvar = dpca_info.expl_var.totalMarginalizedVar/dpca_info.expl_var.totalVar;
                margvar_array{arraynum,bootnum} = table(...
                    {td_boot(1).monkey},{td_boot(1).date_time},spikes_in_td(arraynum),bootnum,dpca_info.marg_names,margvar,...
                    'VariableNames',{'monkey','date_time','array','bootID','marg_names','margvar'});
                margvar_array{arraynum,bootnum}.Properties.VariableDescriptions = {'meta','meta','meta','meta','meta','linear'};
                fprintf('Filenum %d, arraynum %d: Finished bootstrap iteration %d of %d at time %f\n',filenum,arraynum,bootnum,num_boots,toc(filetic))
            end
        end
        margvar_cell{filenum} = vertcat(margvar_array{:});
    end
    margvar_table = vertcat(margvar_cell{:});
    
%% plot on learning axis (scratchpad)
    block_colors = linspecer(3);
    figure
    for trialnum = 1:length(td_boot)
        plot(td_boot(trialnum).PMd_spikes_dpca_learning(:,1),'color',block_colors(td_boot(trialnum).learning_block,:))
        hold on
    end

%% Try to decode what phase of adaptation we are in using LDA
    num_repeats = 20;
    num_folds = 5;
    file_results = cell(length(trial_data_cell),1);
    fprintf('Starting LDA classification analysis...\n')
    filetic = tic;
    for filenum = 1:length(trial_data_cell)
        %% load and preprocess
        td = trial_data_cell{filenum};
        
        % trim from go cue to end time (skip bump)
        spikes_in_td = getTDfields(td,'spikes');
        td = smoothSignals(td,struct('signals',{spikes_in_td}));
        if td(1).bin_size == 0.005
            td = binTD(td,2);
        end
        td = trimTD(td,struct(...
            'idx_start',{{'idx_movement_on',-10}},...
            'idx_end',{{'idx_movement_on',50}},...
            'remove_short',true));
        
        %% run LDA on individual time points
        % get learning conditions
        [~,td_ad] = getTDidx(td,'epoch','AD','range',[0 0.75]);
        learning_blocks = {...
            getTDidx(td_ad,'range',[0 0.33]),...
            getTDidx(td_ad,'range',[0.33 0.67]),...
            getTDidx(td_ad,'range',[0.67 1]),...
            };
        for blocknum = 1:length(learning_blocks)
            [td_ad(learning_blocks{blocknum}).learning_block] = deal(blocknum);
        end
        
        acc_table = cell(length(spikes_in_td),num_repeats,num_folds);
        arraytic = tic;
        for arraynum = 1:length(spikes_in_td)
            % make meta table columns
            meta_table = makeNeuronTableStarter(td_ad,struct('out_signal_names',{spikes_in_td(arraynum)}));
            
            % crossvalidate...
            repeattic = tic;
            for repeatnum = 1:num_repeats
                fold_inds = crossvalind('Kfold',length(td_ad),num_folds);
                foldtic = tic;
                for foldnum = 1:num_folds
                    td_test = td_ad(fold_inds==foldnum);
                    td_train = td_ad(fold_inds~=foldnum);
                    
                    fr_train = cat(3,td_train.(spikes_in_td{arraynum}));
                    fr_train = permute(fr_train,[3 2 1]);
                    class_train = cat(1,td_train.learning_block);
                    
                    fr_test = cat(3,td_test.(spikes_in_td{arraynum}));
                    fr_test = permute(fr_test,[3 2 1]);
                    class_test = cat(1,td_test.learning_block);
                    
                    % train and test lda
                    acc = zeros(1,size(fr_train,3));
                    for timepoint = 1:size(fr_train,3)
                        mdl = fitcdiscr(fr_train(:,:,timepoint),class_train,'discrimtype','pseudolinear');
                        preds = predict(mdl,fr_test(:,:,timepoint));
                        acc(timepoint) = sum(preds==class_test)/length(class_test);
                    end
                    
                    % put results into table
                    temp = table(repeatnum,foldnum,acc,'VariableNames',{'repeatID','foldID','class_accuracy'});
                    temp.Properties.VariableDescriptions = {'meta','meta','linear'};
                    acc_table{arraynum,repeatnum,foldnum} = horzcat(meta_table,temp);
                    % fprintf('      Finished crossval fold %d of %d at time %f\n',foldnum,num_folds,toc(foldtic))
                end
                fprintf('    Finished crossval repeat %d of %d at time %f\n',repeatnum,num_repeats,toc(repeattic))
            end
            fprintf('  Finished array %d of %d at time %f\n',arraynum,length(spikes_in_td),toc(arraytic))
        end
        file_results{filenum} = vertcat(acc_table{:});
        fprintf('Finished file %d at time %f\n',filenum,toc(filetic))
    end
    lda_class_acc = vertcat(file_results{:});
    
    %% average lda class accuracy at different time points of learning
    avg_class_table = neuronAverage(lda_class_acc,struct(...
        'keycols',{{'monkey','date','signalID'}},...
        'do_ci',true));
    timevec = (-10:50)*0.01;
    figure('defaultaxesfontsize',10)
    
    for arraynum = 1:height(avg_class_table)
        switch avg_class_table.signalID{arraynum}
            case 'M1_spikes'
                color = [102,194,165]/255;
            case 'PMd_spikes'
                color = [252,141,98]/255;
            case 'S1_spikes'
                color = [141,160,203]/255;
            otherwise
                color = 'k';
        end
%         patch(...
%             [timevec fliplr(timevec)],...
%             [avg_class_table.class_accuracyCILo(arraynum,:) fliplr(avg_class_table.class_accuracyCIHi(arraynum,:))],...
%             color,'facealpha',0.1,'edgecolor',color)
        hold on
        plot(timevec,avg_class_table.class_accuracy(arraynum,:),'linewidth',2,'color',color)
    end
    plot([timevec(1) timevec(end)],[0.33 0.33],'--k')
    xlabel('Time from movement onset')
    ylabel('Classification accuracy (among 3 classes)')
    
    set(gca,'box','off','tickdir','out')
    legend(strcat(avg_class_table.monkey,{'_'},avg_class_table.signalID))

