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