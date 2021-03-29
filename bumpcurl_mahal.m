function bumpcurl_mahal
% script to check distance in population space from final learning trials

%% Set up meta info
    if ispc
        homefolder = 'C:\Users\Raeed';
    else
        homefolder = '/home/raeed';
    end
    
    dataroot = fullfile(homefolder,'data','project-data','limblab','s1-adapt','td-library');
    file_info = dir(fullfile(dataroot,'*CO*.mat'));
    filenames = horzcat({file_info.name})';
    
    save_figures = true;
    figsavedir = fullfile(homefolder,'Wiki','professional','projects','s1-adaptation','figures','mahal');

%% Load up data into cell array
    trial_data_cell = load_curl_data(fullfile(dataroot,filenames));
    
%% Loop through trial data files
    num_dims = 16;
    trim_start = -0.1;
    trim_end = 0.5;
    num_control_shuffles = 1000;
    do_shuffle_control = true;
    
    for filenum = 4:length(filenames)
        td = trial_data_cell{filenum};
        [~,td] = getTDidx(td,'epoch','AD');
        
        % trim from go cue to end time (skip bump)
        spikes_in_td = getTDfields(td,'spikes');
        td = smoothSignals(td,struct('signals',{spikes_in_td},'width',0.1));
        td = trimTD(td,struct(...
            'idx_start',{{'idx_movement_on',trim_start/td(1).bin_size}},...
            'idx_end',{{'idx_movement_on',trim_end/td(1).bin_size}},...
            'remove_short',true));

        timevec = trim_start:td(1).bin_size:(length(td(1).vel)*td(1).bin_size+trim_start-0.001);
        for trialnum = 1:length(td)
            td(trialnum).timevec = timevec';
        end
        
        for arraynum = 1:length(spikes_in_td)
            % check to make sure there are enough neurons
            if size(getSig(td,spikes_in_td{arraynum}),2)<=num_dims
                continue;
            end
            
            true_mahal_curve = get_mahal_curve(td,struct(...
                'signals',spikes_in_td{arraynum},...
                'num_dims',num_dims));
            
            % get shuffle control
            if do_shuffle_control
                shuffle_mahal_curve = zeros(length(td),num_control_shuffles);
                for shufflenum = 1:num_control_shuffles
                    td_shuffle = shuffle_td_labels(td,'learning_block');
                    td_shuffle = td_shuffle(randperm(length(td_shuffle)));
                    try
                        shuffle_mahal_curve(:,shufflenum) = get_mahal_curve(td_shuffle,struct(...
                            'signals',spikes_in_td{arraynum},...
                            'num_dims',num_dims));
                    catch
                        continue
                    end
                end
            end
            
            fig = figure('defaultaxesfontsize',15,'position',[100 100 1000 500]);
            % plot out dPC space
            subplot(2,2,[1 3])
            plot_dpc_space(td,struct('signals',spikes_in_td{arraynum}))
            
            % plot out mahal distance over trials
            subplot(2,2,2)
            if do_shuffle_control
                shuffle_curve_bounds = prctile(shuffle_mahal_curve,[2.5 97.5],2);
                patch(...
                    [1:length(shuffle_curve_bounds) length(shuffle_curve_bounds):-1:1],...
                    [shuffle_curve_bounds(:,1); flipud(shuffle_curve_bounds(:,2))],...
                    [0.5 0.5 0.5],...
                    'facealpha',0.5)
                hold on
            end
            smoothing_window_size = 20;
            smoothing_kernel = ones(1,smoothing_window_size)/smoothing_window_size;
            smoothed_dist = conv(true_mahal_curve,smoothing_kernel,'same');
            plot(1:length(true_mahal_curve),smoothed_dist,'k','linewidth',2)
            ylabel({'Mahalanobis distance';'to final 40 trials'})
            set(gca,'box','off','tickdir','out','xlim',[0 length(td)+1])

            % plot out learning metric
            subplot(2,2,4)
            metric = cat(1,td.learning_metric);
            plot(conv(metric,smoothing_kernel,'same'),'k-')
            hold on
            plot([0 length(td)],[0 0],'k-','linewidth',2)
            ylabel('Learning metric')
            xlabel('Trial number')
            set(gca,'box','off','tickdir','out','xlim',[0 length(td)+1])
            
            suptitle(sprintf('%s %s %s',td(1).monkey,td(1).date_time,spikes_in_td{arraynum}))
            
            if save_figures
                sanitized_datetime = regexp(td(1).date_time,'^\d+[\/-]\d+[\/-]\d+','match');
                sanitized_datetime = strrep(sanitized_datetime{1},'/','');
                saveas(fig,fullfile(figsavedir,sprintf(...
                    '%s_%s_%s_dpc_mahal_dist.pdf',...
                    td(1).monkey,...
                    sanitized_datetime,...
                    strrep(spikes_in_td{arraynum},'_spikes',''))));
            end
        end
    end
end

%%%% subfunctions %%%%
function mahal_curve = get_mahal_curve(td,params)
% function to get a mahalanobis distace curve from the learning dimensions
% for a given trial_data structure, assuming there are target_block and
% learning_block labels
% inputs -
%   td - trial data structure
%   params - params struct
%       .signals - signals cell array or char array
%       .num_dims - how many dPCA dimensions to use
%       .plot_comp_dpcs - whether to plot dpc dots
%       .comp_to_plot - which components to plot
% Outputs - 
%   neural_tensor - tensor of neural activity in learning dimensions
%       (third order tensor, shape = time_points x learning_dims x trials)

    % params
    num_dims = 16;
    signals = '';
    assignParams(who,params);
    
    % run dPCA on data with target and learning blocks
    [td, ~] = runDPCA(td,'target_block','learning_block',struct(...
        'signals',signals,...
        'combined_params',{{{1}, {2,[1 2]}, {3,[1 3],[2 3],[1 2 3]}}},...
        'marg_names',{{'time','target','learning'}},...
        'marg_colors',[150 150 150; 23 100 171; 187 20 25]/256,... % grey, blue, red
        'do_plot',false,'num_dims',num_dims,'out_sig_prefix',strcat(signals,'_','dpca')));

    % get mahalanobis distance of each trial from final trials
    full_data = cat(3,td.(sprintf('%s_dpca_learning',signals)));
    % marginalize and subselect dims
    full_data = squeeze(mean(full_data,1))';
    full_data = full_data(:,1:4);

    % get ref data
    ref_data = full_data(end-40:end,:);

    mahal_curve = mahal(full_data,ref_data);
end

function plot_dpc_space(td,params)
% function to plot out learning space of average neural activity per trial, separated by learning block

    num_dims = 16;
    signals = '';
    comp_to_plot = 'learning';
    colorfunc = @viridis;
    assignParams(who,params);
    
    % run dPCA on data with target and learning blocks
    [td, ~] = runDPCA(td,'target_block','learning_block',struct(...
        'signals',signals,...
        'combined_params',{{{1}, {2,[1 2]}, {3,[1 3],[2 3],[1 2 3]}}},...
        'marg_names',{{'time','target','learning'}},...
        'marg_colors',[150 150 150; 23 100 171; 187 20 25]/256,... % grey, blue, red
        'do_plot',false,'num_dims',num_dims,'out_sig_prefix',strcat(signals,'_','dpca')));
    
    assert(isfield(td,sprintf('%s_dpca_%s',signals,comp_to_plot)),'Something went wrong')
        
    block_labels = unique([td.(sprintf('%s_block',comp_to_plot))]);
    num_blocks = length(block_labels);
    block_colors = colorfunc(num_blocks);

    % plot one time point per trial in dPC space
    td_avg = binTD(td,'average');
    data = get_vars(td_avg,{sprintf('%s_dpca_%s',signals,comp_to_plot),1:2});
    scatter(data(:,1),data(:,2),[],colorfunc(length(data)),'filled')
    hold on

    for blocknum = 1:num_blocks
        [~,td_block] = getTDidx(td,sprintf('%s_block',comp_to_plot),block_labels(blocknum));
        td_block = binTD(td_block,'average');

        data = get_vars(...
            td_block,...
            {sprintf('%s_dpca_%s',signals,comp_to_plot),1:2});

        plotErrorEllipse(mean(data,1),cov(data),0.95,'color',block_colors(blocknum,:),'linewidth',2)
    end
    xlabel(sprintf('%s dim 1',comp_to_plot))
    ylabel(sprintf('%s dim 2',comp_to_plot))
    title(sprintf('Projection into %s dims',comp_to_plot))
end

function plot_mahal_curve(mahal_table)
% function to plot out the mahalnobis distance curve

    smoothing_window_size = 20;
    smoothing_kernel = ones(1,smoothing_window_size)/smoothing_window_size;
    
    smoothed_dist = conv(true_mahal_curve,smoothing_kernel,'same');
    figure
    subplot(2,1,1)
    if do_shuffle_control
        shuffle_curve_bounds = prctile(shuffle_mahal_curve,[2.5 97.5],2);
        patch(...
            [1:length(shuffle_curve_bounds) length(shuffle_curve_bounds):-1:1],...
            [shuffle_curve_bounds(:,1); flipud(shuffle_curve_bounds(:,2))],...
            [0.5 0.5 0.5],...
            'facealpha',0.5)
        hold on
    end
    plot(1:length(true_mahal_curve),smoothed_dist,'k','linewidth',2)
    set(gca,'box','off','tickdir','out')
    ylabel('Mahalanobis distance to final 40 trials in dPC space')
    title(sprintf('%s %s %s',td(1).monkey,td(1).date_time,spikes_in_td{arraynum}))

    subplot(2,1,2)
    metric = cat(1,td.learning_metric);
    plot(conv(metric,smoothing_kernel,'same'),'k-')
    hold on
    plot([0 length(td)],[0 0],'k-','linewidth',2)
    xlabel('Trial number')
end