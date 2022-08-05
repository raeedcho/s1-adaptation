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
    
    save_figures = false;
    figsavedir = fullfile(homefolder,'Wiki','professional','projects','s1-adaptation','figures','mahal');

%% Load up data into cell array
    trial_data_cell = load_curl_data(fullfile(dataroot,filenames));
    
%% Loop through trial data files
    num_dims = 16;
    num_mahal_dims = 4;
%     alignment_event = 'idx_movement_on';
%     trim_start = -0.1;
%     trim_end = 0.5;
    trim_start = 0;
    trim_end = 0.3;
    alignment_event = 'idx_bumpTime';
    num_control_shuffles = 1000;
    do_shuffle_control = false;
    arrays_to_plot = {'M1','S1'};
    smoothing_window_size = 10;
    smoothing_kernel = ones(1,smoothing_window_size)/smoothing_window_size;
    
    mahal_curve_file=cell(length(filenames),1);
    
    for filenum = [3 6]%1:length(filenames)
        td = trial_data_cell{filenum};
        [~,td] = getTDidx(td,'epoch','AD');
        
        % trim from go cue to end time (skip bump)
        spikes_in_td = getTDfields(td,'spikes');
        td = smoothSignals(td,struct('signals',{spikes_in_td},'width',0.05));

        valid_trials = ~isnan(cat(1,td.(alignment_event)));
        td = trimTD(td(valid_trials),struct(...
            'idx_start',{{alignment_event,trim_start/td(1).bin_size}},...
            'idx_end',{{alignment_event,trim_end/td(1).bin_size}},...
            'remove_short',true));

        timevec = trim_start:td(1).bin_size:(length(td(1).vel)*td(1).bin_size+trim_start-0.001);
        for trialnum = 1:length(td)
            td(trialnum).timevec = timevec';
        end
        
        mahal_curve_arr = cell(1,length(spikes_in_td));
        for arraynum = 1:length(spikes_in_td)
            % check to make sure there are enough neurons
            if size(getSig(td,spikes_in_td{arraynum}),2)<=num_dims
                continue;
            end
            
            % check to only plot the arrays we want
            if ~contains(spikes_in_td{arraynum},arrays_to_plot)
                continue
            end
            
            true_mahal_curve = get_mahal_curve(td,struct(...
                'signals',spikes_in_td{arraynum},...
                'num_dims',num_dims,...
                'num_mahal_dims',num_mahal_dims,...
                'use_dpca',true ...
            ));
            
            initial_mahal_dist = mean(true_mahal_curve(1:floor(length(true_mahal_curve)*0.1)));
            
            % get shuffle control
            if do_shuffle_control
                shuffle_mahal_curve = zeros(length(td),num_control_shuffles);
                for shufflenum = 1:num_control_shuffles
                    td_shuffle = td;
%                     td_shuffle = shuffle_td_labels(td_shuffle,'learning_block');
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
            
            fig = figure(...
                'defaultaxesfontsize',10,...
                'units','inches',...
                'position',[3 3 6 3],...
                'paperunits','inches',...
                'papersize',[6 3]);
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
            else % use a simple chi-squared distribution to estimate 95% CI
                mahal_thresh = chi2inv(0.95,num_mahal_dims);
                plot([1 length(true_mahal_curve)],repmat(mahal_thresh,1,2),'--','linewidth',2,'color',[0.5 0.5 0.5])
                text(length(true_mahal_curve)/2,mahal_thresh,'95% CI bound','FontSize',10,'color',[0.5 0.5 0.5])
            end
            hold on
            
            smoothed_dist = conv(true_mahal_curve,smoothing_kernel,'same');
            plot(1:length(true_mahal_curve),smoothed_dist,'k','linewidth',2)
            text_idx = floor(length(true_mahal_curve)/10);
            text(text_idx,smoothed_dist(text_idx),'Mahalanobis distance to final 40 trials','fontsize',10)
%             ylabel({'Mahalanobis distance';'to final 40 trials'})
            set(gca,'box','off','tickdir','out','xlim',[0 length(td)+1],'xticklabel',{})

            % plot out learning metric
            subplot(2,2,4)
            metric = cat(1,td.learning_metric)*180/pi;
            plot(conv(metric,smoothing_kernel,'same'),'k-','linewidth',2)
            hold on
            plot([0 length(td)],[0 0],'k-','linewidth',2)
            ylabel('Takeoff angle error (deg)')
            xlabel('Trial number')
            set(gca,'box','off','tickdir','out','xlim',[0 length(td)+1])
            
            sgtitle(sprintf('%s %s %s',td(1).monkey,td(1).date_time,strrep(spikes_in_td{arraynum},'_spikes','')))
            
            if save_figures
                sanitized_datetime = regexp(td(1).date_time,'^\d+[\/-]\d+[\/-]\d+','match');
                sanitized_datetime = strrep(sanitized_datetime{1},'/','');
                saveas(fig,fullfile(figsavedir,sprintf(...
                    '%s_%s_%s_dpc_mahal_dist.pdf',...
                    td(1).monkey,...
                    sanitized_datetime,...
                    strrep(spikes_in_td{arraynum},'_spikes',''))));
            end
            
            mahal_curve_arr{arraynum} = table(...
                {td(1).monkey},...
                {td(1).date_time},...
                {strrep(spikes_in_td{arraynum},'_spikes','')},...
                {true_mahal_curve'},...
                initial_mahal_dist,...
                'VariableNames',{'monkey','date_time','array','mahal_curve','initial_mahal_dist'});
        end
        mahal_curve_file{filenum} = vertcat(mahal_curve_arr{:});
    end
    mahal_curve_table = vertcat(mahal_curve_file{:});
    
    
%% plot out all mahal curves, colored by array
    compiled_mahal_fig = figure(...
        'defaultaxesfontsize',10,...
        'units','inches',...
        'position',[3 3 6 3],...
        'paperunits','inches',...
        'papersize',[6 3]);
    ax(1) = subplot(1,3,1:2);
    
    patch([0 0.1 0.1 0],[0 0 60 60],[0.8 0.8 0.8],'edgecolor','none')
    hold on
    
    mahal_thresh = chi2inv(0.95,num_mahal_dims);
    plot([0 1],repmat(mahal_thresh,1,2),'--','linewidth',2,'color',[0.5 0.5 0.5])
    text(0.7,1.4*mahal_thresh,'95% CI bound','FontSize',10,'color',[0.5 0.5 0.5])
    
    array_colors = linspecer(length(arrays_to_plot));
    
    for arraynum = 1:length(arrays_to_plot)
        [~,mahal_curve_array] = getNTidx(mahal_curve_table,'array',arrays_to_plot{arraynum});
        for sessionnum = 1:height(mahal_curve_array)
            temp_mahal_curve = mahal_curve_array{sessionnum,'mahal_curve'}{1};
            plot(...
                (1:length(temp_mahal_curve))/length(temp_mahal_curve),...
                conv(temp_mahal_curve,smoothing_kernel,'same'),...
                'linewidth',2,'color',array_colors(arraynum,:))
        end
        text(0.25*arraynum,(4-arraynum)*mahal_thresh,arrays_to_plot{arraynum},'FontSize',10,'color',array_colors(arraynum,:))
    end
    xlabel('Fraction into adaptation epoch')
    ylabel('Mahal dist to final 40 trials')
    set(gca,'box','off','tickdir','out','xtick',0.2:0.2:1,'ytick',20:20:60)
    
    ax(2) = subplot(1,3,3);
    [~,array_idx] = ismember(mahal_curve_table{:,'array'},arrays_to_plot);
    plot([0.5 0.5+length(arrays_to_plot)],repmat(mahal_thresh,1,2),'--','linewidth',2,'color',[0.5 0.5 0.5])
    hold on
    text(0.7,1.4*mahal_thresh,'95% CI bound','FontSize',10,'color',[0.5 0.5 0.5])
    scatter(array_idx,mahal_curve_table{:,'initial_mahal_dist'},[],array_colors(array_idx,:),'filled')
    ylabel('Distance in first 10% of adaptation')
    set(gca,'box','off','tickdir','out',...
        'xlim',[0.5 0.5+length(arrays_to_plot)],'xtick',1:length(arrays_to_plot),'xticklabel',arrays_to_plot,...
        'ytick',20:20:60,'yticklabel',{})
    linkaxes(ax,'y')
    
    if save_figures
        saveas(compiled_mahal_fig,fullfile(figsavedir,'dpc_mahal_dist_compiled.pdf'));
    end
    
end

%%%% subfunctions %%%%
function mahal_curve = get_mahal_curve(td,params)
% function to get a squared mahalanobis distace curve from the learning dimensions
% for a given trial_data structure, assuming there are target_block and
% learning_block labels
% inputs -
%   td - trial data structure
%   params - params struct
%       .signals - signals cell array or char array
%       .num_dims - how many dPCA dimensions to use
%       .num_mahal_dims - how many dimensions to calculate mahalanobis distance
%       .plot_comp_dpcs - whether to plot dpc dots
%       .comp_to_plot - which components to plot
% Outputs - 
%   mahal_curve - array of squared mahalanobis distances (one for each trial)

    % params
    num_dims = 16;
    num_mahal_dims = 4;
    signals = '';
    use_dpca = true;
    assignParams(who,params);
    
    if use_dpca
        % run dPCA on data with target and learning blocks
        [td, ~] = runDPCA(td,'target_block','learning_block',struct(...
            'signals',signals,...
            'combined_params',{{{1}, {2,[1 2]}, {3,[1 3],[2 3],[1 2 3]}}},...
            'marg_names',{{'time','target','learning'}},...
            'marg_colors',[150 150 150; 23 100 171; 187 20 25]/256,... % grey, blue, red
            'do_plot',false,'num_dims',num_dims,'out_sig_prefix',strcat(signals,'_','dpca')));
    
        full_data = cat(3,td.(sprintf('%s_dpca_learning',signals)));
    else
        [td,~] = dimReduce(td,struct(...
            'algorithm','pca',...
            'signals',signals,...
            'num_dims',num_dims ...
        ));

        full_data = cat(3,td.(strrep(signals,'_spikes','_pca')));
    end

    % get mahalanobis distance of each trial from final trials
    % marginalize and subselect dims
    full_data = squeeze(mean(full_data,1))';
    full_data = full_data(:,1:num_mahal_dims);

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
    
    xlims = get(gca,'xlim');
    ylims = get(gca,'ylim');
    plot([xlims(1) xlims(1)+0.2],[ylims(1) ylims(1)],'k','linewidth',2)
    plot([xlims(1) xlims(1)],[ylims(1) ylims(1)+0.2],'k','linewidth',2)
    text(xlims(1)+0.05,ylims(1)+0.05,'0.2','FontSize',10)
    
    set(gca,'box','off','tickdir','out')
    axis equal
    axis off
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