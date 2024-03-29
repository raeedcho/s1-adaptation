function mahal_curve_table = iterate_mahal_curves(trial_data_cell,params)

    alignment_event = 'idx_movement_on';
    arrays_to_plot = {'M1','S1'};
    num_dims = 16;
    num_mahal_dims = 4;
    save_figures = false;
    figsavedir = '';
    assignParams(who,params)

    if strcmpi(alignment_event,'idx_movement_on')
        trim_start = -0.1;
        trim_end = 0.5;
    elseif strcmpi(alignment_event,'idx_bumpTime')
        trim_start = 0;
        trim_end = 0.3;
    end
    
    mahal_curve_file=cell(length(trial_data_cell),1);
    for filenum = 1:length(trial_data_cell)
        td = trial_data_cell{filenum};
        [~,td] = getTDidx(td,'epoch','AD');

        if ~isfield(td,alignment_event)
            continue
        end
        
        spikes_in_td = getTDfields(td,'spikes');
        td = smoothSignals(td,struct('signals',{spikes_in_td},'width',0.05));

        valid_trials = ~isnan(cat(1,td.(alignment_event)));
        if ~any(valid_trials)
            continue
        end
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
            % check to only plot the arrays we want
            if ~contains(spikes_in_td{arraynum},arrays_to_plot)
                continue
            end

            mahal_curve_arr{arraynum} = get_mahal_curve(td,struct(...
                'signals',spikes_in_td{arraynum},...
                'num_dims',num_dims,...
                'num_mahal_dims',num_mahal_dims,...
                'do_plot',true,...
                'save_figures',save_figures,...
                'figsavedir',figsavedir...
            ));
        end
        mahal_curve_file{filenum} = vertcat(mahal_curve_arr{:});
    end
    mahal_curve_table = vertcat(mahal_curve_file{:});
end


%%%% subfunctions %%%%

function mahal_curve_row = get_mahal_curve(td,params)
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
    do_plot = true;
    save_figures = false;
    figsavedir = '';
    assignParams(who,params);

    % check to make sure there are enough neurons
    if size(getSig(td,signals),2)<=num_dims
        mahal_curve_row = table();
        return
    end
    
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

    % calculate mahal distance and starting point
    mahal_curve = mahal(full_data,ref_data);
    initial_mahal_dist = mean(mahal_curve(1:floor(length(mahal_curve)*0.1)));

    % package output
    mahal_curve_row = table(...
        {td(1).monkey},...
        {td(1).date_time},...
        {strrep(signals,'_spikes','')},...
        {mahal_curve'},...
        initial_mahal_dist,...
        'VariableNames',{'monkey','date_time','array','mahal_curve','initial_mahal_dist'});

    % if we want to plot
    if do_plot
        fig = plot_session_mahal(td,mahal_curve,params);
        if save_figures
            sanitized_datetime = regexp(td(1).date_time,'^\d+[\/-]\d+[\/-]\d+','match');
            sanitized_datetime = strrep(sanitized_datetime{1},'/','');
            saveas(fig,fullfile(figsavedir,sprintf(...
                '%s_%s_%s_dpc_mahal_dist.pdf',...
                td(1).monkey,...
                sanitized_datetime,...
                strrep(signals,'_spikes',''))));
        end
    end
end

function fig = plot_session_mahal(td,mahal_curve,params)

    signals = '';
    num_mahal_dims = 4;
    smoothing_window_size = 10;
    assignParams(who,params);

    fig = figure(...
        'defaultaxesfontsize',10,...
        'units','inches',...
        'position',[3 3 6 3],...
        'paperunits','inches',...
        'papersize',[6 3]);
    % plot out dPC space
    subplot(2,2,[1 3])
    plot_dpc_space(td,params)
    
    % plot out mahal distance over trials
    subplot(2,2,2)
    plot_mahal_curve(mahal_curve,params)
    xlabel('')

    % plot out learning metric
    subplot(2,2,4)
    plot_learning_metric(td,params)
    
    sgtitle(sprintf('%s %s %s',td(1).monkey,td(1).date_time,strrep(signals,'_spikes','')))
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

function plot_mahal_curve(mahal_curve,params)
% function to plot out the mahalnobis distance curve

    num_mahal_dims = 4;
    smoothing_window_size = 10;
    assignParams(who,params)

    smoothing_kernel = ones(1,smoothing_window_size)/smoothing_window_size;

    % use a simple chi-squared distribution to estimate 95% CI
    mahal_thresh = chi2inv(0.95,num_mahal_dims);
    plot([1 length(mahal_curve)],repmat(mahal_thresh,1,2),'--','linewidth',2,'color',[0.5 0.5 0.5])
    text(length(mahal_curve)/2,mahal_thresh,'95% CI bound','FontSize',10,'color',[0.5 0.5 0.5])

    hold on
    
    smoothed_dist = conv(mahal_curve,smoothing_kernel,'same');
    plot(1:length(mahal_curve),smoothed_dist,'k','linewidth',2)
    text_idx = floor(length(mahal_curve)/10);
    text(text_idx,smoothed_dist(text_idx),'Mahalanobis distance to final 40 trials','fontsize',10)
    xlabel('Trial number')
    set(gca,'box','off','tickdir','out','xlim',[0 length(mahal_curve)+1],'xticklabel',{})
end

function plot_learning_metric(td,params)
    
    smoothing_window_size = 10;
    assignParams(who,params)

    smoothing_kernel = ones(1,smoothing_window_size)/smoothing_window_size;

    metric = cat(1,td.learning_metric)*180/pi;
    plot(conv(metric,smoothing_kernel,'same'),'k-','linewidth',2)
    hold on
    plot([0 length(td)],[0 0],'k-','linewidth',2)
    ylabel('Takeoff angle error (deg)')
    xlabel('Trial number')
    set(gca,'box','off','tickdir','out','xlim',[0 length(td)+1])
end


