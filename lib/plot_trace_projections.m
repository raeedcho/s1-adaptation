function plot_trace_projections(trial_data,params)
% plots 3d and 2d projections of traces and their averages
    % parameter setup
    signals = [];
    linestyle = '-';
    color = 'k';
    trials_to_plot = 1:length(trial_data);
    trials_to_use = 1:length(trial_data);
    alpha = 1;
    assignParams(who,params)
    
    % sanity check
    signals = check_signals(trial_data(1),signals);
    assert(size(signals,1)==1,'Can only use one signal field at a time currently')
    
    subplot_signals = repmat(signals,3,1);
    plot_permutes = {[1 2];[2 3];[1 3]};
    subplot_signals(:,2) = plot_permutes;
    for subplotnum = 1:3
        subplot(2,2,subplotnum)
        hold on
        plot_traces(trial_data,struct(...
            'signals',{subplot_signals(subplotnum,:)},...
            'plot_dim',2,...
            'linestyle',linestyle,...
            'color',color,...
            'alpha',alpha,...
            'trials_to_use',trials_to_use,...
            'trials_to_plot',trials_to_plot))
        axis equal
        set(gca,'box','off','tickdir','out')
        xlabel(sprintf('Component %d',plot_permutes{subplotnum}(1)))
        ylabel(sprintf('Component %d',plot_permutes{subplotnum}(2)))
    end
    
    subplot(2,2,4)
    hold on
    plot_traces(trial_data,struct(...
        'signals',{signals},...
        'plot_dim',3,...
        'linestyle',linestyle,...
        'color',color,...
        'alpha',alpha,...
        'trials_to_use',trials_to_use,...
        'trials_to_plot',trials_to_plot))
    view([-37.5 30])
    axis equal
    set(gca,'box','off','tickdir','out')
    xlabel('Component 1')
    ylabel('Component 2')
    zlabel('Component 3')
