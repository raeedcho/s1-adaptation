function plotTensorDecomp(M, params)
    
    % default params
    timevec = [];
    trial_colors = 'k';

    if nargin>1
        assignParams(who,params);
    end

    figure('position',[50 50 800 1200])
    ncol = 3;
    nrow = size(M.U{1},2);

    % Look at signal factors
    signal_factors = M.U{2};
    markervec = 1:size(signal_factors,1);
    for i = 1:size(signal_factors,2)
        subplot(nrow,ncol,(i-1)*ncol+1)
        bar(markervec,signal_factors(:,i))
        set(gca,'box','off','tickdir','out')
    end

    % Plot temporal factors
    temporal_factors = M.U{1};
    if isempty(timevec)
        timevec = 1:size(temporal_factors,1);
    end
    for i = 1:size(temporal_factors,2)
        subplot(nrow,ncol,(i-1)*ncol+2)
        plot(timevec,temporal_factors(:,i),'-k','linewidth',3)
        hold on
        plot(timevec([1 end]),[0 0],'-k','linewidth',2)
        plot([0 0],ylim,'--k','linewidth',2)
        set(gca,'box','off','tickdir','out')
    end

    % Look at trial factors
    trial_factors = M.U{3};
    trialvec = 1:size(trial_factors,1);
    for i = 1:size(trial_factors,2)
        subplot(nrow,ncol,(i-1)*ncol+3)
        plot(trialvec([1 end]),[0 0],'-k','linewidth',2)
        hold on
        scatter(trialvec,trial_factors(:,i),10,trial_colors,'filled')
        plot(trialvec,smooth_data(trial_factors(:,i),1,15),'k')
        set(gca,'box','off','tickdir','out')
    end
