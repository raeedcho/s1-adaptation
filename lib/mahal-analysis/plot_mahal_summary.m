function plot_mahal_summary(mahal_curve_table,params)

    smoothing_window_size = 10;
    num_mahal_dims = 4;
    arrays_to_plot = {'M1','S1'};
    save_figures = false;
    figsavedir = '';
    assignParams(who,params)

    smoothing_kernel = ones(1,smoothing_window_size)/smoothing_window_size;

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
        if height(mahal_curve_array)==0
            continue
        end
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
