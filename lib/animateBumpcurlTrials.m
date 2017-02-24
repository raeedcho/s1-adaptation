function animateBumpcurlTrials(trial_data,field,bump_colors,waitperiod,tail_length)
td = trial_data;
epoch_names = {'BL','AD','WO'};

min_field = min(cat(1,td.(field)));
max_field = max(cat(1,td.(field)));

figure;

for timeCtr = 1:size(td(1).(field),1)
    for epochCtr = 1:length(epoch_names)
        subplot(1,3,epochCtr)
        cla
        hold all;
        for i = getTDidx(td,'epoch',epoch_names(epochCtr))
            bumpDir_idx = td(i).bumpDir/90+1; % only works for 4 bump files

            %extract field
            plot_field = td(i).(field);

            tail = timeCtr-tail_length:timeCtr;
            tail = max(1,tail);
    %         plot3(td(i).emg_pca(:,1),td(i).emg_pca(:,2),td(i).emg_pca(:,3),'k-','linewidth',2,'Color',bump_colors(bumpDir_idx,:))
            plot(plot_field(tail,1),plot_field(tail,3),'-','linewidth',2,'Color',bump_colors(bumpDir_idx,:))
            plot(plot_field(tail(end),1),plot_field(tail(end),3),'o','linewidth',2,'Color',bump_colors(bumpDir_idx,:))
        end
%         axis equal
        lims = [min(min_field(1:3)) max(max_field(1:3))];
        set(gca,'xlim',lims,'ylim',lims)
    end
    
    pause(waitperiod)
    
end
