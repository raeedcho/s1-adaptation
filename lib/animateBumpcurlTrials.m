function animateBumpcurlTrials(trial_data,field,bump_colors,waitperiod,tail_length)
td = trial_data;
epoch_names = {'BL','AD','WO'};
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
            plot(plot_field(tail,1),plot_field(tail,2),'-','linewidth',2,'Color',bump_colors(bumpDir_idx,:))
            plot(plot_field(tail(end),1),plot_field(tail(end),2),'o','linewidth',2,'Color',bump_colors(bumpDir_idx,:))
        end
        axis equal
    end
    
    pause(waitperiod)
    
end