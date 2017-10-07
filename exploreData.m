function exploreData(trial_data)


%% get reward trials only for now
[~,td] = getTDidx(trial_data,'result','R');

%% Assign colors
% get block names and bump directions in order presented
epoch_names = {'AD','BL','WO'};
bumpDir_names = {0,90,180,270};

colors = linspecer(length(epoch_names)+length(bumpDir_names));
epoch_colors = colors(1:length(epoch_names),:);
bump_colors = colors(length(epoch_names)+1:end,:);

%% Trial average by block and bump direction (100 ms before bump and 300 ms after bump)
% % truncate and average
% td = truncateAndBin(trial_data_R,{'idx_bumpTime',-10},{'idx_bumpTime',30});
% td = trialAverage(td,{'bumpDir','epoch'});
% timevec = (-0.1:0.01:0.29)';
% 
% % set up bump direction subplot numbers
% subplot_nums = [6 2 4 8];
% 
% for emgCtr = 1:length(trial_data(1).emg_names)
%     figure('Name',trial_data(1).emg_names{emgCtr}); % meta data seems to be missing from the trial averaged data
%     min_act = inf;
%     max_act = 0;
%     for epochCtr = 1:length(epoch_names)
%         for dirCtr = 1:length(bumpDir_names)
%             td_idx = (epochCtr-1)*length(bumpDir_names)+dirCtr;
%             td(td_idx);
%             subplot(3,3,subplot_nums(dirCtr))
%             hold all
%             emg_act = td(td_idx).emg(:,emgCtr);
%             plot(timevec,emg_act,'Color',epoch_colors(epochCtr,:),'linewidth',2);
%             min_act = min(min(emg_act),min_act);
%             max_act = max(max(emg_act),max_act);
%         end
%     end
%     
%     for subplotCtr = 1:length(subplot_nums)
%         subplot(3,3,subplot_nums(subplotCtr))
%         plot([0;0],[min_act;max_act],'k--','linewidth',3)
%         set(gca,'ylim',[min_act max_act])
%     end
%     legend(epoch_names)
% end
%% Try PCA on EMGs
[~,td] = getTDidx(trial_data,'result','R');
td = truncateAndBin(td,{'idx_bumpTime',-10},{'idx_bumpTime',30});
% td = trialAverage(td,'bumpDir');

% we want the means for each neurons when we project the data later
goodEMG = ~contains(td(1).emg_names,'EMG_FDS');
[td,pca_info] = getPCA(td,struct('signals','emg','signal_idx',goodEMG));

% trial average
td = trialAverage(td,{'bumpDir','epoch'});

% set up plot
epoch_names = {'BL','AD','WO'};
figure('Name','Average EMG')

for epochCtr = 1:length(epoch_names)
    subplot(1,3,epochCtr)
    hold all;
    for i = getTDidx(td,'epoch',epoch_names(epochCtr))
        bumpDir_idx = td(i).bumpDir/90+1; % only works for 4 bump files
        plot3(td(i).emg_pca(:,1),td(i).emg_pca(:,2),td(i).emg_pca(:,3),'k-','linewidth',2,'Color',bump_colors(bumpDir_idx,:))
%         plot(td(i).emg_pca(:,1),td(i).emg_pca(:,2),'k-','linewidth',2,'Color',bump_colors(bumpDir_idx,:))
    end
    axis equal
end

%% Take a look at kinematics
[~,td] = getTDidx(trial_data,'result','R');
td = truncateAndBin(td,{'idx_bumpTime',-10},{'idx_bumpTime',30});
td = trialAverage(td,{'bumpDir','epoch'});

% animateBumpcurlTrials(td,'pos',bump_colors,0.1,5);


% set up plot
epoch_names = {'BL','AD','WO'};
figure('Name','Average Velocity')

for epochCtr = 1:length(epoch_names)
    subplot(1,3,epochCtr)
    hold all;
    for i = getTDidx(td,'epoch',epoch_names(epochCtr))
        bumpDir_idx = td(i).bumpDir/90+1; % only works for 4 bump files
        plot(td(i).vel(:,1),td(i).vel(:,2),'k-','linewidth',2,'Color',bump_colors(bumpDir_idx,:))
    end
    axis equal
end

%% S1 pca
[~,td] = getTDidx(trial_data,'result','R');
td = smoothSignals(td,struct('signals','S1_spikes','sqrt_transform',true));
td = truncateAndBin(td,{'idx_bumpTime',-10},{'idx_bumpTime',15});
% td = trialAverage(td,{'bumpDir','epoch'});

% rename epochs
epoch_list = {'BL','AD1','AD2','AD3','AD4','WO1','WO2','WO3'};
num_epochs = length(epoch_list);

% split up into finer epochs
blocks{2} = getTDidx(td,'epoch','AD','range',[0 0.25]);
blocks{3} = getTDidx(td,'epoch','AD','range',[0.25 0.5]);
blocks{4} = getTDidx(td,'epoch','AD','range',[0.5 0.75]);
blocks{5} = getTDidx(td,'epoch','AD','range',[0.75 1]);
blocks{6} = getTDidx(td,'epoch','WO','range',[0 0.33]);
blocks{7} = getTDidx(td,'epoch','WO','range',[0.33 0.67]);
blocks{8} = getTDidx(td,'epoch','WO','range',[0.67 1]);

for i = 2:8
    [td(blocks{i}).epoch] = deal(epoch_list{i});
end

% find idx of sorted units
sort_idx = find(td(1).S1_unit_guide(:,2));
[td,pca_info] = getPCA(td,struct('signals',{{'S1_spikes',sort_idx}}));

% trial average
td = trialAverage(td,{'bumpDir','epoch'});

% set up plot
figure('Name','Average S1 PC')

for epochCtr = 1:length(epoch_list)
    subplot(1,length(epoch_list),epochCtr)
    hold all;
    for i = getTDidx(td,'epoch',epoch_list(epochCtr))
        bumpDir_idx = td(i).bumpDir/90+1; % only works for 4 bump files
        %plot3(td(i).S1_pca(20,1),td(i).S1_pca(20,2),td(i).S1_pca(20,3),'o-','linewidth',2,'Color',bump_colors(bumpDir_idx,:))
        plot3(td(i).S1_pca(:,1),td(i).S1_pca(:,2),td(i).S1_pca(:,3),'-','linewidth',2,'Color',bump_colors(bumpDir_idx,:))
%         plot(td(i).S1_pca(:,1),td(i).S1_pca(:,2),'k-','linewidth',2,'Color',bump_colors(bumpDir_idx,:))
    end
    axis equal
end

% animateBumpcurlTrials(td,'S1_pca',bump_colors,0.1,60);

%% Try dPCA on S1 spikes with conditions of targets and learning block
[~,td] = getTDidx(trial_data,'result','R');
% get truncated and smoothed signals around bump time (100 ms before bump to 300 ms after)
td = smoothSignals(td,struct('signals',{'S1_spikes'},'sqrt_transform',true,'do_smoothing',true,'kernel_SD',0.1));
td = truncateAndBin(td,{'idx_bumpTime',-10},{'idx_bumpTime',30});

% get dPCA conditions
blocks{1} = getTDidx(td,'epoch','BL','range',[0 0.25]);
blocks{2} = getTDidx(td,'epoch','BL','range',[0.25 0.5]);
blocks{3} = getTDidx(td,'epoch','BL','range',[0.5 0.75]);
blocks{4} = getTDidx(td,'epoch','AD','range',[0 0.25]);
blocks{5} = getTDidx(td,'epoch','AD','range',[0.25 0.5]);
blocks{6} = getTDidx(td,'epoch','AD','range',[0.5 0.75]);
% blocks{5} = getTDidx(td,'epoch','AD','range',[0.75 1]);
% blocks{5} = getTDidx(td,'epoch','WO','range',[0 0.33]);
% blocks{6} = getTDidx(td,'epoch','WO','range',[0.33 0.67]);
% blocks{7} = getTDidx(td,'epoch','WO','range',[0.67 1]);

% get actual dPCA
sort_idx = find(td(1).S1_unit_guide(:,2));
td = getDPCA(td,blocks,'bumpDir',struct('signals',{{'S1_spikes',sort_idx}},'do_plot',true,'num_dims',10));
% 
%% check out behavior
[~,td] = getTDidx(trial_data,'result','R');

td = getSpeed(td);
td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime'));
% td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','which_method','thresh','s_thresh',20));

td = removeBadTrials(td);

metric = getLearningMetrics(td,struct('which_metric','angle','vel_or_pos','pos','time_window',{{'idx_movement_on',0;'idx_peak_speed',0}}));
% metric = getLearningMetrics(td,struct('which_metric','curvature'));

bl_idx = getTDidx(td,'epoch','BL');
ad_idx = getTDidx(td,'epoch','AD');
wo_idx = getTDidx(td,'epoch','WO');

figure
% plot([metric(bl_idx);metric(wo_idx)])
plot(metric)
hold on
plot(repmat(ad_idx(1),2,1),[-1 1],'k--','linewidth',3)
plot(repmat(wo_idx(1),2,1),[-1 1],'k--','linewidth',3)
plot([0 wo_idx(end)],[0 0],'k-','linewidth',2)

%% Plot reaches
[~,td] = getTDidx(trial_data,'result','R');

td = getSpeed(td);
td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime'));

td1 = trimTD(td,{'idx_movement_on',0},{'idx_peak_speed',0});
td2 = trimTD(td,{'idx_movement_on',0},{'idx_endTime',0});

epochs = {'BL','AD','WO'};

% [~,td1bl] = getTDidx(td1,'epoch','BL');
% [~,td1ad] = getTDidx(td1,'epoch','AD');
% [~,td1wo] = getTDidx(td1,'epoch','WO');
% [~,td2bl] = getTDidx(td2,'epoch','BL');
% [~,td2ad] = getTDidx(td2,'epoch','AD');
% [~,td2wo] = getTDidx(td2,'epoch','WO');
% 
% pos1{1} = cat(1,td1bl.pos);
% pos2{1} = cat(1,td2bl.pos);
% pos1{2} = cat(1,td1ad.pos);
% pos2{2} = cat(1,td2ad.pos);
% pos1{3} = cat(1,td1wo.pos);
% pos2{3} = cat(1,td2wo.pos);
colors = {'b','r','g'};

figure
for i=1:3
    subplot(1,2,1)
    [~,td_temp] = getTDidx(td1,'epoch',epochs{i});
    for j=1:length(td_temp)
        plot(td_temp(j).pos(:,1),td_temp(j).pos(:,2),colors{i})
        hold on
    end
    set(gca,'box','off','tickdir','out','xlim',[-15 15],'ylim',[-45,-15])
    subplot(1,2,2)
    [~,td_temp] = getTDidx(td2,'epoch',epochs{i});
    for j=1:length(td_temp)
        plot(td_temp(j).pos(:,1),td_temp(j).pos(:,2),colors{i})
        hold on
    end
    set(gca,'box','off','tickdir','out','xlim',[-15 15],'ylim',[-45,-15])
end

%% Try mahalanobis distance on time-spread PCA
% base_clusters = cell(8,1);
% for i = 1:8
%     base_clusters{i} = base_scores(base_table.bumpDir/45+1==i,1:12);
% end
% 
% adapt_clusters = cell(8,1);
% for i = 1:8
%     adapt_clusters{i} = adapt_scores(adapt_table.bumpDir/45+1==i,1:12);
% end
% 
% wash_clusters = cell(8,1);
% for i = 1:8
%     wash_clusters{i} = wash_scores(wash_table.bumpDir/45+1==i,1:12);
% end
% 
% min_mahal = inf;
% max_mahal = 0;
% 
% base_mahal = cell(8,1);
% for j = 1:8
%     for trialCtr = 1:size(base_clusters{j},1)
%         for i = 1:8
%             if i==j
%                 base_mahal{j}(trialCtr,i) = mahal(base_clusters{j}(trialCtr,1:12),base_clusters{i}([1:trialCtr-1 trialCtr+1:end],:));
%             else
%                 base_mahal{j}(trialCtr,i) = mahal(base_clusters{j}(trialCtr,1:12),base_clusters{i});
%             end
%             if min_mahal>base_mahal{j}(trialCtr,i)
%                 min_mahal = base_mahal{j}(trialCtr,i);
%             end
%             if max_mahal<base_mahal{j}(trialCtr,i)
%                 max_mahal = base_mahal{j}(trialCtr,i);
%             end
%         end
%     end
% end
% 
% adapt_mahal = cell(8,1);
% for j = 1:8
%     for trialCtr = 1:size(adapt_clusters{j},1)
%         for i = 1:8
%             adapt_mahal{j}(trialCtr,i) = mahal(adapt_clusters{j}(trialCtr,1:12),base_clusters{i});
%         end
%         if min_mahal>adapt_mahal{j}(trialCtr,i)
%             min_mahal = adapt_mahal{j}(trialCtr,i);
%         end
%         if max_mahal<adapt_mahal{j}(trialCtr,i)
%             max_mahal = adapt_mahal{j}(trialCtr,i);
%         end
%     end
% end
% wash_mahal = cell(8,1);
% for j = 1:8
%     for trialCtr = 1:size(wash_clusters{j},1)
%         for i = 1:8
%             wash_mahal{j}(trialCtr,i) = mahal(wash_clusters{j}(trialCtr,1:12),base_clusters{i});
%         end
%         if min_mahal>wash_mahal{j}(trialCtr,i)
%             min_mahal = wash_mahal{j}(trialCtr,i);
%         end
%         if max_mahal<wash_mahal{j}(trialCtr,i)
%             max_mahal = wash_mahal{j}(trialCtr,i);
%         end
%     end
% end
% 
% figure
% for i = 1:8
%     subplot(8,3,3*(i-1)+1)
%     imagesc(1./base_mahal{i}',[1/max_mahal 1/min_mahal])
%     colorbar
%     subplot(8,3,3*(i-1)+2)
%     imagesc(1./adapt_mahal{i}',[1/max_mahal 1/min_mahal])
%     colorbar
%     subplot(8,3,3*(i-1)+3)
%     imagesc(1./wash_mahal{i}',[1/max_mahal 1/min_mahal])
%     colorbar
% end
% 
% 
% 
% %% Calculate principle angles
% [base_coeff,base_scores,base_scree] = pca(base_bumps);
% [adapt_coeff,adapt_scores,~] = pca(adapt_bumps);
% [wash_coeff,wash_scores,~] = pca(wash_bumps);
% 
% adapt_base_angles = principle_angles(base_coeff(:,1:12),wash_coeff(:,1:12));
% 
% figure
% plot(adapt_base_angles*180/pi)
