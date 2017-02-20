%% Add stuff to path
datecode = '20170206';

if(ispc)
    homeFolder = 'C:\Users\rhc307\';
else
    homeFolder = '/home/raeed/';
end
addpath(genpath([homeFolder filesep 'Projects' filesep 'limblab' filesep 'ClassyDataAnalysis']))
addpath([homeFolder filesep 'Projects' filesep 'limblab' filesep 'proc-raeed' filesep 'BumpCurl' filesep 'lib' filesep])
addpath(genpath([homeFolder filesep 'Projects' filesep 'limblab' filesep 'proc-raeed' filesep 'lib']))
cd([homeFolder 'Projects' filesep 'limblab' filesep 'data-raeed' filesep 'BumpCurl' filesep 'Han' filesep datecode filesep])

clear homeFolder


%% Load CDS files

% [base_cds,adapt_cds,wash_cds] = loadsave_bumpcurl_cds(datecode);

%% Extract data from cds files
clear params;
params.array_alias = {'LeftS1Area2','S1'};
params.exclude_units = 255;
params.event_list = {'bumpTime';'bumpDir'};
meta = struct('force_direction',1.48,'block','BL');
params.meta = meta;
trial_data_base = parseFileByTrial(base_cds,params);

params.meta.block = 'AD';
trial_data_adapt = parseFileByTrial(adapt_cds,params);

params.meta.block = 'WO';
trial_data_wash = parseFileByTrial(wash_cds,params);

trial_data = cat(2,trial_data_base,trial_data_adapt,trial_data_wash);

%% Check trial tables
% check reward tables
% base_table = base_cds.trials(base_cds.trials.result=='R',:);
% adapt_table = adapt_cds.trials(adapt_cds.trials.result=='R',:);
% wash_table = wash_cds.trials(wash_cds.trials.result=='R',:);
% base_table = base_table(base_cds.trials.result=='R',:);
% adapt_table = adapt_table(adapt_cds.trials.result=='R',:);
% wash_table = wash_table(wash_cds.trials.result=='R',:);

% clear all trials with no center hold bump
base_table = base_table(base_table.ctrHoldBump==1 & ~isnan(base_table.bumpTime),:);
adapt_table = adapt_table(adapt_table.ctrHoldBump==1 & ~isnan(adapt_table.bumpTime),:);
wash_table = wash_table(wash_table.ctrHoldBump==1 & ~isnan(wash_table.bumpTime),:);


%% Check if extra stuff is lined up with neural stuff
figure
subplot(311)
plot(base_data.t,base_data.KinectSyncPulse)
hold on
plot(extra_base_data.t,extra_base_data.Sync)
subplot(312)
plot(adapt_data.t,adapt_data.KinectSyncPulse)
hold on
plot(extra_adapt_data.t,extra_adapt_data.Sync)
subplot(313)
plot(wash_data.t,wash_data.KinectSyncPulse)
hold on
plot(extra_wash_data.t,extra_wash_data.Sync)


%% Concatenate data tables
interp_extra_base = array2table(interp1(extra_base_data.t,extra_base_data{:,2:end},base_data.t),'VariableNames',extra_base_data.Properties.VariableNames(2:end));
interp_extra_adapt = array2table(interp1(extra_adapt_data.t,extra_adapt_data{:,2:end},adapt_data.t),'VariableNames',extra_adapt_data.Properties.VariableNames(2:end));
interp_extra_wash = array2table(interp1(extra_wash_data.t,extra_wash_data{:,2:end},wash_data.t),'VariableNames',extra_wash_data.Properties.VariableNames(2:end));

base_data = [base_data interp_extra_base];
adapt_data = [adapt_data interp_extra_adapt];
wash_data = [wash_data interp_extra_wash];

clear interp_*
clear extra_*_data

%% Concatenate EMGs
interp_extra_base = array2table(interp1(extra_base_emg.t,extra_base_emg{:,2:end},base_emg.t),'VariableNames',extra_base_emg.Properties.VariableNames(2:end));
interp_extra_adapt = array2table(interp1(extra_adapt_emg.t,extra_adapt_emg{:,2:end},adapt_emg.t),'VariableNames',extra_adapt_emg.Properties.VariableNames(2:end));
interp_extra_wash = array2table(interp1(extra_wash_emg.t,extra_wash_emg{:,2:end},wash_emg.t),'VariableNames',extra_wash_emg.Properties.VariableNames(2:end));

base_emg = [base_emg interp_extra_base];
adapt_emg = [adapt_emg interp_extra_adapt];
wash_emg = [wash_emg interp_extra_wash];

clear interp_*
clear extra_*_emg

%% Try PCA on data
neurIDX = contains(base_data.Properties.VariableNames,'LeftS1Area2CH');
base_bumps = [];
for trialCtr=1:height(base_table)
    base_bumps = [base_bumps; base_data{base_data.t>=base_table.bumpTime(trialCtr) & base_data.t<=base_table.bumpTime(trialCtr)+0.3,neurIDX}];
end

neurIDX = contains(adapt_data.Properties.VariableNames,'LeftS1Area2CH');
adapt_bumps = [];
for trialCtr=1:height(adapt_table)
    adapt_bumps = [adapt_bumps; adapt_data{adapt_data.t>=adapt_table.bumpTime(trialCtr) & adapt_data.t<=adapt_table.bumpTime(trialCtr)+0.3,neurIDX}];
end

neurIDX = contains(wash_data.Properties.VariableNames,'LeftS1Area2CH');
wash_bumps = [];
for trialCtr=1:height(wash_table)
    wash_bumps = [wash_bumps; wash_data{wash_data.t>=wash_table.bumpTime(trialCtr) & wash_data.t<=wash_table.bumpTime(trialCtr)+0.3,neurIDX}];
end

[coeff,scores,latent_eigs] = pca(sqrt([base_bumps;adapt_bumps;wash_bumps]));
% [base_coeff,base_scores,base_scree] = pca([base_bumps]);

% project adapt and wash bumps into pc space
% center
% adapt_bumps_ctr = adapt_bumps - mean(adapt_bumps);
% wash_bumps_ctr = wash_bumps - mean(wash_bumps);
% 
% % project
% adapt_scores = adapt_bumps_ctr*base_coeff;
% wash_scores = wash_bumps_ctr*base_coeff;
% 


base_idx = (1:size(base_bumps,1))';
adapt_idx = (base_idx(end)+1:base_idx(end)+size(adapt_bumps,1))';
wash_idx = (adapt_idx(end)+1:adapt_idx(end)+size(wash_bumps,1))';

base_scores = scores(base_idx,:);
adapt_scores = scores(adapt_idx,:);
wash_scores = scores(wash_idx,:);

base_scree = var(base_scores);
adapt_scree = var(adapt_scores);
wash_scree = var(wash_scores);

figure
subplot(311)
bar(base_scree(1:20))
subplot(312)
bar(adapt_scree(1:20))
subplot(313)
bar(wash_scree(1:20))


%% Look at PCs time traces


%% Plot color bumps
num_bumps = 4;

bump_sep = 360/num_bumps;

colors = linspecer(num_bumps);
base_colors = colors(floor(base_table.bumpDir/bump_sep+1),:);
adapt_colors = colors(floor(adapt_table.bumpDir/bump_sep+1),:);
wash_colors = colors(floor(wash_table.bumpDir/bump_sep+1),:);

figure
S = linspace(1,200,size(base_scores,1));
for j = 1:6
    subplot(2,3,j)
    for i = 1:size(base_colors,1)
%         idx = (i-1)*6+1:i*6;
%         idx = i*6-3;
        idx = (i-1)*6+j;
%         plot3(base_scores(idx,1),base_scores(idx,2),base_scores(idx,3),'Color',base_colors(i,:))
%         plot(base_scores(idx,1),base_scores(idx,2),'o','Color',base_colors(i,:))
        hold on
        scatter3(base_scores(idx,1),base_scores(idx,2),base_scores(idx,3),S(i),base_colors(i,:))
    end
    grid off
    axis square
end


figure
S = linspace(1,200,size(adapt_scores,1));
for j = 1:6
    subplot(2,3,j)
    for i = 1:size(adapt_colors,1)
        
%         idx = (i-1)*6+1:i*6;
%         idx = i*6-3;
        idx = (i-1)*6+j;
%         plot3(adapt_scores(idx,1),adapt_scores(idx,2),adapt_scores(idx,3),'Color',adapt_colors(i,:))
%         plot(adapt_scores(idx,1),adapt_scores(idx,2),'o','Color',adapt_colors(i,:))
        hold on
%         scatter(adapt_scores(idx,1),adapt_scores(idx,2),S(i),adapt_colors(i,:))
        scatter3(adapt_scores(idx,1),adapt_scores(idx,2),adapt_scores(idx,3),S(i),adapt_colors(i,:))
        
    end
    grid off
    axis square
end

figure
S = linspace(1,200,size(wash_scores,1));
for j = 1:6
    subplot(2,3,j)
    for i = 1:size(wash_colors,1)
%         idx = (i-1)*6+1:i*6;
%         idx = i*6-3;
        idx = (i-1)*6+j;
%         plot3(wash_scores(idx,1),wash_scores(idx,2),wash_scores(idx,3),'Color',wash_colors(i,:))
%         plot(wash_scores(idx,1),wash_scores(idx,2),'o','Color',wash_colors(i,:))
        hold on
%         scatter(wash_scores(idx,1),wash_scores(idx,2),S(i),wash_colors(i,:))
        scatter3(wash_scores(idx,1),wash_scores(idx,2),wash_scores(idx,3),S(i),wash_colors(i,:))
    end
    grid off
    axis square
end



%% Try mahalanobis distance on time-spread PCA
base_clusters = cell(8,1);
for i = 1:8
    base_clusters{i} = base_scores(base_table.bumpDir/45+1==i,1:12);
end

adapt_clusters = cell(8,1);
for i = 1:8
    adapt_clusters{i} = adapt_scores(adapt_table.bumpDir/45+1==i,1:12);
end

wash_clusters = cell(8,1);
for i = 1:8
    wash_clusters{i} = wash_scores(wash_table.bumpDir/45+1==i,1:12);
end

min_mahal = inf;
max_mahal = 0;

base_mahal = cell(8,1);
for j = 1:8
    for trialCtr = 1:size(base_clusters{j},1)
        for i = 1:8
            if i==j
                base_mahal{j}(trialCtr,i) = mahal(base_clusters{j}(trialCtr,1:12),base_clusters{i}([1:trialCtr-1 trialCtr+1:end],:));
            else
                base_mahal{j}(trialCtr,i) = mahal(base_clusters{j}(trialCtr,1:12),base_clusters{i});
            end
            if min_mahal>base_mahal{j}(trialCtr,i)
                min_mahal = base_mahal{j}(trialCtr,i);
            end
            if max_mahal<base_mahal{j}(trialCtr,i)
                max_mahal = base_mahal{j}(trialCtr,i);
            end
        end
    end
end

adapt_mahal = cell(8,1);
for j = 1:8
    for trialCtr = 1:size(adapt_clusters{j},1)
        for i = 1:8
            adapt_mahal{j}(trialCtr,i) = mahal(adapt_clusters{j}(trialCtr,1:12),base_clusters{i});
        end
        if min_mahal>adapt_mahal{j}(trialCtr,i)
            min_mahal = adapt_mahal{j}(trialCtr,i);
        end
        if max_mahal<adapt_mahal{j}(trialCtr,i)
            max_mahal = adapt_mahal{j}(trialCtr,i);
        end
    end
end
wash_mahal = cell(8,1);
for j = 1:8
    for trialCtr = 1:size(wash_clusters{j},1)
        for i = 1:8
            wash_mahal{j}(trialCtr,i) = mahal(wash_clusters{j}(trialCtr,1:12),base_clusters{i});
        end
        if min_mahal>wash_mahal{j}(trialCtr,i)
            min_mahal = wash_mahal{j}(trialCtr,i);
        end
        if max_mahal<wash_mahal{j}(trialCtr,i)
            max_mahal = wash_mahal{j}(trialCtr,i);
        end
    end
end

figure
for i = 1:8
    subplot(8,3,3*(i-1)+1)
    imagesc(1./base_mahal{i}',[1/max_mahal 1/min_mahal])
    colorbar
    subplot(8,3,3*(i-1)+2)
    imagesc(1./adapt_mahal{i}',[1/max_mahal 1/min_mahal])
    colorbar
    subplot(8,3,3*(i-1)+3)
    imagesc(1./wash_mahal{i}',[1/max_mahal 1/min_mahal])
    colorbar
end



%% Calculate principle angles
[base_coeff,base_scores,base_scree] = pca(base_bumps);
[adapt_coeff,adapt_scores,~] = pca(adapt_bumps);
[wash_coeff,wash_scores,~] = pca(wash_bumps);

adapt_base_angles = principle_angles(base_coeff(:,1:12),wash_coeff(:,1:12));

figure
plot(adapt_base_angles*180/pi)


%% Look at behavior
% kinIDX = [4 5];
kinIDX = [10 11];
base_bumps = [];
for trialCtr=1:height(base_table)
    base_bumps = [base_bumps; base_data{base_data.t>=base_table.bumpTime(trialCtr) & base_data.t<=base_table.bumpTime(trialCtr)+0.3,kinIDX}];
end

adapt_bumps = [];
for trialCtr=1:height(adapt_table)
    adapt_bumps = [adapt_bumps; adapt_data{adapt_data.t>=adapt_table.bumpTime(trialCtr) & adapt_data.t<=adapt_table.bumpTime(trialCtr)+0.3,kinIDX}];
end

wash_bumps = [];
for trialCtr=1:height(wash_table)
    wash_bumps = [wash_bumps; wash_data{wash_data.t>=wash_table.bumpTime(trialCtr) & wash_data.t<=wash_table.bumpTime(trialCtr)+0.3,kinIDX}];
end

num_bumps = 4;

bump_sep = 360/num_bumps;

colors = linspecer(num_bumps);
base_colors = colors(floor(base_table.bumpDir/bump_sep+1),:);
adapt_colors = colors(floor(adapt_table.bumpDir/bump_sep+1),:);
wash_colors = colors(floor(wash_table.bumpDir/bump_sep+1),:);

figure
S = linspace(1,40,size(base_bumps,1));
for j = 1:6
    subplot(2,3,j)
    for i = 1:size(base_colors,1)
%         idx = (i-1)*6+1:i*6;
%         idx = i*6-3;
        idx = (i-1)*6+j;
%         plot3(base_scores(idx,1),base_scores(idx,2),base_scores(idx,3),'Color',base_colors(i,:))
%         plot(base_scores(idx,1),base_scores(idx,2),'o','Color',base_colors(i,:))
        hold on
        scatter(base_bumps(idx,1),base_bumps(idx,2),S(i),base_colors(i,:))
    end
    grid off
%     axis([-30 30 -30 30])
%     axis([-5 5 -40 -30])
end


figure
S = linspace(1,40,size(adapt_bumps,1));
for j = 1:6
    subplot(2,3,j)
    for i = size(adapt_colors,1)-150:size(adapt_colors,1)
        
%         idx = (i-1)*6+1:i*6;
%         idx = i*6-3;
        idx = (i-1)*6+j;
%         plot3(adapt_scores(idx,1),adapt_scores(idx,2),adapt_scores(idx,3),'Color',adapt_colors(i,:))
%         plot(adapt_scores(idx,1),adapt_scores(idx,2),'o','Color',adapt_colors(i,:))
        hold on
        scatter(adapt_bumps(idx,1),adapt_bumps(idx,2),S(i),adapt_colors(i,:))
%         scatter3(adapt_scores(idx,1),adapt_scores(idx,2),adapt_scores(idx,3),S(i),adapt_colors(i,:))
        
    end
    grid off
%     axis([-30 30 -30 30])
%     axis([-5 5 -40 -30])
end

figure
S = linspace(1,40,size(wash_bumps,1));
for j = 1:6
    subplot(2,3,j)
    for i = 1:size(wash_colors,1)
%         idx = (i-1)*6+1:i*6;
%         idx = i*6-3;
        idx = (i-1)*6+j;
%         plot3(wash_scores(idx,1),wash_scores(idx,2),wash_scores(idx,3),'Color',wash_colors(i,:))
%         plot(wash_scores(idx,1),wash_scores(idx,2),'o','Color',wash_colors(i,:))
        hold on
        scatter(wash_bumps(idx,1),wash_bumps(idx,2),S(i),wash_colors(i,:))
        
    end
    grid off
%     axis([-30 30 -30 30])
%     axis([-5 5 -40 -30])
end


%% Look at EMG PEH
base_times = cell(1,4);
adapt_times = cell(1,4);
wash_times = cell(1,4);
for bumpCtr = 1:4
    base_times{bumpCtr} = base_table.bumpTime(base_table.bumpDir==(bumpCtr-1)*90);
    adapt_times{bumpCtr} = adapt_table.bumpTime(adapt_table.bumpDir==(bumpCtr-1)*90);
    wash_times{bumpCtr} = wash_table.bumpTime(wash_table.bumpDir==(bumpCtr-1)*90);
end


num_bumps = 4;

bump_sep = 360/num_bumps;

colors = linspecer(num_bumps);

for emgCtr = 2:width(base_emg)
    figure
    title(base_emg.Properties.VariableNames{emgCtr});
    xwindow = [-0.1 0.3];
    for i = 1:4
        [base_av,baseCIlow,baseCIhigh] = getPEH(base_emg.t,base_emg{:,emgCtr},base_times{i},xwindow);
        [adapt_av,adaptCIlow,adaptCIhigh] = getPEH(adapt_emg.t,adapt_emg{:,emgCtr},adapt_times{i},xwindow);
        [wash_av,washCIlow,washCIhigh] = getPEH(wash_emg.t,wash_emg{:,emgCtr},wash_times{i},xwindow);
        
        minimax(i,:) = [min([baseCIlow adaptCIlow washCIlow]) max([baseCIhigh adaptCIhigh washCIhigh])];
        
        timevec = linspace(-0.1,0.3,length(base_av));
        subplot(1,3,1)
        patch([timevec fliplr(timevec)],[baseCIlow fliplr(baseCIhigh)],colors(i,:))
%         plot(timevec,base_av,'color',colors(i,:))
        hold on
        set(gca,'ylim',ywindow)
        subplot(1,3,2)
        patch([timevec fliplr(timevec)],[adaptCIlow fliplr(adaptCIhigh)],colors(i,:))
%         plot(timevec,adapt_av,'color',colors(i,:))
        title(base_emg.Properties.VariableNames{emgCtr});
        hold on
        set(gca,'ylim',ywindow)
        subplot(1,3,3)
        patch([timevec fliplr(timevec)],[washCIlow fliplr(washCIhigh)],colors(i,:))
%         plot(timevec,wash_av,'color',colors(i,:))
        hold on
        
    end
    ywindow(1) = min(minimax(:,1));
    ywindow(2) = max(minimax(:,2));
    for i = 1:3
        subplot(1,3,i)
        set(gca,'ylim',ywindow)
    end
end

%% Look at EMG PCA
EMGIDX = contains(base_emg.Properties.VariableNames,'EMG');
base_bumps = [];
for trialCtr=1:height(base_table)
    base_bumps = [base_bumps; base_emg{base_emg.t>=base_table.bumpTime(trialCtr)-0.1 & base_emg.t<=base_table.bumpTime(trialCtr)+0.3,EMGIDX}];
end

adapt_bumps = [];
for trialCtr=1:height(adapt_table)
    adapt_bumps = [adapt_bumps; adapt_emg{adapt_emg.t>=adapt_table.bumpTime(trialCtr)-0.1 & adapt_emg.t<=adapt_table.bumpTime(trialCtr)+0.3,EMGIDX}];
end

wash_bumps = [];
for trialCtr=1:height(wash_table)
    wash_bumps = [wash_bumps; wash_emg{wash_emg.t>=wash_table.bumpTime(trialCtr)-0.1 & wash_emg.t<=wash_table.bumpTime(trialCtr)+0.3,EMGIDX}];
end

window = [-0.1 0.3];
base_bumps = catDataWindow(base_emg.t,base_emg{:,EMGIDX},base_table.bumpTime,window);
adapt_bumps = catDataWindow(adapt_emg.t,adapt_emg{:,EMGIDX},adapt_table.bumpTime,window);
wash_bumps = catDataWindow(wash_emg.t,wash_emg{:,EMGIDX},wash_table.bumpTime,window);


[coeff,scores,latent_eigs] = pca([base_bumps;adapt_bumps;wash_bumps]);
% [base_coeff,base_scores,base_scree] = pca([base_bumps]);

% project adapt and wash bumps into pc space
% center
% adapt_bumps_ctr = adapt_bumps - mean(adapt_bumps);
% wash_bumps_ctr = wash_bumps - mean(wash_bumps);
% 
% % project
% adapt_scores = adapt_bumps_ctr*base_coeff;
% wash_scores = wash_bumps_ctr*base_coeff;
% 

% do pca
base_idx = (1:size(base_bumps,1))';
adapt_idx = (base_idx(end)+1:base_idx(end)+size(adapt_bumps,1))';
wash_idx = (adapt_idx(end)+1:adapt_idx(end)+size(wash_bumps,1))';

base_scores = scores(base_idx,:);
adapt_scores = scores(adapt_idx,:);
wash_scores = scores(wash_idx,:);

base_scree = var(base_scores);
adapt_scree = var(adapt_scores);
wash_scree = var(wash_scores);

figure
subplot(311)
bar(base_scree(1:20))
subplot(312)
bar(adapt_scree(1:20))
subplot(313)
bar(wash_scree(1:20))

% plot bumps
num_bumps = 4;

bump_sep = 360/num_bumps;

colors = linspecer(num_bumps);
base_colors = colors(floor(base_table.bumpDir/(pi/2)+1),:);
adapt_colors = colors(floor(adapt_table.bumpDir/bump_sep+1),:);
wash_colors = colors(floor(wash_table.bumpDir/bump_sep+1),:);



% figure
% S = linspace(1,40,size(base_scores,1));
% for j = 1:6
%     subplot(2,3,j)
%     for i = 1:size(base_colors,1)
% %         idx = (i-1)*6+1:i*6;
% %         idx = i*6-3;
%         idx = (i-1)*6+j;
% %         plot3(base_scores(idx,1),base_scores(idx,2),base_scores(idx,3),'Color',base_colors(i,:))
% %         plot(base_scores(idx,1),base_scores(idx,2),'o','Color',base_colors(i,:))
%         hold on
%         scatter(base_scores(idx,1),base_scores(idx,2),S(i),base_colors(i,:))
%     end
%     grid off
% %     axis([-30 30 -30 30])
% %     axis([-5 5 -40 -30])
% end
% 
% 
% figure
% S = linspace(1,40,size(adapt_scores,1));
% for j = 1:6
%     subplot(2,3,j)
%     for i = size(adapt_colors,1)-150:size(adapt_colors,1)
%         
% %         idx = (i-1)*6+1:i*6;
% %         idx = i*6-3;
%         idx = (i-1)*6+j;
% %         plot3(adapt_scores(idx,1),adapt_scores(idx,2),adapt_scores(idx,3),'Color',adapt_colors(i,:))
% %         plot(adapt_scores(idx,1),adapt_scores(idx,2),'o','Color',adapt_colors(i,:))
%         hold on
%         scatter(adapt_scores(idx,1),adapt_scores(idx,2),S(i),adapt_colors(i,:))
% %         scatter3(adapt_scores(idx,1),adapt_scores(idx,2),adapt_scores(idx,3),S(i),adapt_colors(i,:))
%         
%     end
%     grid off
% %     axis([-30 30 -30 30])
% %     axis([-5 5 -40 -30])
% end
% 
% figure
% S = linspace(1,40,size(wash_scores,1));
% for j = 1:6
%     subplot(2,3,j)
%     for i = 1:size(wash_colors,1)
% %         idx = (i-1)*6+1:i*6;
% %         idx = i*6-3;
%         idx = (i-1)*6+j;
% %         plot3(wash_scores(idx,1),wash_scores(idx,2),wash_scores(idx,3),'Color',wash_colors(i,:))
% %         plot(wash_scores(idx,1),wash_scores(idx,2),'o','Color',wash_colors(i,:))
%         hold on
%         scatter(wash_scores(idx,1),wash_scores(idx,2),S(i),wash_colors(i,:))
%         
%     end
%     grid off
% %     axis([-30 30 -30 30])
% %     axis([-5 5 -40 -30])
% end

