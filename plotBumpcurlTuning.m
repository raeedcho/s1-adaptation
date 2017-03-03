function plotBumpcurlTuning(trial_data)

[~,td] = getTDidx(trial_data,'result','R');

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


% get only times after bumps (120 ms)
% go to double bins to avoid edge case in trialAverage
nbins = 12;
td_unav = truncateAndBin(td,4,{'idx_bumpTime',0},{'idx_bumpTime',12});
td = truncateAndBin(td,nbins,{'idx_bumpTime',0},{'idx_bumpTime',nbins*2});

% trial average
% td_unbin = trialAverage(td_unbin,{'bumpDir','epoch'});
td = trialAverage(td,{'bumpDir','epoch'});

% plot kinematics
% figure
% hold all
% for trialidx = 1:length(td_unbin)
%     plot(td_unbin(trialidx).vel(:,1),td_unbin(trialidx).vel(:,2))
% end

% tuning curves are the first bin of S1 data divided by the bin size
for epochctr = 1:num_epochs
    [~,td_epoch] = getTDidx(td,'epoch',epoch_list{epochctr});
    temp = cat(1,td_epoch.S1_spikes)/0.15;
    temp_idx = 1:2:size(temp,1); % only index first bins
    tuning_curves{epochctr} = temp(temp_idx,:);
    bins{epochctr} = pi/180*cat(1,td_epoch.bumpDir);
end

% get tuning curves from unbinned data
for epochctr = 1:num_epochs
    [~,td_epoch] = getTDidx(td_unav,'epoch',epoch_list{epochctr});
    [~,glm_info] = getModel(td_epoch,struct('model_type','glm','model_name','linmodel',...
                                            'in_signals',{'vel'},'out_signals',{'S1_spikes'}));
    PD(epochctr,:) = atan2(glm_info.b(3,:),glm_info.b(2,:));
    offset(epochctr,:) = glm_info.b(1,:);
    moddepth(epochctr,:) = sqrt(sum(glm_info.b(2:3,:).^2,1));
end

% plot dPD curves
dPD = PD-repmat(PD(1,:),num_epochs,1);
dPD(dPD>pi) = dPD(dPD>pi)-2*pi;
dPD(dPD<-pi) = dPD(dPD<-pi)+2*pi;
figure(1)
hold all

% plot results
figure(2)
bl_color = colormap('autumn');
ad_color = colormap('summer');
wo_color = colormap('winter');
bl_idx = 1;
ad_idx = floor(linspace(1,length(ad_color),4));
wo_idx = floor(linspace(1,length(wo_color),3));
epoch_colors = [bl_color(bl_idx,:);...
               ad_color(ad_idx,:);...
               wo_color(wo_idx,:)];
           
% num_epochs = 5; % redefine for now to only take the adaptation part
for neurctr = 1:length(td(1).S1_spikes)
    figure(2)
    clf
    max_rad = 0;
%     for epochctr = 1:num_epochs
%         max_rad = max(max(tuning_curves{epochctr}(:,neurctr)),max_rad);
%     end
    max_rad = max(exp(offset(:,neurctr)+moddepth(:,neurctr)));
    for epochctr = 1:num_epochs
%         plot(bins{epochctr},tuning_curves{epochctr}(:,neurctr),'linewidth',2,'color',epoch_colors(epochctr,:));
%         set(gca,'box','off','tickdir','out','xtick',[0 pi/2 pi 3*pi/2 2*pi],'xticklabel',{'0','','\pi','','2\pi'})
%         hold on
        h = polar(0,max_rad,'.');
        set(h,'markersize',1)
        hold on
%         h = polar(repmat(bins{epochctr},2,1),repmat(tuning_curves{epochctr}(:,neurctr),2,1));
        theta = linspace(-pi,pi,100);
        h = polar(theta,exp(offset(epochctr,neurctr)+moddepth(epochctr,neurctr)*cos(theta-PD(epochctr,neurctr))));
        set(h,'color',epoch_colors(epochctr,:),'linewidth',2,'markersize',10)
        h = polar(repmat(PD(epochctr,neurctr),1,2),[0 max_rad*0.8],'-o');
        set(h,'color',epoch_colors(epochctr,:),'linewidth',2,'markersize',15)
    end
%     legend(epoch_list)
    title(['Channel ' num2str(td(1).S1_unit_guide(neurctr,1)) ', Unit ' num2str(td(1).S1_unit_guide(neurctr,2))])
    
    % dPD
    figure(1)
    clf
    plot(180/pi*dPD(1:num_epochs,neurctr),'linewidth',2)
    set(gca,'box','off','tickdir','out','xtick',1:num_epochs,'xticklabel',{epoch_list{1:num_epochs}})
    title(['Channel ' num2str(td(1).S1_unit_guide(neurctr,1)) ', Unit ' num2str(td(1).S1_unit_guide(neurctr,2))])
    
    waitforbuttonpress;
end