% %% load & process td file
clear, clc

% pathname = fullfile('~','Koofr','Limblab','Data','S1-adaptation','td');
filename = 'Duncan_20190911_CObumpcurl_all_5ms_spindle_paper.mat';
load([filename]);

%% PD Calculations
td_bump = td(~isnan([td.idx_bumpTime]));

td_pd = td;

trial_types = {td_pd.epoch};

blIdx = strcmpi(trial_types,'bl');
adIdx = strcmpi(trial_types,'ad');
woIdx = strcmpi(trial_types,'wo');


td_pd = trimTD(td_pd,{'idx_movement_on',0},{'idx_movement_on',50});

% td_pd = trimTD(td_pd,{'idx_bumpTime',0},{'idx_bumpTime',50});

td_bl = td_pd(blIdx);
td_ad = td_pd(adIdx);
td_wo = td_pd(woIdx);

nTrials = 150;

td_ad_early = td_ad(1:nTrials);
td_ad_late = td_ad(end-nTrials:end);
td_wo_early = td_wo(1:nTrials);
td_wo_late = td_wo(end-nTrials:end);

params.out_signals = {'S1_spikes_bins','all'};
params.out_signal_names = td(1).S1_unit_guide;
params.in_signals = {'vel',1:2};
params.num_boots = 100;




blRand = randperm(length(td_bl));
bl_n = length(td_bl);
idx1 = blRand(1:round(bl_n/2));
idx2 = blRand(round(bl_n/2)/2:bl_n);

bl1_pdTable = getTDPDs(td_bl(idx1),params);
bl2_pdTable = getTDPDs(td_bl(idx2),params);
bl_pdTable = getTDPDs(td_bl(:),params);
ade_pdTable = getTDPDs(td_ad_early,params);
adl_pdTable = getTDPDs(td_ad_late,params);
woe_pdTable = getTDPDs(td_wo,params);
wol_pdTable = getTDPDs(td_wo_late,params);



% ad_pdTable = getTDPDs(td_ad,params);
% wo_pdTable = getTDPDs(td_wo,params);
%% Separate PDs by which population to analyze
% popIdx = 1:size(td_pd(1).S1_spindle,2);
popIdx = 1:size(td_pd(1).(params.out_signals{1}),2);


popIdx = popIdx(bl_pdTable.velTuned ...
    & ade_pdTable.velTuned ...
    & adl_pdTable.velTuned ...
    & woe_pdTable.velTuned);



% popIdx = popIdx(...
%     abs(bl_pdTable.velPDCI(popIdx,2) - bl_pdTable.velPDCI(popIdx,1)) < pi/4 ...
%     & abs(ade_pdTable.velPDCI(popIdx,2) - ade_pdTable.velPDCI(popIdx,1)) < pi/4 ...
%     & abs(adl_pdTable.velPDCI(popIdx,2) - adl_pdTable.velPDCI(popIdx,1)) < pi/4 ...
%     & abs(woe_pdTable.velPDCI(popIdx,2) - woe_pdTable.velPDCI(popIdx,1)) < pi/4);




bl1PDs = rad2deg(bl1_pdTable.velPD(popIdx));
bl1PDs = wrapTo360(bl1PDs);
% bl1Boots = wrapTo360(rad2deg(cat(1,bl1_pdTable.velbootstraps{:})));

bl2PDs = rad2deg(bl2_pdTable.velPD(popIdx));
bl2PDs = wrapTo360(bl2PDs);
% bl2Boots = wrapTo360(rad2deg(cat(1,bl2_pdTable.velbootstraps{:})));

blPDs = rad2deg(bl_pdTable.velPD(popIdx));
blPDs = wrapTo360(blPDs);
% blBoots = wrapTo360(rad2deg(cat(1,bl_pdTable.velbootstraps{:})));

adePDs = rad2deg(ade_pdTable.velPD(popIdx));
adePDs = wrapTo360(adePDs);
% adeBoots = wrapTo360(rad2deg(cat(1,ade_pdTable.velbootstraps{:})));

adlPDs = rad2deg(adl_pdTable.velPD(popIdx));
adlPDs = wrapTo360(adlPDs);
% adlBoots = wrapTo360(rad2deg(cat(1,adl_pdTable.velbootstraps{:})));

woePDs = rad2deg(woe_pdTable.velPD(popIdx));
woePDs = wrapTo360(woePDs);
% woeBoots = wrapTo360(rad2deg(cat(1,woe_pdTable.velbootstraps{:})));

wolPDs = rad2deg(wol_pdTable.velPD(popIdx));
wolPDs = wrapTo360(wolPDs);
% wolBoots = wrapTo360(rad2deg(cat(1,wol_pdTable.velbootstraps{:})));

% woPDs = rad2deg(wo_pdTable.velPD(popIdx));
% woPDs(woPDs<0) = woPDs(woPDs<0)+360;
%
% adPDs = rad2deg(ad_pdTable.velPD(popIdx));
% adPDs(adPDs<0) = adPDs(adPDs<0)+360;
%% Make figures a la London and Miller 2013
bl1CI = wrapTo180(-bl1PDs + rad2deg(bl1_pdTable.velPDCI(popIdx,:)));
bl2CI = wrapTo180(-bl2PDs + rad2deg(bl2_pdTable.velPDCI(popIdx,:)));
blCI = wrapTo180(-blPDs + rad2deg(bl_pdTable.velPDCI(popIdx,:)));
adeCI = wrapTo180(-adePDs + rad2deg(ade_pdTable.velPDCI(popIdx,:)));
adlCI = wrapTo180(-adlPDs + rad2deg(adl_pdTable.velPDCI(popIdx,:)));
woeCI = wrapTo180(-woePDs + rad2deg(woe_pdTable.velPDCI(popIdx,:)));
% wolCI = wrapTo180(-wolPDs + rad2deg(wol_pdTable.velPDCI(popIdx,:)));
% woCI = wrapTo180(-woPDs + rad2deg(wo_pdTable.velPDCI(popIdx,:)));
% adCI = wrapTo180(-adPDs + rad2deg(ad_pdTable.velPDCI(popIdx,:)));





%bl2bl PD scatter
hfig = figure;
hold on; set(gcf,'Color','White','Units','Normalized','Position',[0.25 0.25 0.5 0.32])
% subplot(2,4,1), hold on;
% shadeX = rad2deg([0 2*pi pi 0]);
% shadeY = rad2deg([0 2*pi 2*pi pi]);
% fill(shadeX,shadeY,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
% shadeX2 = rad2deg([pi 2*pi 2*pi pi]);
% shadeY2 = rad2deg([0 0 pi 0]);
% fill(shadeX2,shadeY2,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
% 
% plot(bl1PDs,bl2PDs,'.','Color',[157 138 142]/256,'MarkerSize',15);
% errorbar(bl1PDs,bl2PDs,bl2CI(:,1),bl2CI(:,2),bl1CI(:,1),bl1CI(:,2),...
%     'LineStyle','none','Color',[157 138 142]/256);
ticks = [0 180 360];
% set(gca,'FontName','Helvetica','FontSize',12,'YTick',ticks,'XTick',ticks), axis([0 360 0 360])
% % xlabel('Baseline PDs'), ylabel('Baseline PDs')
% title('Baseline PDs','color',[157 138 142]/256)
% 
% 
%bl2bl differences
bl2blPDdiff = angleDiff(bl2PDs,bl1PDs,false,true);
% bl2blPDdiff = bl2blPDdiff(~any((bl1PDs - (bl2CI+bl2PDs))'<0)');
% bl2blPDdiff = angleDiff(bl2Boots,bl1Boots,false,true);
% subplot(2,5,6), hold on;
% shadeX3 = [0 -180 -180 0];
% shadeY3 = [0 0 length(bl1PDs) length(bl1PDs)];
% fill(shadeX3,shadeY3,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
% 
% histogram(bl2blPDdiff,-180:9:180,'FaceColor',[157 138 142]/256,'FaceAlpha',1)
% line([mean(bl2blPDdiff) mean(bl2blPDdiff)],[0 length(bl1Boots)],'Color',[0 0 0],'LineWidth',2)
% ylabel('Number of neurons')
% set(gca,'FontName','Helvetica','FontSize',12)
% axis([-90 90 0 length(bl1PDs)])


%bl2ade PD scatter
subplot(2,3,1), hold on;
shadeX = rad2deg([0 2*pi pi 0]);
shadeY = rad2deg([0 2*pi 2*pi pi]);
fill(shadeX,shadeY,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
shadeX2 = rad2deg([pi 2*pi 2*pi pi]);
shadeY2 = rad2deg([0 0 pi 0]);
fill(shadeX2,shadeY2,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)

plot(blPDs,adePDs,'.','color',[0 215 157]/256,'MarkerSize',15);
errorbar(blPDs,adePDs,adeCI(:,1),adeCI(:,2),blCI(:,1),blCI(:,2),...
    'LineStyle','none','Color',[0 215 157]/256);
set(gca,'FontName','Helvetica','FontSize',12,'YTick',ticks,'XTick',ticks), axis([0 360 0 360])
ylabel('Early Force PDs','color',[0 215 157]/256)



%bl2ade differences
bl2adePDdiff = angleDiff(adePDs,blPDs,false,true);
bl2adePDdiff_sig = bl2adePDdiff(~any((blPDs - (adeCI+adePDs))'<0)' | ~any((blPDs - (adeCI+adePDs))'>0)');
% bl2adePDdiff_notSig = bl2adePDdiff(any((blPDs - (adeCI+adePDs))'<0)' & any((blPDs - (adeCI+adePDs))'>0)');

subplot(2,3,4), hold on;
shadeX3 = [0 -180 -180 0];
shadeY3 = [0 0 length(blPDs) length(blPDs)];
fill(shadeX3,shadeY3,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)

histogram(bl2adePDdiff_sig,-90:9:90,'FaceColor',[0 215 157]/256,'FaceAlpha',1)
% histogram(bl2blPDdiff,-90:9:90,'FaceColor',[239 68 232]/256,'FaceAlpha',0.5)
% histogram(bl2adePDdiff_notSig,-90:9:90,'FaceColor',[0 0 0]/256,'FaceAlpha',0.2)

line([median(bl2adePDdiff_sig) median(bl2adePDdiff_sig)],[0 length(blPDs)],'Color',[0 215 157]/256,'LineWidth',2)
% line([median(bl2blPDdiff) median(bl2blPDdiff)],[0 length(blPDs)],'Color',[239 68 232]/256,'LineWidth',2)
% line([median(bl2adePDdiff_notSig) median(bl2adePDdiff_notSig)],[0 length(blPDs)],'Color',[0 0 0]/256,'LineWidth',0.5)

set(gca,'FontName','Helvetica','FontSize',12)
axis([-90 90 0 length(blPDs)])
text(-85,10,...
    ['n = ' num2str(length(bl2adePDdiff_sig)) '/' num2str(numel(popIdx))],...
    'FontName','Helvetica','FontSize',12)
ylabel('Number of Neurons')




%bl2adl PD scatter
subplot(2,3,2), hold on;
shadeX = rad2deg([0 2*pi pi 0]);
shadeY = rad2deg([0 2*pi 2*pi pi]);
fill(shadeX,shadeY,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
shadeX2 = rad2deg([pi 2*pi 2*pi pi]);
shadeY2 = rad2deg([0 0 pi 0]);
fill(shadeX2,shadeY2,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)

plot(blPDs,adlPDs,'.','color',[0 215 157]/256,'MarkerSize',15);
errorbar(blPDs,adlPDs,adlCI(:,1),adlCI(:,2),blCI(:,1),blCI(:,2),...
    'LineStyle','none','Color',[0 215 157]/256);
set(gca,'FontName','Helvetica','FontSize',12,'YTick',ticks,'XTick',ticks), axis([0 360 0 360])
xlabel('Baseline PDs','color',[239 68 232]/256)
%, ylabel('Adaptation (late) PDs')
ylabel('Late Force PDs','color',[0 215 157]/256)


%bl2adl differences
bl2adlPDdiff = angleDiff(adlPDs,blPDs,false,true);
bl2adlPDdiff_sig = bl2adlPDdiff(~any((blPDs - (adlCI+adlPDs))'<0)' | ~any((blPDs - (adlCI+adlPDs))'>0)');
% bl2adlPDdiff_notSig = bl2adlPDdiff(any((blPDs - (adlCI+adlPDs))'<0)' & any((blPDs - (adlCI+adlPDs))'>0)');

subplot(2,3,5), hold on;
shadeX3 = [0 -180 -180 0];
shadeY3 = [0 0 length(blPDs) length(blPDs)];
fill(shadeX3,shadeY3,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)

histogram(bl2adlPDdiff_sig,-90:9:90,'FaceColor',[0 215 157]/256,'FaceAlpha',1)
% histogram(bl2blPDdiff,-90:9:90,'FaceColor',[239 68 232]/256,'FaceAlpha',0.5)

% histogram(bl2adlPDdiff_notSig,-90:9:90,'FaceColor',[239 68 232]/256,'FaceAlpha',0.8)

line([median(bl2adlPDdiff_sig) median(bl2adlPDdiff_sig)],[0 length(blPDs)],'Color',[0 215 157]/256,'LineWidth',2)
% line([median(bl2blPDdiff) median(bl2blPDdiff)],[0 length(blPDs)],'Color',[239 68 232]/256,'LineWidth',2)

% line([median(bl2adlPDdiff_notSig) median(bl2adlPDdiff_notSig)],[0 length(blPDs)],'Color',[239 68 232]/256,'LineWidth',2)

xlabel('PD difference')
set(gca,'FontName','Helvetica','FontSize',12),axis([-90 90 0 length(blPDs)])
text(-85,10,...
    ['n = ' num2str(length(bl2adlPDdiff_sig)) '/' num2str(numel(popIdx))],...
    'FontName','Helvetica','FontSize',12)

%bl2woe PD scatter
subplot(2,3,3), hold on;
shadeX = rad2deg([0 2*pi pi 0]);
shadeY = rad2deg([0 2*pi 2*pi pi]);
fill(shadeX,shadeY,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
shadeX2 = rad2deg([pi 2*pi 2*pi pi]);
shadeY2 = rad2deg([0 0 pi 0]);
fill(shadeX2,shadeY2,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)

plot(blPDs,woePDs,'.','color',[213 153 0]/256,'MarkerSize',15);
errorbar(blPDs,woePDs,woeCI(:,1),woeCI(:,2),blCI(:,1),blCI(:,2),...
    'LineStyle','none','Color',[213 153 0]/256);
set(gca,'FontName','Helvetica','FontSize',12,'YTick',ticks,'XTick',ticks), axis([0 360 0 360])
% xlabel('Baseline PDs'), ylabel('Washout (early) PDs')
ylabel('Early Washout PDs','color',[213 153 0]/256)


%bl2woe differences
bl2woePDdiff = angleDiff(woePDs,blPDs,false,true);
bl2woePDdiff_sig = bl2woePDdiff(~any((blPDs - (woeCI+woePDs))'<0)');
% bl2woePDdiff_notSig = bl2woePDdiff(any((blPDs - (woeCI+woePDs))'<0)');


subplot(2,3,6), hold on;

shadeX3 = [0 -180 -180 0];
shadeY3 = [0 0 length(blPDs) length(blPDs)];
fill(shadeX3,shadeY3,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)

histogram(bl2woePDdiff_sig,-90:9:90,'FaceColor',[213 153 0]/256,'FaceAlpha',1)
% histogram(bl2blPDdiff_sig,-90:9:90,'FaceColor',[239 68 232]/256,'FaceAlpha',0.5)

% histogram(bl2woePDdiff_notSig,-90:9:90,'FaceColor',[239 68 232]/256,'FaceAlpha',0.8)

line([median(bl2woePDdiff_sig) median(bl2woePDdiff_sig)],[0 length(blPDs)],'Color',[213 153 0]/256,'LineWidth',2)
% line([median(bl2blPDdiff) median(bl2blPDdiff)],[0 length(blPDs)],'Color',[239 68 232]/256,'LineWidth',2)

% line([median(bl2woePDdiff_notSig) median(bl2woePDdiff_notSig)],[0 length(blPDs)],'Color',[239 68 232]/256,'LineWidth',2)

set(gca,'FontName','Helvetica','FontSize',12),axis([-90 90 0 length(blPDs)])

text(-85,10,...
    ['n = ' num2str(length(bl2woePDdiff_sig)) '/' num2str(numel(popIdx))],...
    'FontName','Helvetica','FontSize',12)

% %bl2wol PD scatter
% subplot(2,5,5), hold on;
% shadeX = rad2deg([0 2*pi pi 0]);
% shadeY = rad2deg([0 2*pi 2*pi pi]);
% fill(shadeX,shadeY,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
% shadeX2 = rad2deg([pi 2*pi 2*pi pi]);
% shadeY2 = rad2deg([0 0 pi 0]);
% fill(shadeX2,shadeY2,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
% 
% plot(blPDs,wolPDs,'.','color',[162 166 188]/256,'MarkerSize',15);
% errorbar(blPDs,wolPDs,wolCI(:,1),wolCI(:,2),blCI(:,1),blCI(:,2),...
%     'LineStyle','none','Color',[162 166 188]/256);
% set(gca,'FontName','Helvetica','FontSize',12,'YTick',ticks,'XTick',ticks), axis([0 360 0 360])
% title('Late Washout PDs','color',[162 166 188]/256)


% %bl2wol differences
% bl2wolPDdiff = angleDiff(wolPDs,blPDs,false,true);
% bl2wolPDdiff = bl2wolPDdiff(~any((blPDs - (wolCI+wolPDs))'<0)');
% 
% subplot(2,5,10), hold on;
% shadeX3 = [0 -180 -180 0];
% shadeY3 = [0 0 length(blPDs) length(blPDs)];
% fill(shadeX3,shadeY3,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
% 
% histogram(bl2wolPDdiff,-180:9:180,'FaceColor',[162 166 188]/256)
% line([median(bl2wolPDdiff) median(bl2wolPDdiff)],[0 length(blPDs)],'Color',[0 0 0],'LineWidth',2)
% set(gca,'FontName','Helvetica','FontSize',12),axis([-180 180 0 length(blPDs)])


% hfig.Position = [-1.0738 1.2152 0.3637 0.32];