function [base_cds,adapt_cds,wash_cds] = loadsave_bumpcurl_cds(datecode)
%% Load CDS files
% set up input structs
root_folder = ['C:\Users\rhc307\Projects\limblab\data-raeed\BumpCurl\Han\' datecode filesep];
fname_prefix = ['Han_' datecode '_CObumpcurl'];

lab=6;
ranBy='ranByRaeed';
monkey='monkeyHan';
task='taskCObump';
array='arrayLeftS1Area2';
folder=[root_folder 'preCDS'];
mapfile='mapFileC:\Users\rhc307\Projects\limblab\data-raeed\ForceKin\OutOutReach\Han\mapfile\left_S1\SN 6251-001459.cmp';

base_fname= '_baseline_area2EMG_001';
adapt_fname='_adaptation_area2EMG_003';
wash_fname= '_washout_area2EMG_004';

extra_base_fname= '_baseline_EMGextra_001';
extra_adapt_fname='_adaptation_EMGextra_003';
extra_wash_fname= '_washout_EMGextra_004';

base_cds = commonDataStructure();
base_cds.file2cds([folder filesep fname_prefix base_fname],      ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);
base_cds.file2cds([folder filesep fname_prefix extra_base_fname],ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);

adapt_cds = commonDataStructure();
adapt_cds.file2cds([folder filesep fname_prefix adapt_fname],      ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);
adapt_cds.file2cds([folder filesep fname_prefix extra_adapt_fname],ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);

wash_cds = commonDataStructure();
wash_cds.file2cds([folder filesep fname_prefix wash_fname],      ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);
wash_cds.file2cds([folder filesep fname_prefix extra_wash_fname],ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);

%% Save CDS
save_folder = [root_folder 'CDS\'];

if(~isdir(save_folder))
    mkdir(save_folder);
end

save([save_folder '_CDS.mat'],'base_cds','adapt_cds','wash_cds','-v7.3');


