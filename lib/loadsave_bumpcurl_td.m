function trial_data = loadsave_bumpcurl_td(datecode,do_save)

%% Take either 3 cds files or datecode
if(~ischar(datecode))
    error('Input type error: expected datecode string')
end

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

%% Extract trial data from cds files
params.array_alias = {'LeftS1Area2','S1'};
params.exclude_units = 255;
params.event_list = {'bumpTime';'bumpDir'};
params.trial_results = {'R','A','F','I'};
meta = struct('force_direction',1.48,'epoch','BL');
params.meta = meta;
trial_data_base = parseFileByTrial(base_cds,params);

params.meta.epoch = 'AD';
trial_data_adapt = parseFileByTrial(adapt_cds,params);

params.meta.epoch = 'WO';
trial_data_wash = parseFileByTrial(wash_cds,params);

trial_data = cat(2,trial_data_base,trial_data_adapt,trial_data_wash);

%% Save files
if(do_save)
    save_folder = [root_folder 'CDS\'];

    if(~isdir(save_folder))
        mkdir(save_folder);
    end

    save([save_folder '_CDS.mat'],'base_cds','adapt_cds','wash_cds','-v7.3');

    save_folder = [root_folder 'TD\'];

    if(~isdir(save_folder))
        mkdir(save_folder);
    end

    save([save_folder '_TD.mat'],'base_cds','adapt_cds','wash_cds','-v7.3');
end

