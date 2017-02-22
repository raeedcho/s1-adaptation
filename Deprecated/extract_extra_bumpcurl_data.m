function [data,emg] = extract_extra_bumpcurl_data(cds)

% Make CDS file
% cds = commonDataStructure();
% cds.file2cds([input_struct.folder filesep input_struct.fname],input_struct.ranBy,input_struct.array,input_struct.monkey,input_struct.lab,'ignoreJumps',input_struct.task,input_struct.mapfile);
% 
% if(isfield(input_struct,'save_folder'))
%     save(input_struct.save_folder,'cds','-v7.3');
% end

%% Move CDS into Experiment class
ex = experiment();
% set variables to load from cds
ex.meta.hasLfp=false;
ex.meta.hasKinematics=false;
ex.meta.hasForce=false;
ex.meta.hasUnits=false;
ex.meta.hasTrials=false;
ex.meta.hasAnalog=true;
ex.meta.hasEmg=true;
% add session to experiment
ex.addSession(cds);

ex.emg.processDefault();

%% extract emg
emg = ex.emg.data;

%% Bin experiment data
ex.binConfig.include(1).field='emg';
ex.binConfig.include(1).which={};
ex.binConfig.include(2).field='analog';
ex.binConfig.include(2).which={'Sync'};
ex.firingRateConfig.cropType='tightCrop';
% ex.firingRateConfig.offset=-0.015;
ex.firingRateConfig.kw = 0.01; % 10 ms binning

ex.binData()

data = ex.bin.data;