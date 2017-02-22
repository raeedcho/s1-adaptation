function [binned_data,emg,trial_table] = extract_bumpcurl_data(cds)

% % Make CDS file
% cds = commonDataStructure();
% cds.file2cds([input_struct.folder filesep input_struct.fname],input_struct.ranBy,input_struct.array,input_struct.monkey,input_struct.lab,'ignoreJumps',input_struct.task,input_struct.mapfile);
% 
% if(isfield(input_struct,'save_folder'))
%     save(input_struct.save_folder,'cds','-v7.3');
% end


%% Extract trial table
trial_table = cds.trials;

%% Move CDS into Experiment class
% DL stuff
ex = experiment();
% set variables to load from cds
ex.meta.hasLfp=false;
ex.meta.hasKinematics=true;
ex.meta.hasForce=true;
ex.meta.hasUnits=true;
ex.meta.hasTrials=true;
ex.meta.hasAnalog=true;
ex.meta.hasEmg=true;
% add session to experiment
ex.addSession(cds);

ex.emg.processDefault();

%% extract emg
emg = ex.emg.data;

%% Bin experiment data
ex.binConfig.include(1).field='units';
ex.binConfig.include(1).which=find([ex.units.data.ID]==0);
ex.binConfig.include(2).field='kin';
ex.binConfig.include(2).which={};
ex.binConfig.include(3).field='force';
ex.binConfig.include(3).which={};
ex.binConfig.include(4).field='emg';
ex.binConfig.include(4).which={};
ex.binConfig.include(5).field='analog';
ex.binConfig.include(5).which={'KinectSyncPulse'};
% base_ex.binConfig.include(3).which=base_ex.analog(2).data.Properties.VariableNames(2:end);%kinect data
ex.firingRateConfig.cropType='tightCrop';
% ex.firingRateConfig.offset=-0.015;
ex.firingRateConfig.kw = 0.01; %10 ms binning

ex.binData()

binned_data = ex.bin.data;