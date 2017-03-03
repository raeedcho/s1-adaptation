function pathSetup
% PATHSETUP sets up path for data analysis with CDS and TrialData code, along with proc library and bumpcurl library

if(ispc)
    homeFolder = 'C:\Users\rhc307\';
else
    homeFolder = '/home/raeed/';
end
addpath(genpath([homeFolder filesep 'Projects' filesep 'limblab' filesep 'TrialData']))
addpath(genpath([homeFolder filesep 'Projects' filesep 'limblab' filesep 'ClassyDataAnalysis']))
addpath([homeFolder filesep 'Projects' filesep 'limblab' filesep 'proc-raeed' filesep 'BumpCurl' filesep 'lib' filesep])
addpath(genpath([homeFolder filesep 'Projects' filesep 'limblab' filesep 'proc-raeed' filesep 'lib']))
% cd([homeFolder 'Projects' filesep 'limblab' filesep 'data-raeed' filesep 'BumpCurl' filesep 'Han' filesep datecode filesep])

clear homeFolder
