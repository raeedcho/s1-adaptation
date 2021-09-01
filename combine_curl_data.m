function combine_curl_data
% script to combine separate BL,AD,WO files from Kyle into one trial_data file and save

    %%
    if ispc
        homefolder = 'C:\Users\Raeed';
    else
        homefolder = '/home/raeed';
    end
    
    dataroot = fullfile(homefolder,'data','project-data','limblab','s1-adapt','td-library','split');
    
    %%
    sessions = {...
        'Chips_20171006',...
        'Duncan_20190911',...
        'Duncan_20191120',...
        'Han_20200113',...
        'Han_20200708'};
    
    epochs = {'bl','ad','wo'};
    
    for sessionnum = 1:length(sessions)
        trial_data = cell(1,length(epochs));
        for epochnum = 1:length(epochs)
            file_info = dir(fullfile(dataroot,sprintf('%s_CObumpcurl_%s*.mat',sessions{sessionnum},epochs{epochnum})));

            temp = load(fullfile(dataroot,file_info.name),'trial_data');
            trial_data{epochnum} = temp.trial_data;
            trial_data{epochnum}.epoch = upper(epochs{epochnum});
            trial_data{epochnum} = reorderTDfields(trial_data{epochnum});
        end
        trial_data = cat(2,trial_data{:});
        save(fullfile(dataroot,sprintf('%s_CObumpcurl_TD.mat',sessions{sessionnum})),'trial_data')
    end
    