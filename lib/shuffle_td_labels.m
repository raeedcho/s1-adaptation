function trial_data = shuffle_td_labels(trial_data,shuffle_conds)
% function to shuffle labels of trial_data
% Inputs:
%   trial_data - trial data struct
%   shuffle_conds - cell array of condition labels to shuffle (e.g. shuffle_conds={'learning_block'} will shuffle learning block labels)

    % shuffle conditions if we want to
    if ~isempty(shuffle_conds)
        if ischar(shuffle_conds)
            shuffle_conds = {shuffle_conds};
        end
        if iscell(shuffle_conds)
            for shufflenum = 1:length(shuffle_conds)
                if ~ischar(shuffle_conds{shufflenum})
                    error('shuffle_conds must be a char array or a cell array of char arrays')
                else
                    if isfield(trial_data,shuffle_conds{shufflenum})
                        current_labels = {trial_data.(shuffle_conds{shufflenum})};
                        new_labels = current_labels(randperm(length(current_labels)));
                        [trial_data.(shuffle_conds{shufflenum})] = deal(new_labels{:});
                    end
                end
            end
        end
    end
    