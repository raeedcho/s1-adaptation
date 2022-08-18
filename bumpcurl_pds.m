function bumpcurl_pds
% script to check distance in population space from final learning trials

%% Set up meta info
    if ispc
        homefolder = 'C:\Users\Raeed';
    else
        homefolder = '/home/raeed';
    end
    
    dataroot = fullfile(homefolder,'data','project-data','limblab','s1-adapt','td-library');
    file_info = dir(fullfile(dataroot,'*CO*.mat'));
    filenames = horzcat({file_info.name})';

    % Load up data into cell array
    trial_data_cell = load_curl_data(fullfile(dataroot,filenames));

%% calculate PDs

    active_params = struct(...
        'save_figures', true,...
        'figsavedir', fullfile(homefolder,'Wiki','0-projects','s1-adaptation','figures','pds','active'),...
        'alignment_event', 'idx_movement_on',...
        'arrays_to_plot', {{'M1','S1'}},...
        'num_boots', 10 ...
    );
    
    pd_table = iterate_pds(trial_data_cell,active_params);

    % plot out all mahal curves, colored by array
%     plot_pd_summary(pd_table,active_params)

end

%%%% Subfunctions %%%%
function pd_table = iterate_pds(trial_data_cell,params)
% calculates a PD table for each block of learning (including fake blocks
% assigned in BL and WO). Returns a table where each row is the PD information
% of a single neuron of an array during one learning block.

    alignment_event = 'idx_movement_on';
    arrays_to_plot = {'M1','S1'};
    num_boots = 10;
    trim_start = 0;
    trim_end = 0.25;
    assignParams(who,params)

%     if strcmpi(alignment_event,'idx_movement_on')
%         trim_start = -0.1;
%         trim_end = 0.5;
%     elseif strcmpi(alignment_event,'idx_bumpTime')
%         trim_start = 0;
%         trim_end = 0.3;
%     end

    pd_file=cell(length(trial_data_cell),1);
    for filenum = 1:length(trial_data_cell)
        td = trial_data_cell{filenum};

        if ~isfield(td,alignment_event)
            continue
        end

        valid_trials = ~isnan(cat(1,td.(alignment_event)));
        if ~any(valid_trials)
            continue
        end
        td = trimTD(td(valid_trials),struct(...
            'idx_start',{{alignment_event,trim_start/td(1).bin_size}},...
            'idx_end',{{alignment_event,trim_end/td(1).bin_size}},...
            'remove_short',true));

        spikes_in_td = getTDfields(td,'spikes');
        pd_arr = cell(1,length(spikes_in_td));
        for arraynum = 1:length(spikes_in_td)
            % check to only plot the arrays we want
            if ~contains(spikes_in_td{arraynum},arrays_to_plot)
                continue
            end

            arrayname = strrep(spikes_in_td{arraynum},'_spikes','');
            unit_guide = td(1).(strcat(arrayname,'_unit_guide'));
            unit_guide_char = strsplit(strip(sprintf('ch%du%d ',unit_guide')),' ')';

            % for each learning block, calculate PDs
            learning_blocks = sort(unique([td.learning_block]));
            pd_block = cell(1,length(learning_blocks));
            for blocknum = 1:length(learning_blocks)
                [~,td_block] = getTDidx(td,'learning_block',learning_blocks(blocknum));
                pd_params = struct(...
                    'out_signals',strcat(arrayname,'_spikes'),...
                    'out_signal_names',unit_guide,...
                    'in_signals','vel',...
                    'num_boots',num_boots,...
                    'verbose',false,...
                    'meta',struct('array',arrayname,'epoch',td_block(1).epoch,'learning_block',learning_blocks(blocknum))...
                );
                pd_block{blocknum} = getTDPDs(td_block,pd_params);
            end
            pd_arr{arraynum} = vertcat(pd_block{:});
        end
        pd_file{filenum} = vertcat(pd_arr{:});
    end
    pd_table = vertcat(pd_file{:});
end

function plot_pd_summary(pd_table,params)

    % replace non-tuned rows with NaNs
    pd_table_nans = pd_table;
    pd_table_nans{~pd_table_nans.velTuned,{'velPD','velPDCI'}} = NaN;

    % get ranges for common sets of learning blocks across sessions
    learning_block_ranges = groupsummary(pd_table_nans,{'monkey','date'},{'min','max'},'learning_block');
    min_common_learning_block = max(learning_block_ranges.min_learning_block);
    max_common_learning_block = min(learning_block_ranges.max_learning_block);

    % remove learning blocks outside of common ranges
    pd_table_nans(pd_table_nans.learning_block<min_common_learning_block,:) = [];
    pd_table_nans(pd_table_nans.learning_block>max_common_learning_block,:) = [];

    % pivot table to get learning blocks
    pd_table_nans.signalID = strsplit(strip(sprintf('ch%du%d ',pd_table_nans{:,'signalID'}')),' ')';
    pd_curves = unstack(pd_table_nans,'velPD','learning_block','groupingvariable',{'date','signalID'},'constantvariables',{'monkey','task','array'},'aggregationfunction',@circ_mean);
end