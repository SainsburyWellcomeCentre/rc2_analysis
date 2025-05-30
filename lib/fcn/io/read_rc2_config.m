function config = read_rc2_config(cfg_fname) %#ok<*STOUT>
% READ_RC2_CONFIG Reads the .cfg file associated with a .bin recording.
% 
%   CONFIG = read_rc2_config(FILENAME)
%  
%   Convert the config file, FILENAME, saved by rc2 into a structure.
%   FILENAME is a full path to a .cfg file.
%
%   The config file is a sequence of key/value pairs. It only performs a
%   very simple parse, so if there are new lines etc. in any of the filenames
%   etc. we have a problem...
%
%   See also: read_rc2_bin, read_bin, data_transform



% open, read, close config file
fid = fopen(cfg_fname, 'r');
str = fread(fid, 'char');
fclose(fid);

% convert to character row array
str = char(str(:)');

% split blocks on four newline characters
str_blocks = strsplit(str, '\n\n\n\n');

% remove last entry, which is empty
str_blocks(end) = [];

% parse first block
[keys, vals] = split_key_value(str_blocks{1});

% convert string values into cell arrays, double arrays, scalar vals
scalar = {'saving.ai_min_voltage', 'saving.ai_max_voltage', ...
    'nidaq.ai.rate', 'nidaq.ao.rate', 'nidaq.ao.idle_offset', ...
    'nidaq.co.init_delay', 'nidaq.co.low_samps', 'nidaq.co.high_samps', 'nidaq.ai.scale', 'nidaq.ai.offset'};
cell_array = {'nidaq.ai.channel_names', 'nidaq.ai.channel_ids', ...
    'nidaq.ao.channel_names', 'nidaq.ao.channel_ids', ...
    'nidaq.co.channel_names', 'nidaq.co.channel_ids', ...
    'nidaq.do.channel_names', 'nidaq.do.channel_ids', ...
    'nidaq.di.channel_names', 'nidaq.di.channel_ids'};

% do the conversion
vals = convert_values(keys, vals, scalar, cell_array); %#ok<NASGU>

% enter keys 
for i = 1 : length(keys)
    cmd = sprintf('config.%s = vals{i};', keys{i});
    eval(cmd);
end

% remove first (nidaq) block, which we just 
str_blocks(1) = [];

% for conversion of protocol fields
scalar = {'prot.start_pos', 'prot.stage_pos', 'prot.switch_pos', ...
    'prot.back_limit', 'prot.forward_limit', 'prot.start_dwell_time', ...
    'prot.handle_acquisition', 'prot.wait_for_reward', 'prot.log_trial', ...
    'prot.follow_previous_protocol', 'prot.reward.randomize', ...
    'prot.reward.min_time', 'prot.reward.max_time', 'prot.reward_duration', ...
    'prot.enable_vis_stim', 'prot.initial_trial', 'prot.initiation_speed'};


% parse remaining protocols
for i = 1 : length(str_blocks)
    
    % parse protocol blocks
    [keys, vals] = split_key_value(str_blocks{i});
    
    % convert some entries from strings to numbers
    vals = convert_values(keys, vals, scalar, {}); %#ok<NASGU>
    
    % store in the config structure
    for j = 1 : length(keys)
        subfield = strsplit(keys{j}, '.');
        cmd = sprintf('config.prot(i).%s = vals{j};', subfield{2});
        eval(cmd);
    end
end





function [keys, vals] = split_key_value(str)
%%[keys, vals] = split_key_value(str)
%   Takes individual blocks and parses into keys and values

% split key/value pairs on newlines
str_lines = strsplit(str, '\n');

% remove empty lines??
keys = cell(1, length(str_lines));
vals = cell(1, length(str_lines));


for i = 1 : length(str_lines)
    
    % split line on the first equals character
    this_key = regexp(str_lines{i}, ' = ', 'split', 'once');
    
    % store keys and values
    keys{i} = this_key{1};
    vals{i} = this_key{2};
end



function vals = convert_values(keys, vals, scalar, cell_array)

% convert scalar and cell string arrays
for i = 1 : length(keys)
    if any(strcmp(keys{i}, scalar))
        if strcmp(vals{i}, '---')
            vals{i} = nan;
        else
            vals{i} = eval(sprintf('[%s]', vals{i}));
        end
    end
    if any(strcmp(keys{i}, cell_array))
        vals{i} = strsplit(vals{i}, ',');
    end
end
