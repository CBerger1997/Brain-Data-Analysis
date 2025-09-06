%% PREPROCESSING
% Loading data - 
root_dir = uigetdir (pwd, 'select folder'); 

% here we tell the toolbox about a hierarchy in our data set, for example
% the structure has folder for a group, session and then subject on the
% lowest level
raw = nirs.io.loadDirectory((root_dir), {'subject'});

for i = 1:length(raw)
    [~, name, ~] = fileparts(raw(i).description);  % Extract name from full path
    raw(i).demographics('subject') = name;         % Set subject name
end

%%
% we look at the stimuli
raw=nirs.viz.StimUtil(raw);
raw.draw;

%%
% renemane stimuli labels
job= nirs.modules.RenameStims();
job.listOfChanges = {
     'r' 'R'
     'g' 'G'
     'j' 'J'
     'b' 'B'
     'c' 'C'
     };
% changes labels to ensure consistency across participants

% Discard labels that we do not need
job = nirs.modules.DiscardStims (job);
job.listOfStims = {'video', 'I', 'P', 'stim_aux1'};

raw = job.run(raw);

% renemane stimuli labels
job= nirs.modules.RenameStims();
job.listOfChanges = {
     'G' 'Girl_running'
     'J' 'Jump_scare'
     'B' 'Bedroom'
     'C' 'Chainsaw'
     'R' 'Rest'
     'V' 'Video'
     };
% change labels now to their relevant video or condition name

raw = job.run(raw);

%% SPLIT DATA INTO SECTIONS
% Choose which section to analyze
section_choice = 'calibration'; % Change to 'experience' as needed

% Split the data
[split_raw, section_info] = split_nirs_data(raw, section_choice);

% Display information about the split
fprintf('\n=== DATA SPLIT SUMMARY ===\n');
if iscell(section_info)
    % Multiple subjects
    for i = 1:length(section_info)
        fprintf('Subject %d - %s section:\n', i, section_info{i}.section);
        fprintf('  Duration: %.1f seconds\n', section_info{i}.split_duration);
        fprintf('  Events: %d\n', section_info{i}.events_included);
    end
else
    % Single subject
    fprintf('Section: %s\n', section_info.section);
    fprintf('Duration: %.1f seconds\n', section_info.split_duration);
    fprintf('Events: %d\n', section_info.events_included);
end

%% CONTINUE WITH STANDARD PREPROCESSING
% Use split_raw for the rest of the analysis
raw = split_raw;

%%
% change stimuli duration 
raw=nirs.design.change_stimulus_duration(raw, {} ,30);

%%
tbl = nirs.createStimulusTable(raw);

%%
j = nirs.modules.ChangeStimulusInfo();
j.ChangeTable = tbl;
raw = j.run(raw);
%%
% trim baseline (5 sec before first stim onset and 5 sec after the last one)
job = nirs.modules.TrimBaseline(job);
job.preBaseline = 5;
job.postBaseline = 5;

% label short separation channel
job = nirs.modules.LabelShortSeperation(job); 
job.max_distance = 15; 

% convert to OD
job = nirs.modules.OpticalDensity(job);

% resample to 4hz to deal with autocorrelation 
job = nirs.modules.Resample(job);
job.Fs = 4; 

% TDDR fixes motion artifacts
job = nirs.modules.TDDR(job); 

%1 (use) or 0 (don't use)
%job.usePCA = 0;

%convert to HB
job = nirs.modules.BeerLambertLaw(job);

raw = job.run(raw); 

participant_raw = raw; 

save participant_raw;

%% 1st level Analysis

job1 = nirs.modules.GLM();
job1. type = 'AR-IRLS';
basis=nirs.design.basis.BoxCar;
% The choice of hemodynamic function depends on the study design and population
job1.AddShortSepRegressors = true;
SubjStats = job1.run(participant_raw);
save SubjStats;

%% 2nd level

job2=nirs.modules.MixedEffects();
job2.formula='beta ~ -1 + cond + (1|subject)';
% The formula depends on the study design, available Matlab documentation on how to specify the formula - https://www.mathworks.com/help/stats/fitlme.html#btyabbf
GroupStats=job2.run(SubjStats);
%change q and p, increase q < 0.07
GroupStats.draw('beta',[],'q < 0.05');
StatsTable=GroupStats.table;

%%
nirs.viz.nirsviewer(SubjStats);

%% 3rd level contrast analysis will test for differences between conditions

% the difference between the first and second condition must sum up to 0

GroupStats.conditions

gr_vs_r = [0 0 1 0 -1];  % Girl_running - Rest
r_vs_gr = [0 0 -1 0 1];  % Rest - Girl_running
js_vs_r = [0 0 0 1 -1];
r_vs_js = [0 0 0 -1 1];
br_vs_r = [1 0 0 0 -1];
r_vs_br = [-1 0 0 0 1];
cs_vs_r = [0 1 0 0 -1];
r_vs_cs = [0 -1 0 0 1];

ContrastGRvsR = GroupStats.ttest(gr_vs_r);
ContrastRvsGR = GroupStats.ttest(r_vs_gr);
ContrastJSvsR = GroupStats.ttest(js_vs_r);
ContrastRvsJS = GroupStats.ttest(r_vs_js);
ContrastBRvsR = GroupStats.ttest(br_vs_r);
ContrastRvsBR = GroupStats.ttest(r_vs_br);
ContrastCSvsR = GroupStats.ttest(cs_vs_r);
ContrastRvsCS = GroupStats.ttest(r_vs_cs);

ContrastGRvsR.draw('beta', [], 'q < 0.05');
ContrastRvsGR.draw('beta', [], 'q < 0.05');
ContrastJSvsR.draw('beta', [], 'q < 0.05');
ContrastRvsJS.draw('beta', [], 'q < 0.05');
ContrastBRvsR.draw('beta', [], 'q < 0.05');
ContrastRvsBR.draw('beta', [], 'q < 0.05');
ContrastCSvsR.draw('beta', [], 'q < 0.05');
ContrastRvsCS.draw('beta', [], 'q < 0.05');

% Define all contrasts and their names
contrastVectors = {
    [0 0 1 0 -1],  'GR_vs_R';
    [0 0 -1 0 1],  'R_vs_GR';
    [0 0 0 1 -1],  'JS_vs_R';
    [0 0 0 -1 1],  'R_vs_JS';
    [1 0 0 0 -1],  'BR_vs_R';
    [-1 0 0 0 1],  'R_vs_BR';
    [0 1 0 0 -1],  'CS_vs_R';
    [0 -1 0 0 1],  'R_vs_CS';
};

%%
summary = table();  % Table to collect summary stats

for i = 1:size(contrastVectors,1)
    vec = contrastVectors{i,1};
    name = contrastVectors{i,2};
    
    % Run the contrast
    C = GroupStats.ttest(vec);
    
    % Count significant channels
    nSig_q05 = sum(C.q < 0.05);
    nSig_q10 = sum(C.q < 0.10);
    
    % Display
    fprintf('%s: %d significant channels at q<0.05, %d at q<0.10\n', ...
        name, nSig_q05, nSig_q10);
    
    % Save table of significant channels
    rows = C.q < 0.05;
    T = C.table;            % Store the table in a variable
    sigTable = T(rows, :);  % Index into the variable, not into C.table directly

    
    if ~isempty(sigTable)
        writetable(sigTable, [name '_q05.csv']);
    end
    
    % Plot
    C.draw('beta', [], 'q < 0.05');
    title([name ' (q < 0.05)']);
    
    % Add to summary table
    summary = [summary; table({name}, nSig_q05, nSig_q10, ...
        'VariableNames', {'Contrast', 'Sig_q05', 'Sig_q10'})];
end

%%
T = ContrastGRvsR.table;  % The results table

% Find the row for source 8 and detector 2 which is what we are interested in
idx = (T.source == 8) & (T.detector == 2);

disp("Girl Running Vs Rest");

if any(idx)
    channelStats = T(idx, :);
    disp(channelStats(:, {'source', 'detector', 'type', 'beta', 'p', 'q'}));
else
    disp('No data found for source 8 and detector 2.');
end

T = ContrastJSvsR.table;  % Your results table

% Find the row for source 8 and detector 2
idx = (T.source == 8) & (T.detector == 2);

disp("Jump Scare Vs Rest");

if any(idx)
    channelStats = T(idx, :);
    disp(channelStats(:, {'source', 'detector', 'type', 'beta', 'p', 'q'}));
else
    disp('No data found for source 8 and detector 2.');
end

% NIRS Toolbox visualization:
fig = figure;
GroupStats.draw('beta',[], 'q < 0.05');  % This plots significant betas
title('Significant HbO Changes (Beta) - Jump Scare vs Rest');

T = ContrastCSvsR.table;  % Your results table

% Find the row for source 8 and detector 2
idx = (T.source == 8) & (T.detector == 2);

disp("Chainsaw Vs Rest");

if any(idx)
    channelStats = T(idx, :);
    disp(channelStats(:, {'source', 'detector', 'type', 'beta', 'p', 'q'}));
else
    disp('No data found for source 8 and detector 2.');
end

T = ContrastBRvsR.table;  % Your results table

% Find the row for source 8 and detector 2
idx = (T.source == 8) & (T.detector == 2);

disp("Bedroom Vs Rest");

if any(idx)
    channelStats = T(idx, :);
    disp(channelStats(:, {'source', 'detector', 'type', 'beta', 'p', 'q'}));
else
    disp('No data found for source 8 and detector 2.');
end

%%

jsRatings = readtable('Participant Condition Ratings.csv');

% Define mapping from short codes to full names used in GLM
conditionMap = containers.Map(...
    {'G', 'J', 'B', 'C'}, ...
    {'Girl_running', 'Jump_scare', 'Bedroom', 'Chainsaw'});

numSubjects = length(SubjStats);
betas = zeros(numSubjects, 1);
ratings = zeros(numSubjects, 1);

uniqueSubjects = unique(jsRatings.Participant);
for i = 1:length(uniqueSubjects)
    count = sum(jsRatings.Participant == uniqueSubjects(i));
    fprintf('Participant %d has %d ratings\n', uniqueSubjects(i), count);
end

for i = 1:numSubjects
    % Get subject name from fNIRS file
    subjectName = SubjStats(i).demographics('subject');
    
    % Extract numeric part
    subjIDstr = regexp(subjectName, '\d+', 'match');
    subjID = str2double(subjIDstr{1});
    
    % Find matching row in the ratings table
    row = jsRatings.Participant == subjID;
    
    if any(row)
        % Average beta across all channels for this participant
        betas(i) = mean(SubjStats(i).beta);  
        ratings(i) = mean(jsRatings.Rating(row));
    else
        % Mark as NaN if no match found
        betas(i) = NaN;
        ratings(i) = NaN;
        warning('No rating found for subject %s (ID %d)', subjectName, subjID);
    end
end

% Remove NaNs
valid = ~isnan(betas) & ~isnan(ratings);
X = betas(valid);
Y = ratings(valid);

% Correlation
[r,p] = corr(X, Y, 'Type', 'Spearman');

fprintf('Correlation between subjective rating and beta: r = %.3f, p = %.3f\n', r, p);


matchedBetas = X;
matchedRatings = Y;

scatter(matchedRatings, matchedBetas, 'filled');
xlabel('Subjective Ratings');
ylabel('Beta Values (e.g. HbO)');
title('Correlation between Ratings and Brain Activation');
lsline;

%%
ratingsTable = readtable('Participant Condition Ratings.csv');


conditions = {'Bedroom', 'Chainsaw', 'Girl_running', 'Jump_scare', 'Rest'};
chan_source = 2;
chan_detector = 8;

nSubs = length(SubjStats);
nConds = length(conditions);

% Preallocate matrix for HbT betas
betas_HbT = NaN(nSubs, nConds);

for i = 1:nSubs
    varTable = SubjStats(i).variables; % table with source, detector, type, cond, etc
    
    for c = 1:nConds
        cond = conditions{c};
        
        % Find indices for HbO and HbR betas for this cond & channel
        condMask = contains(varTable.cond, cond);
        
        idx_HbO = find(condMask & varTable.source == chan_source & varTable.detector == chan_detector & strcmp(varTable.type, 'hbo'));
        idx_HbR = find(condMask & varTable.source == chan_source & varTable.detector == chan_detector & strcmp(varTable.type, 'hbr'));
        
        if ~isempty(idx_HbO) && ~isempty(idx_HbR)
            betaHbO = mean(SubjStats(i).beta(idx_HbO));
            betaHbR = mean(SubjStats(i).beta(idx_HbR));
            betas_HbT(i,c) = betaHbO + betaHbR;
        else
            betas_HbT(i,c) = NaN;
        end
    end
end

% Correlate with subjective ratings per condition
for c = 1:nConds
    cond = conditions{c};
    ratings_cond = NaN(nSubs,1);
    
    for i = 1:nSubs
        subjName = SubjStats(i).demographics('subject');  % Dictionary key access
        subjNum = str2double(extractAfter(subjName, 'P')); 
        
        % Find rating for this participant and condition
        idx = ratingsTable.Participant == subjNum & strcmp(ratingsTable.Condition, cond);
        if any(idx)
            ratings_cond(i) = ratingsTable.Rating(idx);
        else
            ratings_cond(i) = NaN;
        end
    end
    
    % Remove NaNs for correlation
    valid = ~isnan(betas_HbT(:,c)) & ~isnan(ratings_cond);
    if sum(valid) > 2
        [r, p] = corr(betas_HbT(valid,c), ratings_cond(valid));
        fprintf('%s: r = %.3f, p = %.3f, n = %d\n', cond, r, p, sum(valid));
    else
        fprintf('%s: Not enough valid data for correlation\n', cond);
    end
end



%% Storing SubjStats in a CSV

% create the CSV

filename = 'Per Participant Brain HbO.csv';

fid = fopen(filename, 'a');

fprintf(fid, '%s\n', "Participant, Condition, Brain Beta HbO");

for i  = 1 : 35
    currentStat = SubjStats(1, i);

    newBetaCalmHbO = sprintf('%.4f',currentStat.beta(11));
    %newBetaCalmHbR = sprintf('%.4f',currentStat.beta(12));
    newBetaFearHbO = sprintf('%.4f',currentStat.beta(23));
    %newBetaFearHbR = sprintf('%.4f',currentStat.beta(24));

    currentPartStr = int2str(i);
    combString1 = "P";
    newString = append(combString1, currentPartStr);
    combString2 = ",";
    newString2 = append(newString, combString2);

    fprintf(fid, '%s', newString2);
    fprintf(fid, '%s', combString2);
    fprintf(fid, '%s', "Calm");
    fprintf(fid, '%s', combString2);
    fprintf(fid, '%s\n', newBetaCalmHbO);
    
    fprintf(fid, '%s', newString2);
    fprintf(fid, '%s', combString2);
    fprintf(fid, '%s', "Fear");
    fprintf(fid, '%s', combString2);
    fprintf(fid, '%s\n', newBetaFearHbO);
    %fprintf(fid, '%s', combString2);
    %fprintf(fid, '%s', newBetaCalmHbR);
    %fprintf(fid, '%s', combString2);
    %fprintf(fid, '%s\n', newBetaFearHbR);

end

fclose(fid);

function [split_data, section_info] = split_nirs_data(data, section_choice)
% SPLIT_NIRS_DATA - Split NIRS data into calibration and experience sections
%
% Inputs:
%   data - NIRS data object or array of data objects
%   section_choice - 'calibration' or 'experience'
%
% Outputs:
%   split_data - NIRS data object containing only the selected section
%   section_info - struct with information about the split

    if nargin < 2
        section_choice = 'calibration'; % default
    end
    
    % Handle array of data objects
    if length(data) > 1
        split_data = nirs.core.Data.empty();
        section_info = {};
        for i = 1:length(data)
            fprintf('Processing subject %d/%d...\n', i, length(data));
            [temp_data, temp_info] = split_single_subject(data(i), section_choice);
            split_data(i) = temp_data;
            section_info{i} = temp_info;
        end
        return;
    else
        [split_data, section_info] = split_single_subject(data, section_choice);
    end
end

function [split_data, section_info] = split_single_subject(data, section_choice)
    
    % Get all stimulus events
    stimulus_names = data.stimulus.keys;
    all_events = [];
    
    % Collect all events with their types and timings
    for i = 1:length(stimulus_names)
        name = stimulus_names{i};
        stim = data.stimulus(name);
        
        for j = 1:length(stim.onset)
            event.name = name;
            event.onset = stim.onset(j);
            event.dur = stim.dur(j);
            event.amp = stim.amp(j);
            all_events = [all_events; event];
        end
    end
    
    % Sort events by onset time
    [~, sort_idx] = sort([all_events.onset]);
    all_events = all_events(sort_idx);
    
    % Find the split point based on event counts
    split_point = find_split_point(all_events);
    
    % Determine time boundaries
    if strcmp(section_choice, 'calibration')
        time_start = 0;
        % Add 60 seconds after the last calibration event
        time_end = all_events(split_point).onset + all_events(split_point).dur + 60;
        events_to_keep = all_events(1:split_point);
        fprintf('Extracting CALIBRATION section: 0 to %.1f seconds (includes 60s buffer)\n', time_end);
    else % experience
        time_start = all_events(split_point + 1).onset;
        time_end = data.time(end);
        events_to_keep = all_events(split_point + 1:end);
        fprintf('Extracting EXPERIENCE section: %.1f to %.1f seconds\n', time_start, time_end);
    end
    
    % Make sure we don't exceed the actual data length
    if time_end > data.time(end)
        fprintf('Warning: Requested end time (%.1f s) exceeds data length (%.1f s)\n', time_end, data.time(end));
        time_end = data.time(end);
    end
    
    % Find time indices
    time_indices = data.time >= time_start & data.time <= time_end;
    
    % Create new data object
    split_data = data;
    split_data.data = data.data(time_indices, :);
    split_data.time = data.time(time_indices) - time_start; % Reset time to start at 0
    
    % Rebuild stimulus dictionary with adjusted timing
    split_data.stimulus = Dictionary();
    
    % Group events by name
    event_groups = struct();
    for i = 1:length(events_to_keep)
        name = events_to_keep(i).name;
        if ~isfield(event_groups, name)
            event_groups.(name) = [];
        end
        % Adjust onset time relative to new start
        adjusted_event = events_to_keep(i);
        adjusted_event.onset = adjusted_event.onset - time_start;
        event_groups.(name) = [event_groups.(name), adjusted_event];
    end
    
    % Create stimulus events for each group
    group_names = fieldnames(event_groups);
    for i = 1:length(group_names)
        name = group_names{i};
        events = event_groups.(name);
        
        s = nirs.design.StimulusEvents();
        s.name = name;
        s.onset = [events.onset]';
        s.dur = [events.dur]';
        s.amp = [events.amp]';
        
        split_data.stimulus(name) = s;
    end
    
    % Create section info
    section_info.section = section_choice;
    section_info.original_duration = data.time(end);
    section_info.split_duration = split_data.time(end);
    section_info.time_range = [time_start, time_end];
    section_info.events_included = length(events_to_keep);
    
    % Count events by type
    section_info.event_counts = struct();
    for i = 1:length(group_names)
        name = group_names{i};
        section_info.event_counts.(name) = length(event_groups.(name));
    end
    
    fprintf('Section info:\n');
    fprintf('  Duration: %.1f seconds\n', section_info.split_duration);
    fprintf('  Events included: %d\n', section_info.events_included);
    fprintf('  Event counts:\n');
    for i = 1:length(group_names)
        name = group_names{i};
        fprintf('    %s: %d\n', name, section_info.event_counts.(name));
    end
end

function split_point = find_split_point(all_events)
    % Find where calibration ends based on event counts
    % Calibration: 12 Rest + 3 each of Girl_running, Jump_scare, Bedroom, Chainsaw = 24 total events
    % Experience: 6 Rest + 6 Video = 12 total events
    
    event_counts = struct('Rest', 0, 'Girl_running', 0, 'Jump_scare', 0, 'Bedroom', 0, 'Chainsaw', 0, 'Video', 0);
    
    for i = 1:length(all_events)
        name = all_events(i).name;
        if isfield(event_counts, name)
            event_counts.(name) = event_counts.(name) + 1;
        end
        
        % Check if we've reached the end of calibration
        if event_counts.Rest == 12 && event_counts.Girl_running == 3 && ...
           event_counts.Jump_scare == 3 && event_counts.Bedroom == 3 && event_counts.Chainsaw == 3
            split_point = i;
            fprintf('Calibration section ends at event %d (time: %.1f sec)\n', i, all_events(i).onset);
            return;
        end
    end
    
    % If we didn't find the exact pattern, estimate based on total events
    % Calibration should have 24 events, so split there
    if length(all_events) > 24
        split_point = 24;
        fprintf('Warning: Using estimated split point at event 24\n');
    else
        error('Not enough events found to identify calibration/experience split');
    end
end
