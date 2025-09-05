%% PREPROCESSING
% Loading data - 
root_dir = uigetdir (pwd, 'select folder'); 

% here we tell the toolbox about a hierarchy in our data set, for example
% the structure has folder for a group, session and then subject on the
% lowest level
raw = nirs.io.loadDirectory((root_dir), {'subject'});

for i = 1:length(raw)
    [~, name, ~] = fileparts(raw(i).description);  % e.g. "P10"
    raw(i).demographics('subject') = upper(strtrim(name));  % Set as 'P10'
end

%%
% we can look at stimili (for example triggers)
raw=nirs.viz.StimUtil(raw);
raw.draw;

%%
% renemane stimuli labels
job= nirs.modules.RenameStims();
job.listOfChanges = {
     'v' 'V'
     };
% the label "fix" will be changed to "F" and the label "e" will be changed to "E"

% we can discard events that we do not need 
job = nirs.modules.DiscardStims (job);
job.listOfStims = {'R', 'G', 'J', 'B', 'C'};

raw = job.run(raw);

% renemane stimuli labels
job= nirs.modules.RenameStims();
job.listOfChanges = {
     'V' 'Video'
     };
% the label "fix" will be changed to "F" and the label "e" will be changed to "E"

raw = job.run(raw);

%% REMAP VIDEO TO CONDITION

% Load mapping table
map = readtable('condition_mapping.csv');  % Columns: Participant, Condition

% Loop through each participant in raw
for i = 1:length(raw)
    sub = raw(i);
    
    % Get subject ID (make sure it's like 'P10', 'P11', etc.)
    pid = sub.demographics('subject');  % adjust key name if needed
    
    % Extract conditions for this participant
    idx = strcmp(map.Participant, pid);
    conds = map.Condition(idx);
    
    if length(conds) ~= 6
        disp(sub.demographics)
        warning(['Expected 6 mappings for ' pid ', found ' num2str(length(conds))]);
        continue;
    end

    % Get 'Video' stimulus object
    if ~any(strcmp('Video', sub.stimulus.keys))
        warning(['No "Video" stimulus found for ' pid]);
        continue;
    end

    videoStim = sub.stimulus('Video');

    % Sort by onset
    [~, sortIdx] = sort(videoStim.onset);
    
    for k = 1:min(6, length(sortIdx))
        newStim = nirs.design.StimulusEvents;
        newStim.onset = videoStim.onset(sortIdx(k));
        newStim.dur = videoStim.dur(sortIdx(k));
        newStim.amp = videoStim.amp(sortIdx(k));
        
        label = conds{k};

        % If this label already exists, append to it
        if ismember(label, sub.stimulus.keys)
            existingStim = sub.stimulus(label);

            existingStim.onset(end+1) = newStim.onset;
            existingStim.dur(end+1)   = newStim.dur;
            existingStim.amp(end+1)   = newStim.amp;

            sub.stimulus(label) = existingStim;
        else
            sub.stimulus(label) = newStim;
        end
    end

    % Save back updated subject
    raw(i) = sub;
end

job = nirs.modules.DiscardStims (job);
job.listOfStims = {'Video'};

raw = job.run(raw);

%% CHANGE CONDITION DURATIONS

% Scene durations in order (in seconds)
sceneDurations = [249, 242, 298, 469, 189, 142];  % Refugee to Supper

% Go through each subject
for i = 1:length(raw)
    sub = raw(i);
    
    % Collect all C1, C2, C3 stim events together
    allEvents = struct('label', {}, 'onset', {}, 'idx', {});
    
    for label = ["Refugee_LF", "Refugee_HF", "Refugee_A", "Camp_LF", "Camp_HF", "Camp_A", "Plants_LF", "Plants_HF", "Plants_A", "Raid_LF", "Raid_HF", "Raid_A", "Processing_LF", "Processing_HF", "Processing_A", "Supper_LF", "Supper_HF", "Supper_A"]
        label = char(label);
        if ~any(strcmp(sub.stimulus.keys, label))
            continue;
        end
        stim = sub.stimulus(label);
        for j = 1:length(stim.onset)
            allEvents(end+1).label = label;
            allEvents(end).onset = stim.onset(j);
            allEvents(end).idx = j;
        end
    end

    % Sort all events by onset time
    [~, sortIdx] = sort([allEvents.onset]);
    allEvents = allEvents(sortIdx);

    % Assign new durations
    for k = 1:min(6, length(allEvents))
        label = allEvents(k).label;
        idx = allEvents(k).idx;
        stim = sub.stimulus(label);
        stim.dur(idx) = sceneDurations(k);
        sub.stimulus(label) = stim;
    end

    raw(i) = sub;  % Save back
end

%%
% Loop over each subject individually
for i = 1:length(raw)
    
    % Subject-specific data
    r = raw(i);

    % Begin pipeline for this subject
    job = nirs.modules.TrimBaseline();
    job.preBaseline = 5;
    job.postBaseline = 5;
    r = job.run(r);

    job = nirs.modules.LabelShortSeperation();
    job.max_distance = 15;
    r = job.run(r);

    job = nirs.modules.OpticalDensity();
    r = job.run(r);

    job = nirs.modules.Resample();
    job.Fs = 4;
    r = job.run(r);

    job = nirs.modules.BeerLambertLaw();
    r = job.run(r);

    % Optional: save preprocessed per-subject
    preprocessed(i) = r;

end