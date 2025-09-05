%% BLOCK AVERAGING ANALYSIS (BAA)
% Continue from your preprocessed data

%% OPTION 1: TRADITIONAL BLOCK AVERAGING
% This does actual trial averaging (preserves temporal dynamics)

% Run the block averaging
pre_stim = 5;  % seconds before stimulus
post_stim = 20; % seconds after stimulus
baa_results = perform_block_averaging(preprocessed, pre_stim, post_stim);

%% VISUALIZATION OF BAA RESULTS
% Plot averaged responses for each condition

% Select a subject to visualize (e.g., subject 1)
subj_idx = 2;
subj_data = baa_results(subj_idx);

% Select a channel to plot (e.g., channel 1)
channel_idx = 15;

% Check data dimensions first
fprintf('Subject %d data dimensions:\n', subj_idx);
fprintf('Time vector length: %d\n', length(subj_data.time));
fprintf('Data dimensions: %s\n', mat2str(size(subj_data.data)));
fprintf('Number of conditions: %d\n', length(subj_data.conditions));
fprintf('Conditions: %s\n', strjoin(subj_data.conditions, ', '));

figure;

% Get scenes and conditions
scenes = {'Refugee', 'Camp', 'Plants', 'Raid', 'Processing', 'Supper'};
conditions = {'LF', 'HF', 'A'};
colors = {'blue', 'red', 'green'};

% Plot each condition
for scene_idx = 1:length(scenes)
    scene = scenes{scene_idx};
    
    subplot(2, 3, scene_idx);
    hold on;
    
    legend_entries = {};
    plot_count = 0;
    
    for cond_idx = 1:length(conditions)
        cond_name = sprintf('%s_%s', scene, conditions{cond_idx});
        
        % Find this condition in the data
        cond_data_idx = find(strcmp(subj_data.conditions, cond_name));
        
        if ~isempty(cond_data_idx)
            % Extract data for this condition and channel
            condition_data = subj_data.data(:, channel_idx, cond_data_idx);
            
            % Debug: check dimensions
            fprintf('Condition %s: time=%d, data=%d\n', cond_name, length(subj_data.time), length(condition_data));
            
            % Ensure vectors are the same length
            min_length = min(length(subj_data.time), length(condition_data));
            time_vec = subj_data.time(1:min_length);
            data_vec = condition_data(1:min_length);
            
            % Plot
            plot(time_vec, data_vec, 'Color', colors{cond_idx}, ...
                 'LineWidth', 2);
            plot_count = plot_count + 1;
            legend_entries{plot_count} = conditions{cond_idx};
        else
            fprintf('Condition %s not found for %s\n', cond_name, scene);
        end
    end
    
    title(sprintf('%s - Channel %d', scene, channel_idx));
    xlabel('Time (s)');
    ylabel('Concentration Change');
    if plot_count > 0
        legend(legend_entries, 'Location', 'best');
    end
    grid on;
end

sgtitle(sprintf('Block Averaged Responses - Subject %s', subj_data.demographics('subject')));

%% STATISTICAL ANALYSIS OF BAA RESULTS

scenes = {'Refugee', 'Camp', 'Plants', 'Raid', 'Processing', 'Supper'};
conds  = {'LF', 'HF', 'A'};

% Step 1: Initialize
comparison_results = struct();
for s = 1:length(scenes)
    for c = 1:length(conds)
        key = [scenes{s} '_' conds{c}];
        comparison_results.(key).values = [];
        comparison_results.(key).subjects = {};
    end
end

% Step 2: Populate from BAA
for subj_idx = 1:length(baa_results)
    subj_data = baa_results(subj_idx);
    subj_id = subj_data.demographics('subject');  % dictionary access
    conds = subj_data.conditions;
    time = subj_data.time;
    data = subj_data.data;  % time x channels x conditions

    for cond_idx = 1:length(conds)
        cond_name = conds{cond_idx};

        % Get data for this condition: time x channels
        cond_data = squeeze(data(:,:,cond_idx));

        % Peak extraction: max HbO (channel-wise average)
        % You can choose to average channels first or extract per ROI
        hbo_data = cond_data(:,15);  % assuming HbO = odd channels (typical)
        mean_hbo = mean(hbo_data, 2);     % average across HbO channels

        peak_val = max(mean_hbo);
          
        auc_val = trapz(subj_data.time, hbo_data);
        
        % Store in comparison_results
        if ~isfield(comparison_results, cond_name)
            comparison_results.(cond_name).values = [];
            comparison_results.(cond_name).subjects = {};
        end

        comparison_results.(cond_name).values(end+1) = peak_val;
        comparison_results.(cond_name).subjects{end+1} = subj_id;
    end
end

%% SCENE-BY-CONDITION COMPARISON

scenes = {'Refugee', 'Camp', 'Plants', 'Raid', 'Processing', 'Supper'};
conds = {'LF', 'HF', 'A'};

for s = 1:length(scenes)
    scene = scenes{s};

    group_vals = [];
    group_labels = {};

    for c = 1:length(conds)
        tag = [scene '_' conds{c}];

        if isfield(comparison_results, tag)
            vals = comparison_results.(tag).values;
            group_vals = [group_vals, vals];
            group_labels = [group_labels, repmat({conds{c}}, 1, length(vals))];
        end
    end

    if numel(unique(group_labels)) == 3
        [p, ~, stats] = anova1(group_vals, group_labels, 'off');
        fprintf('%.4f\n', p);
    else
        fprintf('%s: Missing conditions, skipping.\n', scene);
    end
end


%%

scenes = {'Refugee', 'Camp', 'Plants', 'Raid', 'Processing', 'Supper'};
condA = 'HF';
condB = 'A';

for s = 1:length(scenes)
    scene = scenes{s};

    keyA = [scene '_' condA];
    keyB = [scene '_' condB];

    valsA = comparison_results.(keyA).values;
    valsB = comparison_results.(keyB).values;

    if ~isempty(valsA) && ~isempty(valsB)
        % Run unpaired t-test
        [h, p, ci, stats] = ttest2(valsA, valsB);
        fprintf('%s: %s vs %s → p = %.4f (t=%.2f, df=%d)\n', ...
            scene, condA, condB, p, stats.tstat, stats.df);
        mean_diff = mean(valsA) - mean(valsB);
        pooled_sd = sqrt(((std(valsA)^2) + (std(valsB)^2)) / 2);
        cohen_d = mean_diff / pooled_sd;
        fprintf('Cohen''s d = %.2f\n', cohen_d);
        boxplot([valsA, valsB], ...
        [repmat({'LF'}, 1, length(valsA)), repmat({'HF'}, 1, length(valsB))]);
        title('Refugee Scene - HbO AUC: LF vs HF');
        ylabel('HbO AUC');
    else
        fprintf('%s: Not enough data for %s and/or %s\n', scene, condA, condB);
    end
end

%%

% Step 1: Load subjective ratings CSV
ratings = readtable('Participant Ratings (APFE).csv');

% Rename columns manually
ratings.Properties.VariableNames = {'Participant', 'Scene_Name', 'Condition', ...
    'Involved', 'Immersed', 'Fear'};

% Metrics and signals to analyze
metric_types = {'auc', 'peak'};
signal_types = {'hbo', 'hbr'};

% For each signal and metric
metric_type = 'peak'
signal_type = 'hbo'

% Step 2: Match ratings with baa_results and extract AUC
results = [];

for subj_idx = 1:length(baa_results)
    subj_data = baa_results(subj_idx);
    subj_id = subj_data.demographics('subject');

    for c = 1:length(subj_data.conditions)
        cond_name = subj_data.conditions{c};
        scene = extractBefore(cond_name, '_');

        row_idx = find(strcmp(ratings.Participant, subj_id) & strcmp(ratings.Scene_Name, scene));

        if ~isempty(row_idx)
            
            cond_data = squeeze(subj_data.data(:, :, c));
            
            hbo_data = cond_data(:,15);
            
            if strcmpi(signal_type, 'hbr')
                hbo_data = cond_data(:,16);
            end
            
            if strcmpi(metric_type, 'auc')
                val = trapz(subj_data.time, hbo_data);
            elseif strcmpi(metric_type, 'peak')
                mean_hbo = mean(hbo_data, 2);
                val = max(mean_hbo);
            else
                error('Unknown metric type');
            end
            
            % Use a temporary struct, then append
            temp.subject = subj_id;
            temp.scene = scene;
            temp.condition = cond_name;
            temp.auc = val;
            temp.involved = ratings.Involved(row_idx(1));   % ensure scalar
            temp.immersed = ratings.Immersed(row_idx(1));   % ensure scalar
            temp.fear = ratings.Fear(row_idx(1));           % ensure scalar
            
            results = [results; temp];  % append to results array
        end
    end
end

% Now convert to table
results_tbl = struct2table(results);

% Assumes `results_tbl` contains the following columns:
% - condition (e.g., 'Refugee_LF')
% - auc
% - fear
% - immersed
% - involved

% Step 1: Initialize
sceneconds = unique(results_tbl.condition);
summary = struct();

fprintf('\nCorrelation of %s %s vs Ratings per Scene-Condition:\n', upper(metric_type), upper(signal_type));
fprintf('------------------------------------------------------------\n');
fprintf('%-18s | %-5s | %-7s | %-7s | %-7s | %-7s | %-7s | %-7s\n', ...
    'Condition', 'n', 'Fear_r', 'Fear_p', 'Imm_r', 'Imm_p', 'Invol_r');
fprintf('%s\n', repmat('-',1,90));

% Step 2: Loop through each condition
for i = 1:length(sceneconds)
    sc = sceneconds{i};
    idx = strcmp(results_tbl.condition, sc);
    tbl = results_tbl(idx, :);
    n = height(tbl);

    if n < 5
        fprintf('%-18s | %-5d | %s\n', sc, n, 'skipped (n < 5)');
        continue;
    end 
     
    % Calculate correlations
    [r_fear, p_fear] = corr(tbl.auc, tbl.fear, 'Type', 'Spearman');
    [r_imm,  p_imm]  = corr(tbl.auc, tbl.immersed, 'Type', 'Spearman');
    [r_inv,  p_inv]  = corr(tbl.auc, tbl.involved, 'Type', 'Spearman');


    % Print to console
    fprintf('%-18s | %-5d | %-7.2f | %-7.4f | %-7.2f | %-7.4f | %-7.2f\n', ...
        sc, n, r_fear, p_fear, r_imm, p_imm, r_inv);

    % Store results
    summary(i).condition = sc;
    summary(i).n = n;
    summary(i).r_fear = r_fear;
    summary(i).p_fear = p_fear;
    summary(i).r_imm = r_imm;
    summary(i).p_imm = p_imm;
    summary(i).r_invol = r_inv;
    summary(i).p_invol = p_inv;
end

% Step 3: Export to CSV (optional)
summary_tbl = struct2table(summary);

filename = sprintf('results_%s_%s.csv', lower(metric_type), lower(signal_type));
writetable(summary_tbl, filename);
fprintf('Exported: %s\n', filename);

%% EXPORT RESULTS
% Save BAA results to file
save('BAA_results.mat', 'baa_results', 'comparison_results', 'group_LF', 'group_HF', 'group_A');

% Create summary table
summary_table = table();
for subj_idx = 1:n_subjects
    subj_id = comparison_results(subj_idx).subject;
    
    for scene_idx = 1:n_scenes
        scene = scenes{scene_idx};
        
        % Average across channels
        lf_avg = mean(group_LF(subj_idx, scene_idx, :));
        hf_avg = mean(group_HF(subj_idx, scene_idx, :));
        a_avg = mean(group_A(subj_idx, scene_idx, :));
        
        % Add to table
        new_row = table({subj_id}, {scene}, lf_avg, hf_avg, a_avg, ...
                       'VariableNames', {'Subject', 'Scene', 'LF_Peak', 'HF_Peak', 'A_Peak'});
        summary_table = [summary_table; new_row];
    end
end

writetable(summary_table, 'BAA_summary.csv');

fprintf('\nBlock averaging analysis complete!\n');
fprintf('Results saved to: BAA_results.mat\n');
fprintf('Summary saved to: BAA_summary.csv\n');

%% FUNCTION DEFINITION
function baa_data = perform_block_averaging(data, pre_stim, post_stim)
    baa_data = [];
    
    for subj = 1:length(data)
        subj_data = data(subj);
        
        % Get conditions for this subject
        conditions = subj_data.stimulus.keys;
        
        % Initialize storage for this subject
        subj_baa = struct();
        subj_baa.demographics = subj_data.demographics;
        subj_baa.probe = subj_data.probe;
        subj_baa.conditions = {};
        subj_baa.data = [];
        
        % Time vector for averaged trials
        time_vec = (-pre_stim:1/subj_data.Fs:post_stim)';
        subj_baa.time = time_vec;
        
        fprintf('Processing subject %d: %s\n', subj, subj_data.demographics('subject'));
        
        % For each condition
        for c = 1:length(conditions)
            cond_name = conditions{c};
            stim_events = subj_data.stimulus(cond_name);
            
            if isempty(stim_events.onset)
                continue;
            end
            
            % Extract trials for this condition
            trials = [];
            valid_trials = 0;
            
            for trial = 1:length(stim_events.onset)
                onset_time = stim_events.onset(trial);
                [~, onset_sample] = min(abs(subj_data.time - onset_time));                
                
                pre_samples = round(pre_stim * subj_data.Fs);
                post_samples = round(post_stim * subj_data.Fs);
                                
                start_idx = onset_sample - pre_samples;
                end_idx = onset_sample + post_samples;
                                
                if start_idx < 1 || end_idx > size(subj_data.data, 1)
                    fprintf('    Skipping trial %d (index %d–%d out of bounds)\n', trial, start_idx, end_idx);
                    continue;
                end

                
                if start_idx > 0 && end_idx <= size(subj_data.data, 1)
                    trial_data = subj_data.data(start_idx:end_idx, :);
                    
                    % Store trial
                    if isempty(trials)
                        trials = zeros(length(trial_data), size(subj_data.data, 2), length(stim_events.onset));
                    end
                    trials(:, :, trial) = trial_data;
                    valid_trials = valid_trials + 1;
                end
            end
            
            % Average across trials
            if valid_trials > 0
                avg_trial = mean(trials(:, :, 1:valid_trials), 3);
                
                % Store averaged data
                subj_baa.data = cat(3, subj_baa.data, avg_trial);
                subj_baa.conditions{end+1} = cond_name;
                
                fprintf('  %s: %d trials averaged\n', cond_name, valid_trials);
            end
        end
        
        baa_data = [baa_data; subj_baa];
    end
end