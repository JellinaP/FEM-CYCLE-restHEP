% =========================================================================
% ASR PARAMETER TUNING & QUALITY CONTROL
% =========================================================================
% Interactive tool for optimizing ASR (Artifact Subspace Reconstruction)
% parameters before applying to full dataset.
%
% USAGE:
%   1. Run this script after RS1_to_RS4_combined_preASR.m
%   2. Select a representative file (or random sample)
%   3. Test different ASR parameters
%   4. Visualize before/after comparisons
%   5. Save optimal parameters when satisfied
%
% WORKFLOW:
%   - Loads _preASR.set files
%   - Tests ASR with various BurstCriterion and ChannelCriterion values
%   - Shows data quality metrics and visualizations
%   - Exports recommended settings for batch processing
% =========================================================================

clear; clc;

%% =========================== PARAMETERS =================================

% --- Paths ---
DATA_PATH = '/data/pt_02938/rsHEP_analysis/0_data/';
EEGLAB_PATH = '/data/u_prinsen_software/eeglab2024.2.1';
OUTPUT_PATH = '/data/pt_02938/rsHEP_analysis/1_scripts/MINT_RS_EEG_scripts/';

% --- ASR Parameter Ranges to Test ---
BURST_CRITERIA = [5, 10, 20, 40, 100];  % Standard deviation thresholds
CHANNEL_CRITERIA = [0.7, 0.8, 0.85, 0.9];  % Correlation thresholds

% --- Other ASR Settings (typically kept constant) ---
FLATLINE_CRITERION = 5;     % seconds (remove flatlined channels)
LINE_NOISE_CRITERION = 4;   % z-score (remove line-noise contaminated channels)
WINDOW_CRITERION = 'off';   % 'off' = don't remove time windows (we'll do separately)

% --- Channel Configuration ---
VEOG_CHANNEL_IDX = 21;  % vEOG channel should be excluded from ASR

% --- Sampling for Quick QC ---
N_SAMPLES_TO_TEST = 3;  % Number of random files to test initially
RANDOM_SEED = 42;       % For reproducibility

% --- Visualization Settings ---
PLOT_TIME_RANGE = [0 30];  % seconds to plot for comparison
MAX_FREQ_PLOT = 70;        % Hz for PSD plots

%% =========================== INITIALIZE =================================

addpath(EEGLAB_PATH);
addpath(genpath(fullfile(EEGLAB_PATH, 'plugins', 'clean_rawdata')));
eeglab('nogui');

% Find all pre-ASR files
files = dir(fullfile(DATA_PATH, '**', '*_preASR.set'));

if isempty(files)
    error('No *_preASR.set files found! Run RS1_to_RS4_combined_preASR.m first.');
end

fprintf('\n========================================\n');
fprintf('ASR PARAMETER TUNING - QUALITY CONTROL\n');
fprintf('========================================\n');
fprintf('Found %d pre-ASR files\n\n', length(files));

%% =========================== FILE SELECTION =============================

% User choice: test specific file or random sample
fprintf('Select option:\n');
fprintf('  1. Test specific file (you choose)\n');
fprintf('  2. Test random sample of %d files\n', N_SAMPLES_TO_TEST);
fprintf('  3. Test first file only\n');

choice = input('Enter choice (1-3): ');

switch choice
    case 1
        % List all files
        fprintf('\nAvailable files:\n');
        for i = 1:length(files)
            fprintf('  %d. %s\n', i, files(i).name);
        end
        file_idx = input(sprintf('Select file (1-%d): ', length(files)));
        files_to_test = files(file_idx);

    case 2
        % Random sample
        rng(RANDOM_SEED);
        n_to_sample = min(N_SAMPLES_TO_TEST, length(files));
        rand_idx = randperm(length(files), n_to_sample);
        files_to_test = files(rand_idx);
        fprintf('\nRandomly selected %d files for testing\n', n_to_sample);

    case 3
        % First file only
        files_to_test = files(1);
        fprintf('\nTesting first file: %s\n', files(1).name);

    otherwise
        error('Invalid choice');
end

%% =========================== ASR TESTING ================================

results = struct();
result_idx = 0;

for file_i = 1:length(files_to_test)

    fprintf('\n========================================\n');
    fprintf('Testing: %s\n', files_to_test(file_i).name);
    fprintf('========================================\n');

    % Load pre-ASR data
    EEG_orig = pop_loadset(files_to_test(file_i).name, files_to_test(file_i).folder);
    EEG_orig = eeg_checkset(EEG_orig);

    % Separate vEOG channel before ASR
    veogIdx = find(strcmpi({EEG_orig.chanlocs.type}, 'EOG'));
    if ~isempty(veogIdx)
        vEOG_data = EEG_orig.data(veogIdx, :);
        vEOG_chanloc = EEG_orig.chanlocs(veogIdx);
        EEG_for_ASR = pop_select(EEG_orig, 'nochannel', veogIdx);
    else
        vEOG_data = [];
        EEG_for_ASR = EEG_orig;
    end

    % Store original for later comparison
    original_chanlocs = EEG_for_ASR.chanlocs;

    % Test each parameter combination
    fprintf('\nTesting %d parameter combinations...\n', ...
        length(BURST_CRITERIA) * length(CHANNEL_CRITERIA));

    for burst_i = 1:length(BURST_CRITERIA)
        for chan_i = 1:length(CHANNEL_CRITERIA)

            burst_crit = BURST_CRITERIA(burst_i);
            chan_crit = CHANNEL_CRITERIA(chan_i);

            fprintf('  Testing: Burst=%d, Channel=%.2f... ', burst_crit, chan_crit);

            try
                % Apply ASR with current parameters
                EEG_clean = EEG_for_ASR;
                [EEG_clean, HP, BUR, removed_channels] = clean_artifacts(EEG_clean, ...
                    'FlatlineCriterion', FLATLINE_CRITERION, ...
                    'ChannelCriterion', chan_crit, ...
                    'LineNoiseCriterion', LINE_NOISE_CRITERION, ...
                    'BurstCriterion', burst_crit, ...
                    'WindowCriterion', WINDOW_CRITERION, ...
                    'Highpass', 'off');  % Already filtered

                % Calculate quality metrics
                result_idx = result_idx + 1;

                results(result_idx).filename = files_to_test(file_i).name;
                results(result_idx).burst_criterion = burst_crit;
                results(result_idx).channel_criterion = chan_crit;

                % Channels removed
                results(result_idx).n_channels_removed = sum(removed_channels);
                results(result_idx).removed_channel_labels = ...
                    {original_chanlocs(removed_channels).labels};
                results(result_idx).pct_channels_removed = ...
                    100 * sum(removed_channels) / length(removed_channels);

                % Data modified by ASR
                data_diff = EEG_for_ASR.data - EEG_clean.data;
                results(result_idx).pct_data_modified = ...
                    100 * sum(data_diff(:) ~= 0) / numel(data_diff);

                % RMS of modification
                results(result_idx).rms_change = rms(data_diff(:));

                % Variance retained
                var_original = var(EEG_for_ASR.data(:));
                var_clean = var(EEG_clean.data(:));
                results(result_idx).pct_variance_retained = ...
                    100 * var_clean / var_original;

                % Power in high frequencies (muscle artifact indicator)
                [psd_orig, freqs] = pwelch(EEG_for_ASR.data', [], [], [], EEG_for_ASR.srate);
                [psd_clean, ~] = pwelch(EEG_clean.data', [], [], [], EEG_clean.srate);

                muscle_band = freqs >= 20 & freqs <= 70;
                results(result_idx).muscle_power_reduction_db = ...
                    10 * log10(mean(psd_clean(muscle_band,:), 'all') / ...
                               mean(psd_orig(muscle_band,:), 'all'));

                fprintf('Done. Removed %d ch, %.1f%% data modified\n', ...
                    results(result_idx).n_channels_removed, ...
                    results(result_idx).pct_data_modified);

            catch ME
                fprintf('FAILED: %s\n', ME.message);
                results(result_idx).error = ME.message;
            end

        end
    end
end

%% =========================== RESULTS SUMMARY ============================

fprintf('\n========================================\n');
fprintf('RESULTS SUMMARY\n');
fprintf('========================================\n');

% Convert results to table
results_table = struct2table(results);

% Remove error column if no errors
if isfield(results, 'error')
    results_table.error = [];
end

% Sort by burst criterion
results_table = sortrows(results_table, 'burst_criterion');

% Display key metrics
fprintf('\nKey Metrics by Parameter Setting:\n');
fprintf('%-15s %-10s | %s | %s | %s | %s\n', ...
    'Burst Crit', 'Chan Crit', 'Ch Removed', '% Data Mod', '% Var Kept', 'Muscle -dB');
fprintf('%s\n', repmat('-', 1, 90));

for i = 1:height(results_table)
    fprintf('%-15d %-10.2f | %10d | %10.1f | %10.1f | %10.1f\n', ...
        results_table.burst_criterion(i), ...
        results_table.channel_criterion(i), ...
        results_table.n_channels_removed(i), ...
        results_table.pct_data_modified(i), ...
        results_table.pct_variance_retained(i), ...
        results_table.muscle_power_reduction_db(i));
end

%% =========================== VISUALIZATION ==============================

% Create comparison figure
fprintf('\nGenerating comparison plots...\n');

% Select a few representative parameter combinations to visualize
viz_params = [
    5, 0.85;    % Aggressive
    20, 0.85;   % Moderate
    100, 0.85;  % Lenient
];

figure('Position', [100 100 1400 800], 'Name', 'ASR Parameter Comparison');

for viz_i = 1:size(viz_params, 1)

    burst_val = viz_params(viz_i, 1);
    chan_val = viz_params(viz_i, 2);

    % Find matching result
    match_idx = find(results_table.burst_criterion == burst_val & ...
                     results_table.channel_criterion == chan_val, 1);

    if isempty(match_idx)
        continue;
    end

    % Reload and reprocess for visualization
    EEG_orig = pop_loadset(files_to_test(1).name, files_to_test(1).folder);

    % Remove vEOG
    veogIdx = find(strcmpi({EEG_orig.chanlocs.type}, 'EOG'));
    if ~isempty(veogIdx)
        EEG_for_viz = pop_select(EEG_orig, 'nochannel', veogIdx);
    else
        EEG_for_viz = EEG_orig;
    end

    % Apply ASR
    EEG_clean = clean_artifacts(EEG_for_viz, ...
        'FlatlineCriterion', FLATLINE_CRITERION, ...
        'ChannelCriterion', chan_val, ...
        'LineNoiseCriterion', LINE_NOISE_CRITERION, ...
        'BurstCriterion', burst_val, ...
        'WindowCriterion', WINDOW_CRITERION, ...
        'Highpass', 'off');

    % Plot time series
    subplot(3, 2, (viz_i-1)*2 + 1);
    time_range = PLOT_TIME_RANGE(1)*EEG_for_viz.srate : PLOT_TIME_RANGE(2)*EEG_for_viz.srate;
    time_range = time_range(time_range <= size(EEG_for_viz.data, 2));

    plot_channels = [1:min(10, EEG_for_viz.nbchan)];  % Plot first 10 channels
    time_vec = time_range / EEG_for_viz.srate;

    plot(time_vec, EEG_for_viz.data(plot_channels, time_range)' + ...
         (1:length(plot_channels))*100, 'k', 'LineWidth', 0.5);
    hold on;
    plot(time_vec, EEG_clean.data(plot_channels, time_range)' + ...
         (1:length(plot_channels))*100, 'b', 'LineWidth', 0.5);

    title(sprintf('Burst=%d, Chan=%.2f (%.1f%% modified)', ...
        burst_val, chan_val, results_table.pct_data_modified(match_idx)));
    xlabel('Time (s)');
    ylabel('Channel + offset');
    legend({'Original', 'ASR Cleaned'}, 'Location', 'best');
    grid on;

    % Plot PSD
    subplot(3, 2, (viz_i-1)*2 + 2);
    [psd_orig, freqs] = pwelch(EEG_for_viz.data', [], [], [], EEG_for_viz.srate);
    [psd_clean, ~] = pwelch(EEG_clean.data', [], [], [], EEG_clean.srate);

    freq_mask = freqs <= MAX_FREQ_PLOT;

    semilogy(freqs(freq_mask), mean(psd_orig(freq_mask,:), 2), 'k', 'LineWidth', 1.5);
    hold on;
    semilogy(freqs(freq_mask), mean(psd_clean(freq_mask,:), 2), 'b', 'LineWidth', 1.5);

    title(sprintf('PSD: %.1f dB muscle reduction', ...
        results_table.muscle_power_reduction_db(match_idx)));
    xlabel('Frequency (Hz)');
    ylabel('Power (µV²/Hz)');
    legend({'Original', 'ASR Cleaned'}, 'Location', 'best');
    grid on;
    xlim([0 MAX_FREQ_PLOT]);
end

sgtitle(sprintf('ASR Parameter Comparison - %s', files_to_test(1).name), ...
    'Interpreter', 'none');

% Save figure
saveas(gcf, fullfile(OUTPUT_PATH, 'ASR_parameter_comparison.png'));
fprintf('Saved comparison figure to: %s\n', ...
    fullfile(OUTPUT_PATH, 'ASR_parameter_comparison.png'));

%% =========================== RECOMMENDATIONS ============================

fprintf('\n========================================\n');
fprintf('RECOMMENDATIONS\n');
fprintf('========================================\n');

% Find balanced parameters (moderate cleaning)
% Good heuristics:
% - Remove 0-3 channels on average
% - Modify 5-20% of data points
% - Retain >95% of variance
% - Reduce muscle power by 2-5 dB

good_params = results_table( ...
    results_table.n_channels_removed <= 3 & ...
    results_table.pct_data_modified >= 5 & ...
    results_table.pct_data_modified <= 20 & ...
    results_table.pct_variance_retained >= 95, :);

if ~isempty(good_params)
    fprintf('\nRecommended parameter combinations:\n\n');
    disp(good_params(:, {'burst_criterion', 'channel_criterion', ...
        'n_channels_removed', 'pct_data_modified', 'pct_variance_retained'}));

    % Suggest most balanced
    [~, best_idx] = min(abs(good_params.pct_data_modified - 10));
    fprintf('\n** BEST BALANCED CHOICE **\n');
    fprintf('  Burst Criterion: %d\n', good_params.burst_criterion(best_idx));
    fprintf('  Channel Criterion: %.2f\n', good_params.channel_criterion(best_idx));
else
    fprintf('\nNo parameters met the balanced criteria.\n');
    fprintf('Consider these moderate options:\n');

    moderate = results_table(results_table.burst_criterion == 20 & ...
                            results_table.channel_criterion == 0.85, :);
    if ~isempty(moderate)
        disp(moderate(:, {'burst_criterion', 'channel_criterion', ...
            'n_channels_removed', 'pct_data_modified', 'pct_variance_retained'}));
    end
end

fprintf('\nGuidelines:\n');
fprintf('  - Burst Criterion: Lower = more aggressive (typical: 5-20)\n');
fprintf('  - Channel Criterion: Lower = remove more channels (typical: 0.7-0.85)\n');
fprintf('  - For resting-state HEP: aim for moderate cleaning (Burst ~20)\n');
fprintf('  - Excessive cleaning may remove true HEP signals!\n');

%% =========================== SAVE PARAMETERS ============================

fprintf('\n========================================\n');
save_params = input('Save optimal parameters to file? (y/n): ', 's');

if strcmpi(save_params, 'y')
    fprintf('\nEnter optimal parameters:\n');
    optimal_burst = input('  Burst Criterion: ');
    optimal_channel = input('  Channel Criterion: ');

    % Create parameter file
    param_struct = struct();
    param_struct.burst_criterion = optimal_burst;
    param_struct.channel_criterion = optimal_channel;
    param_struct.flatline_criterion = FLATLINE_CRITERION;
    param_struct.line_noise_criterion = LINE_NOISE_CRITERION;
    param_struct.window_criterion = WINDOW_CRITERION;
    param_struct.date_created = datestr(now);
    param_struct.tested_files = {files_to_test.name}';

    % Save as .mat file
    save(fullfile(OUTPUT_PATH, 'optimal_ASR_parameters.mat'), 'param_struct');
    fprintf('\n✓ Parameters saved to: optimal_ASR_parameters.mat\n');

    % Also save as readable text
    fid = fopen(fullfile(OUTPUT_PATH, 'optimal_ASR_parameters.txt'), 'w');
    fprintf(fid, 'OPTIMAL ASR PARAMETERS\n');
    fprintf(fid, '======================\n');
    fprintf(fid, 'Date: %s\n\n', datestr(now));
    fprintf(fid, 'BurstCriterion: %d\n', optimal_burst);
    fprintf(fid, 'ChannelCriterion: %.2f\n', optimal_channel);
    fprintf(fid, 'FlatlineCriterion: %d\n', FLATLINE_CRITERION);
    fprintf(fid, 'LineNoiseCriterion: %d\n', LINE_NOISE_CRITERION);
    fprintf(fid, 'WindowCriterion: %s\n', WINDOW_CRITERION);
    fprintf(fid, '\nTested on files:\n');
    for i = 1:length(files_to_test)
        fprintf(fid, '  - %s\n', files_to_test(i).name);
    end
    fclose(fid);

    fprintf('✓ Parameters saved to: optimal_ASR_parameters.txt\n');
end

fprintf('\n========================================\n');
fprintf('QC Complete!\n');
fprintf('Next: Use optimal parameters in RS5_ASR_and_ICA.m\n');
fprintf('========================================\n\n');
