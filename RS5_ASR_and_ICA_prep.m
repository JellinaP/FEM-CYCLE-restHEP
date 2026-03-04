% =========================================================================
% RS5: ASR CLEANING & ICA PREPARATION
% =========================================================================
% Applies optimized ASR parameters (from ASR_parameter_tuning_QC.m) to
% pre-processed data, then prepares for ICA/AMICA.
%
% INPUT:  Files ending in '_preASR.set'
% OUTPUT: Files ending in '_clean.set' (ready for ICA)
%
% STEPS:
%   1. Load pre-ASR data
%   2. Apply ASR (Artifact Subspace Reconstruction)
%   3. Interpolate removed channels back to full montage
%   4. Re-reference to average
%   5. Add vEOG and ECG channels back
%   6. Save cleaned data ready for ICA
% =========================================================================

clear; clc;

%% =========================== PARAMETERS =================================

% --- Paths ---
DATA_PATH = '/data/pt_02938/rsHEP_analysis/0_data/';
EEGLAB_PATH = '/data/u_prinsen_software/eeglab2024.2.1';
PARAM_FILE = '/data/pt_02938/rsHEP_analysis/1_scripts/MINT_RS_EEG_scripts/optimal_ASR_parameters.mat';

% --- ASR Parameters (will be overridden if parameter file exists) ---
BURST_CRITERION = 20;           % SD threshold for burst rejection (lower = more aggressive)
CHANNEL_CRITERION = 0.85;       % Correlation threshold for bad channels (lower = more aggressive)
FLATLINE_CRITERION = 5;         % Seconds (remove flatlined channels)
LINE_NOISE_CRITERION = 4;       % Z-score (remove line-noise channels)
WINDOW_CRITERION = 'off';       % 'off' = don't remove time windows automatically

% --- Channel Configuration ---
VEOG_CHANNEL_IDX = 21;  % vEOG at TP10 (channel index before ECG removal)

% --- Processing Options ---
FORCE_REPROCESS = false;        % Set to true to reprocess already-processed files
GENERATE_QC_PLOTS = true;       % Generate before/after PSD plots
PSD_MAX_FREQ = 70;              % Maximum frequency for PSD plots (Hz)

%% =========================== INITIALIZE =================================

addpath(EEGLAB_PATH);
addpath(genpath(fullfile(EEGLAB_PATH, 'plugins', 'clean_rawdata')));
eeglab('nogui');

% Load optimal parameters if available
if exist(PARAM_FILE, 'file')
    fprintf('Loading optimal ASR parameters from file...\n');
    param_data = load(PARAM_FILE);
    params = param_data.param_struct;

    BURST_CRITERION = params.burst_criterion;
    CHANNEL_CRITERION = params.channel_criterion;
    FLATLINE_CRITERION = params.flatline_criterion;
    LINE_NOISE_CRITERION = params.line_noise_criterion;
    WINDOW_CRITERION = params.window_criterion;

    fprintf('  Burst Criterion: %d\n', BURST_CRITERION);
    fprintf('  Channel Criterion: %.2f\n', CHANNEL_CRITERION);
else
    fprintf('No optimal parameter file found. Using default parameters.\n');
    fprintf('  Burst Criterion: %d\n', BURST_CRITERION);
    fprintf('  Channel Criterion: %.2f\n', CHANNEL_CRITERION);
    fprintf('  (Run ASR_parameter_tuning_QC.m first to optimize!)\n');
end

% Find all pre-ASR files
files = dir(fullfile(DATA_PATH, '**', '*_preASR.set'));

if isempty(files)
    error('No *_preASR.set files found! Run RS1_to_RS4_combined_preASR.m first.');
end

fprintf('\n========================================\n');
fprintf('ASR CLEANING & ICA PREPARATION\n');
fprintf('========================================\n');
fprintf('Found %d files to process\n\n', length(files));

% Initialize logs
processLog = {};
errorLog = {};
qcMetrics = {};

%% =========================== MAIN LOOP ==================================

for file_i = 1:length(files)

    input_file = files(file_i).name;
    input_folder = files(file_i).folder;

    % Check if already processed
    if ~FORCE_REPROCESS
        output_pattern = strrep(input_file, '_preASR.set', '_clean.set');
        existingFile = dir(fullfile(input_folder, output_pattern));
        if ~isempty(existingFile)
            fprintf('[%d/%d] SKIP: %s (already processed)\n', ...
                file_i, length(files), input_file);
            continue;
        end
    end

    fprintf('\n[%d/%d] Processing: %s\n', file_i, length(files), input_file);

    try
        %% ==================== STEP 1: LOAD DATA =====================
        fprintf('  Step 1/6: Loading pre-ASR data...\n');

        EEG = pop_loadset(input_file, input_folder);
        EEG = eeg_checkset(EEG);

        % Retrieve stored ECG data
        if ~isfield(EEG.etc, 'ECG_data')
            error('ECG data not found in EEG.etc. File may be corrupted.');
        end
        ECG_data = EEG.etc.ECG_data;
        ECG_chanloc = EEG.etc.ECG_chanloc;

        % Store original for QC comparison
        if GENERATE_QC_PLOTS
            EEG_preASR = EEG;
        end

        %% ==================== STEP 2: SEPARATE vEOG =================
        fprintf('  Step 2/6: Separating vEOG channel...\n');

        % Find and remove vEOG channel (should not undergo ASR)
        veogIdx = find(strcmpi({EEG.chanlocs.type}, 'EOG'));

        if ~isempty(veogIdx)
            vEOG_data = EEG.data(veogIdx, :);
            vEOG_chanloc = EEG.chanlocs(veogIdx);
            EEG = pop_select(EEG, 'nochannel', veogIdx);
            EEG = eeg_checkset(EEG);
        else
            warning('vEOG channel not found - skipping vEOG separation');
            vEOG_data = [];
        end

        % Keep original channel locations for interpolation
        if isfield(EEG.etc, 'original_chanlocs')
            original_chanlocs = EEG.etc.original_chanlocs;
            original_chanlocs(veogIdx) = [];  % Remove vEOG from original list
        else
            original_chanlocs = EEG.chanlocs;
        end

        %% ==================== STEP 3: APPLY ASR =====================
        fprintf('  Step 3/6: Applying ASR (Burst=%d, Channel=%.2f)...\n', ...
            BURST_CRITERION, CHANNEL_CRITERION);

        % Apply ASR
        [EEG, ~, ~, removed_channels] = clean_artifacts(EEG, ...
            'FlatlineCriterion', FLATLINE_CRITERION, ...
            'ChannelCriterion', CHANNEL_CRITERION, ...
            'LineNoiseCriterion', LINE_NOISE_CRITERION, ...
            'BurstCriterion', BURST_CRITERION, ...
            'WindowCriterion', WINDOW_CRITERION, ...
            'Highpass', 'off');  % Already filtered in pre-ASR stage

        EEG = eeg_checkset(EEG);

        % Log removed channels
        removed_channel_labels = {original_chanlocs(removed_channels).labels};
        n_removed = sum(removed_channels);

        fprintf('    Removed %d channels: %s\n', n_removed, ...
            strjoin(removed_channel_labels, ', '));

        % Store ASR info
        EEG.etc.asr_info.removed_channels = removed_channel_labels;
        EEG.etc.asr_info.n_removed_channels = n_removed;
        EEG.etc.asr_info.burst_criterion = BURST_CRITERION;
        EEG.etc.asr_info.channel_criterion = CHANNEL_CRITERION;

        %% ==================== STEP 4: INTERPOLATE ===================
        fprintf('  Step 4/6: Interpolating removed channels...\n');

        % Interpolate removed channels back to original montage
        if n_removed > 0
            EEG = pop_interp(EEG, original_chanlocs, 'spherical');
            EEG = eeg_checkset(EEG);
            fprintf('    Interpolated %d channels\n', n_removed);
        else
            fprintf('    No channels to interpolate\n');
        end

        %% ==================== STEP 5: RE-REFERENCE ==================
        fprintf('  Step 5/6: Re-referencing to average...\n');

        % Re-reference EEG channels to average
        EEG = pop_reref(EEG, []);
        EEG = eeg_checkset(EEG);

        %% ==================== STEP 6: ADD CHANNELS BACK =============
        fprintf('  Step 6/6: Adding vEOG and ECG channels back...\n');

        % Add vEOG back
        if ~isempty(vEOG_data)
            % Verify data length matches
            if size(vEOG_data, 2) == size(EEG.data, 2)
                EEG.data(end+1, :) = vEOG_data;
                EEG.nbchan = EEG.nbchan + 1;
                EEG.chanlocs(end+1) = vEOG_chanloc;
            else
                warning('vEOG data length mismatch - NOT adding vEOG back');
            end
        end

        % Add ECG back
        if size(ECG_data, 2) == size(EEG.data, 2)
            EEG.data(end+1, :) = ECG_data;
            EEG.nbchan = EEG.nbchan + 1;
            EEG.chanlocs(end+1) = ECG_chanloc;
        else
            error('ECG data length mismatch - cannot continue');
        end

        EEG = eeg_checkset(EEG);

        %% ==================== QUALITY METRICS =======================

        % Calculate data quality metrics
        qc_idx = size(qcMetrics, 1) + 1;
        qcMetrics{qc_idx, 1} = EEG.subject;
        qcMetrics{qc_idx, 2} = EEG.session;
        qcMetrics{qc_idx, 3} = n_removed;
        qcMetrics{qc_idx, 4} = strjoin(removed_channel_labels, ';');
        qcMetrics{qc_idx, 5} = EEG.pnts / EEG.srate;  % Duration in seconds

        %% ==================== GENERATE QC PLOTS =====================

        if GENERATE_QC_PLOTS
            fprintf('  Generating QC plot...\n');

            % Create before/after PSD comparison
            fig = figure('Visible', 'off', 'Position', [100 100 1000 400]);

            % Exclude vEOG and ECG from PSD
            eeg_channels = 1:(EEG.nbchan-2);

            % Before ASR
            subplot(1, 2, 1);
            [psd_pre, freqs] = pwelch(EEG_preASR.data(eeg_channels,:)', [], [], [], EEG.srate);
            freq_mask = freqs <= PSD_MAX_FREQ;
            semilogy(freqs(freq_mask), mean(psd_pre(freq_mask,:), 2), 'k', 'LineWidth', 1.5);
            title('Before ASR');
            xlabel('Frequency (Hz)');
            ylabel('Power (µV²/Hz)');
            grid on;
            xlim([0 PSD_MAX_FREQ]);
            ylim([1e-2 1e4]);

            % After ASR
            subplot(1, 2, 2);
            [psd_post, ~] = pwelch(EEG.data(eeg_channels,:)', [], [], [], EEG.srate);
            semilogy(freqs(freq_mask), mean(psd_post(freq_mask,:), 2), 'b', 'LineWidth', 1.5);
            title(sprintf('After ASR (%d ch removed)', n_removed));
            xlabel('Frequency (Hz)');
            ylabel('Power (µV²/Hz)');
            grid on;
            xlim([0 PSD_MAX_FREQ]);
            ylim([1e-2 1e4]);

            sgtitle(sprintf('%s - PSD Comparison', EEG.setname), 'Interpreter', 'none');

            % Save plot
            plot_filename = sprintf('%s_PSD_ASR_comparison.png', EEG.setname(1:15));
            saveas(fig, fullfile(input_folder, plot_filename));
            close(fig);
        end

        %% ==================== SAVE CLEANED DATA =====================

        % Create output filename
        output_filename = strrep(input_file, '_preASR.set', '_clean.set');

        % Clean up temporary fields
        if isfield(EEG.etc, 'ECG_data')
            EEG.etc = rmfield(EEG.etc, 'ECG_data');
        end

        % Save
        EEG = pop_saveset(EEG, output_filename, input_folder);

        fprintf('  ✓ Saved: %s\n', output_filename);

        % Log success
        processLog{end+1, 1} = EEG.subject;
        processLog{end, 2} = EEG.session;
        processLog{end, 3} = n_removed;
        processLog{end, 4} = 'SUCCESS';

    catch ME
        % Log error
        errorLog{end+1, 1} = input_file;
        errorLog{end, 2} = ME.identifier;
        errorLog{end, 3} = ME.message;

        fprintf('  ✗ ERROR: %s\n', ME.message);
        continue;
    end

end % file loop

%% =========================== SUMMARY ====================================

fprintf('\n========================================\n');
fprintf('PROCESSING COMPLETE\n');
fprintf('========================================\n');

if ~isempty(processLog)
    processTable = cell2table(processLog, ...
        'VariableNames', {'Subject', 'Session', 'ChRemoved', 'Status'});
    fprintf('\nSuccessfully processed %d files:\n', size(processLog, 1));
    disp(processTable);

    % Summary statistics
    fprintf('\nChannel Removal Statistics:\n');
    fprintf('  Mean channels removed: %.1f\n', mean([processTable.ChRemoved]));
    fprintf('  Median channels removed: %.0f\n', median([processTable.ChRemoved]));
    fprintf('  Max channels removed: %d\n', max([processTable.ChRemoved]));
end

if ~isempty(qcMetrics)
    qcTable = cell2table(qcMetrics, ...
        'VariableNames', {'Subject', 'Session', 'ChRemoved', 'RemovedLabels', 'DurationSec'});

    % Save QC metrics
    writetable(qcTable, fullfile(DATA_PATH, 'ASR_quality_metrics.csv'));
    fprintf('\n✓ QC metrics saved to: ASR_quality_metrics.csv\n');
end

if ~isempty(errorLog)
    errorTable = cell2table(errorLog, ...
        'VariableNames', {'Filename', 'ErrorID', 'Message'});
    fprintf('\n✗ Errors (%d):\n', size(errorLog, 1));
    disp(errorTable);

    writetable(errorTable, fullfile(DATA_PATH, 'ASR_processing_errors.csv'));
end

fprintf('\nFiles ready for ICA/AMICA!\n');
fprintf('Next step: Run AMICA, then use RS6_IClabel.m\n\n');
