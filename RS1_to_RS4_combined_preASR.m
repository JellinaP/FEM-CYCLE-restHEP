% =========================================================================
% COMBINED PREPROCESSING: Import through Pre-ASR Cleaning
% =========================================================================
% This script combines RS1 (Import/Resample), RS3 (Trimming), and
% RS4 (CleanLine/Filtering) into a single pipeline that prepares data
% for ASR parameter tuning.
%
% OUTPUT: Files ending in '_preASR.set' ready for ASR parameter testing
% =========================================================================

clear; clc;

%% =========================== PARAMETERS =================================

% --- Paths ---
RAW_DATA_PATH = '/data/p_02938/FEM-CYCLE/NC_group/';
OUTPUT_PATH = '/data/pt_02938/rsHEP_analysis/0_data/';
EEGLAB_PATH = '/data/u_prinsen_software/eeglab2024.2.1';
CHANLOCS_FILE = '/data/p_02938/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for BrainAmpDC/CACS-64/CACS-64_REF.bvef';
CYCLE_TABLE_FILE = fullfile(RAW_DATA_PATH, 'participant_cycle_phases.xlsx');

% --- Channel Configuration ---
N_EEG_CHANNELS = 64;          % Number of EEG channels
ECG_CHANNEL_IDX = 65;         % ECG is channel 65 in bv data
RESP_CHANNEL_IDX = 66;        % RESP is channel 66 in bv data files
VEOG_CHANNEL_IDX = 21;        % vEOG is at TP10 (channel 21)
VEOG_CHANNEL_LABEL = 'TP10';  % Label of vEOG channel

% --- Output Folder Names ---
ECG_FOLDER_NAME = 'ecg';      % Subfolder for exported ECG files
RESP_FOLDER_NAME = 'resp';    % Subfolder for exported respiration files (if present)

% --- Resampling ---
TARGET_SRATE = 500;  % Hz

% --- Trimming ---
EXPECTED_DURATION = 600;    % seconds (10 minutes)
DURATION_TOLERANCE = 10;    % seconds (warn if duration differs by more than this)
START_MARKERS = {'S200', 'S 50', 'S  8', 'S200', 'S192'};  % Priority order for start markers
STOP_MARKERS  = {'S201', 'S 50', 'S  9', 'S200', 'S192'};  % Priority order for stop markers

% --- CleanLine Parameters ---
CLEANLINE_FREQS = [50 100];      % Line noise frequencies (Hz)
CLEANLINE_BANDWIDTH = 1;         % Bandwidth around each frequency
CLEANLINE_WIN_LENGTH = 4;        % Sliding window length (seconds)
CLEANLINE_WIN_STEP = 2;          % Sliding window step (seconds)
CLEANLINE_PADDING = 2;           % Padding factor

% --- Filtering Parameters ---
HIGHPASS_CUTOFF = 0.5;  % Hz (NOTE: Consider using 1.0 Hz if running ICA after)
LOWPASS_CUTOFF = 30;    % Hz

% --- Processing Options ---
SKIP_VISUAL_CHECK = true;  % Set to false to enable interactive visual checking
FORCE_REPROCESS = false;   % Set to true to reprocess already-processed files

%% =========================== INITIALIZE =================================

addpath(EEGLAB_PATH);
addpath(fullfile(EEGLAB_PATH, 'plugins', 'Cleanline2.1'));
eeglab('nogui');

% Load cycle phase information
cycleTable = readtable(CYCLE_TABLE_FILE);

% Get list of subject folders
cd(RAW_DATA_PATH);
folders = dir;
folders = folders([folders.isdir] & ~ismember({folders.name}, {'.', '..'}));
folders = folders(startsWith({folders.name}, 'NC'));

% Initialize logs
processLog = {};
errorLog = {};
warningLog = {};

fprintf('\n========================================\n');
fprintf('COMBINED PREPROCESSING (Pre-ASR)\n');
fprintf('========================================\n');
fprintf('Processing %d subjects\n', length(folders));
fprintf('Output suffix: _preASR.set\n\n');

%% =========================== MAIN LOOP ==================================

for subj_i = 1:length(folders)

    subject_folder = folders(subj_i).name;
    sessions = dir(fullfile(RAW_DATA_PATH, subject_folder, '*ses*'));

    for sess_i = 1:numel(sessions)

        session_folder = sessions(sess_i).name;
        output_folder = fullfile(OUTPUT_PATH, subject_folder, session_folder);

        % Check if already processed
        if ~FORCE_REPROCESS
            % Check for main EEG file
            existingFile = dir(fullfile(output_folder, '*resampled.set'));
            eegExists = ~isempty(existingFile);

            % Check for ECG export
            ecgFolder = fullfile(output_folder, ECG_FOLDER_NAME);
            ecgFile = dir(fullfile(ecgFolder, '*_ECG.set'));
            ecgExists = ~isempty(ecgFile);

            % Check for RESP export (optional)
            respFolder = fullfile(output_folder, RESP_FOLDER_NAME);
            respFile = dir(fullfile(respFolder, '*_RESP.set'));
            respExists = ~isempty(respFile);

            % Determine what needs processing
            needProcessEEG = ~eegExists;
            needProcessECG = ~ecgExists;
            needProcessRESP = ~respExists;

            % Skip only if ALL outputs exist
            if ~needProcessEEG && ~needProcessECG && ~needProcessRESP
                fprintf('[SKIP] %s/%s - already processed\n', subject_folder, session_folder);
                continue;
            end
        else
            % Force reprocess everything
            needProcessEEG = true;
            needProcessECG = true;
            needProcessRESP = true;
            ecgFolder = fullfile(output_folder, ECG_FOLDER_NAME);
            respFolder = fullfile(output_folder, RESP_FOLDER_NAME);
        end

        try
            fprintf('\n[PROCESSING] %s/%s\n', subject_folder, session_folder);

            % Find .vhdr file
            vhdrFile = dir(fullfile(RAW_DATA_PATH, subject_folder, session_folder, ...
                'task_REST', 'eeg', '*.vhdr'));

            if isempty(vhdrFile)
                error('No .vhdr file found');
            end

            %% ==================== STEP 1: IMPORT & RESAMPLE =============
            fprintf('  Step 1/7: Import and resample...\n');

            % Load BrainVision data (64 EEG + 1 ECG + 1 RESP)
            EEG = pop_loadbv(vhdrFile.folder, vhdrFile.name, [], []);

            % Extract metadata from filename
            nameParts = strsplit(vhdrFile.name, '_');
            EEG.subject = erase(nameParts{1}, 'sub-');
            EEG.session = erase(nameParts{2}, 'ses-');
            EEG.group = 'NC';

            % Get menstrual cycle phase from table
            row = strcmp(cycleTable.Subject, EEG.subject) & ...
                  cycleTable.Session == str2double(EEG.session);

            if any(row)
                EEG.condition = cycleTable.Phase{row};
            else
                warning('No cycle phase info for %s session %s', EEG.subject, EEG.session);
                warningLog{end+1, 1} = subject_folder;
                warningLog{end, 2} = session_folder;
                warningLog{end, 3} = 'Missing cycle phase info';
                EEG.condition = 'Unknown';
            end

            EEG.setname = strjoin({EEG.subject, nameParts{2}, EEG.condition(1:3)}, '_');
            EEG.filename = vhdrFile.name;
            EEG.filepath = vhdrFile.folder;

            % Load and assign channel locations
            % Load 64-channel EEG montage, then extend to 66 channels
            chanlocs = loadbvef(CHANLOCS_FILE);
            chanlocs = chanlocs(3:end); % Remove ref and ground (leaves 64 EEG channels)
            EEG.chanlocs = chanlocs;

            % Set types for EEG channels (1-64)
            for ch = 1:N_EEG_CHANNELS
                EEG.chanlocs(ch).type = 'EEG';
            end
            
            % Add ECG and RESP channel
            EEG.chanlocs(ECG_CHANNEL_IDX).type = 'ECG';
            EEG.chanlocs(ECG_CHANNEL_IDX).labels = 'ECG';
            EEG.chanlocs(RESP_CHANNEL_IDX).type = 'RESP';
            EEG.chanlocs(RESP_CHANNEL_IDX).labels = 'RESP';

            % Mark TP10 (ch 21) as vEOG
            EEG.chanlocs(VEOG_CHANNEL_IDX).type = 'EOG';
            EEG = eeg_checkset(EEG);

            % Resample to target sampling rate
            if EEG.srate ~= TARGET_SRATE
                EEG = pop_resample(EEG, TARGET_SRATE);
            end

            %% ==================== STEP 2: VISUAL CHECK ==================
            % (Optional - can be done manually between processing stages)

            if ~SKIP_VISUAL_CHECK
                fprintf('  Step 2/7: Visual inspection (close window when done)...\n');
                pop_eegplot(EEG, 1, 1, 1, 1, 'winlength', 50);
                uiwait;

                % User may have marked bad segments - this will be in EEG.reject
            end

            %% ==================== STEP 3: TRIMMING ======================
            fprintf('  Step 3/7: Trimming to %d seconds...\n', EXPECTED_DURATION);

            % Find start marker (try in priority order)
            % IMPORTANT: Use 'last' to get the LATEST start marker in case of multiple S200s
            startIdx = [];
            for m = 1:length(START_MARKERS)
                allStartIdx = find(strcmp({EEG.event.type}, START_MARKERS{m}));
                if ~isempty(allStartIdx)
                    startIdx = allStartIdx(end);  % Take the LATEST (last) start marker

                    % Warn if multiple start markers found
                    if length(allStartIdx) > 1
                        latencies = [EEG.event(allStartIdx).latency] / EEG.srate;
                        fprintf('  ⚠ Found %d start markers (%s) at times: %s seconds\n', ...
                                length(allStartIdx), START_MARKERS{m}, ...
                                mat2str(round(latencies, 1)));
                        fprintf('  → Using the LATEST one at %.1f seconds\n', latencies(end));

                        % Log this as a warning
                        warningLog{end+1, 1} = subject_folder;
                        warningLog{end, 2} = session_folder;
                        warningLog{end, 3} = sprintf('Multiple start markers (%d found), used latest', length(allStartIdx));
                    end
                    break;
                end
            end

            % Find stop marker (try in priority order)
            stopIdx = [];
            for m = 1:length(STOP_MARKERS)
                stopIdx = find(strcmp({EEG.event.type}, STOP_MARKERS{m}), 1, 'last');
                if ~isempty(stopIdx)
                    break;
                end
            end

            % Validate markers found
            if isempty(startIdx) || isempty(stopIdx)
                error('Could not find valid start/stop markers');
            end

            % Get timing
            startTime = EEG.event(startIdx).latency;
            stopTime = EEG.event(stopIdx).latency;
            duration = (stopTime - startTime) / EEG.srate;

            % Warn if duration is unexpected
            if abs(duration - EXPECTED_DURATION) > DURATION_TOLERANCE
                warning('Duration is %.1f s (expected %d s)', duration, EXPECTED_DURATION);
                warningLog{end+1, 1} = subject_folder;
                warningLog{end, 2} = session_folder;
                warningLog{end, 3} = sprintf('Duration %.1f s', duration);
            end

            % Trim data
            EEG = pop_select(EEG, 'point', [startTime stopTime]);
            EEG = eeg_checkset(EEG);

            %% ==================== STEP 4: EXTRACT ECG ==================
            fprintf('  Step 4/7: Extracting ECG channel...\n');

            % Find ECG channel
            ecgIdx = find(strcmpi({EEG.chanlocs.labels}, 'ECG'));

            if isempty(ecgIdx)
                error('ECG channel not found');
            end

            % Extract ECG data (trimmed, no filtering applied yet)
            ECG_data = EEG.data(ecgIdx, :);
            ECG_chanloc = EEG.chanlocs(ecgIdx);

            %% ==================== STEP 5: EXPORT ECG TO FILE ============
            if needProcessECG
                fprintf('  Step 5/7: Exporting ECG to separate file...\n');

                % Create ECG EEGLAB structure (full structure, not just array)
                EEG_ECG = EEG;  % Copy structure to preserve events and metadata
                EEG_ECG.data = ECG_data;
                EEG_ECG.nbchan = 1;
                EEG_ECG.chanlocs = ECG_chanloc;
                EEG_ECG.setname = [EEG.setname '_ECG'];
                EEG_ECG = eeg_checkset(EEG_ECG);

                % Create ECG folder if needed
                if ~exist(ecgFolder, 'dir')
                    mkdir(ecgFolder);
                end

                % Save ECG to separate folder
                ecgFilename = [EEG.setname(1:15), '_ECG.set'];
                pop_saveset(EEG_ECG, ecgFilename, ecgFolder);
                fprintf('  ✓ ECG exported to: %s/%s\n', ECG_FOLDER_NAME, ecgFilename);
            else
                fprintf('  Step 5/7: ECG already exported (skipped)\n');
            end

            %% ==================== STEP 6: EXPORT RESP (if exists) =======
            % Check if respiration channel exists (beyond standard EEG+vEOG+ECG)
            respIdx = find(strcmpi({EEG.chanlocs.labels}, 'RESP') | ...
                           strcmpi({EEG.chanlocs.type}, 'RESP'));

            hasRespChannel = ~isempty(respIdx);

            if hasRespChannel && needProcessRESP
                fprintf('  Step 6/7: Exporting RESP to separate file...\n');

                % Extract respiration data
                RESP_data = EEG.data(respIdx, :);
                RESP_chanloc = EEG.chanlocs(respIdx);

                % Create RESP EEGLAB structure
                EEG_RESP = EEG;  % Copy structure to preserve events
                EEG_RESP.data = RESP_data;
                EEG_RESP.nbchan = 1;
                EEG_RESP.chanlocs = RESP_chanloc;
                EEG_RESP.setname = [EEG.setname '_RESP'];
                EEG_RESP = eeg_checkset(EEG_RESP);

                % Create RESP folder if needed
                if ~exist(respFolder, 'dir')
                    mkdir(respFolder);
                end

                % Save RESP to separate folder
                respFilename = [EEG.setname(1:15), '_RESP.set'];
                pop_saveset(EEG_RESP, respFilename, respFolder);
                fprintf('  ✓ RESP exported to: %s/%s\n', RESP_FOLDER_NAME, respFilename);

                % Remove RESP from main EEG structure
                % Don't remove RESP yet - will be removed with ECG later if doing EEG processing
            elseif hasRespChannel && ~needProcessRESP
                fprintf('  Step 6/7: RESP already exported (skipped)\n');
            else
                fprintf('  Step 6/7: No RESP channel detected (typical for resting-state)\n');
            end

            %% ========= CONDITIONAL: CONTINUE WITH EEG PROCESSING? =======
            % If EEG already exists, stop here (ECG/RESP exported, done)
            if ~needProcessEEG
                fprintf('  ECG/RESP export complete. EEG already processed - skipping to next file.\n');

                % Log success
                processLog{end+1, 1} = subject_folder;
                processLog{end, 2} = session_folder;
                processLog{end, 3} = sprintf('%.1f s', duration);
                processLog{end, 4} = 'EEG_EXISTS';

                if needProcessECG
                    processLog{end, 5} = 'ECG_EXPORTED';
                else
                    processLog{end, 5} = 'ECG_EXISTS';
                end

                if hasRespChannel
                    if needProcessRESP
                        processLog{end, 6} = 'RESP_EXPORTED';
                    else
                        processLog{end, 6} = 'RESP_EXISTS';
                    end
                else
                    processLog{end, 6} = 'NO_RESP';
                end

                continue;  % Skip to next subject/session
            end

            %% ========== CONTINUE EEG PREPROCESSING (Steps 7+) ===========
            fprintf('  Continuing with EEG preprocessing...\n');

            % Remove ECG and RESP from EEG structure (won't undergo CleanLine/filtering)
            channelsToRemove = ecgIdx;
            if hasRespChannel
                channelsToRemove = [channelsToRemove, respIdx];
            end
            EEG = pop_select(EEG, 'nochannel', channelsToRemove);
            EEG = eeg_checkset(EEG);

            % Store original channel info for later interpolation
            EEG.etc.original_chanlocs = EEG.chanlocs;
            EEG.etc.veog_channel_idx = VEOG_CHANNEL_IDX;

            %% ==================== STEP 7: CLEANLINE & FILTERING =========
            fprintf('  Step 7/7: CleanLine and filtering...\n');

            % Apply CleanLine to remove line noise (50/100 Hz)
            % NOTE: This is applied to all channels including vEOG
            % Rationale: Line noise affects all channels similarly
            EEG = pop_cleanline(EEG, ...
                'Bandwidth', CLEANLINE_BANDWIDTH, ...
                'ChanCompIndices', 1:EEG.nbchan, ...
                'SignalType', 'Channels', ...
                'LineFrequencies', CLEANLINE_FREQS, ...
                'PaddingFactor', CLEANLINE_PADDING, ...
                'ScanForLines', true, ...
                'VerboseOutput', false, ...
                'SlidingWinLength', CLEANLINE_WIN_LENGTH, ...
                'SlidingWinStep', CLEANLINE_WIN_STEP);
            EEG = eeg_checkset(EEG);

            % Apply band-pass filter
            % NOTE: This is applied to all channels including vEOG
            % For different filter settings per channel type, separate before this step
            EEG = pop_eegfiltnew(EEG, 'locutoff', HIGHPASS_CUTOFF, 'hicutoff', LOWPASS_CUTOFF);
            EEG = eeg_checkset(EEG);

            %% ==================== STEP 8: SAVE PRE-ASR FILE =============

            % Create output folder if needed
            if ~exist(output_folder, 'dir')
                mkdir(output_folder);
            end

            % Store ECG and processing parameters in EEG structure
            EEG.etc.ECG_data = ECG_data;
            EEG.etc.ECG_chanloc = ECG_chanloc;
            EEG.etc.preprocessing_params.resample_rate = TARGET_SRATE;
            EEG.etc.preprocessing_params.highpass = HIGHPASS_CUTOFF;
            EEG.etc.preprocessing_params.lowpass = LOWPASS_CUTOFF;
            EEG.etc.preprocessing_params.cleanline_freqs = CLEANLINE_FREQS;
            EEG.etc.preprocessing_params.duration_seconds = EEG.pnts / EEG.srate;

            % Save with _preASR suffix
            output_filename = [EEG.setname(1:15), '_preASR.set'];
            EEG = pop_saveset(EEG, output_filename, output_folder);
            fprintf('  ✓ Saved: %s\n', output_filename);

            % Log success with detailed status
            % (We only reach here if needProcessEEG was true - early exit handles the other case)
            processLog{end+1, 1} = subject_folder;
            processLog{end, 2} = session_folder;
            processLog{end, 3} = sprintf('%.1f s', EEG.pnts/EEG.srate);
            processLog{end, 4} = 'EEG_SAVED';

            if needProcessECG
                processLog{end, 5} = 'ECG_EXPORTED';
            else
                processLog{end, 5} = 'ECG_EXISTS';
            end

            if hasRespChannel
                if needProcessRESP
                    processLog{end, 6} = 'RESP_EXPORTED';
                else
                    processLog{end, 6} = 'RESP_EXISTS';
                end
            else
                processLog{end, 6} = 'NO_RESP';
            end

        catch ME
            % Log error
            errorLog{end+1, 1} = subject_folder;
            errorLog{end, 2} = session_folder;
            errorLog{end, 3} = ME.message;

            fprintf('  ✗ ERROR: %s\n', ME.message);
            continue;
        end

    end % sessions
end % subjects

%% =========================== SUMMARY ====================================

fprintf('\n========================================\n');
fprintf('PROCESSING COMPLETE\n');
fprintf('========================================\n');

if ~isempty(processLog)
    processTable = cell2table(processLog, ...
        'VariableNames', {'Subject', 'Session', 'Duration', 'EEG_Status', 'ECG_Status', 'RESP_Status'});
    fprintf('\nSuccessfully processed %d files:\n', size(processLog, 1));
    disp(processTable);

    % Summary statistics
    fprintf('\nExport Summary:\n');
    fprintf('  ECG exported: %d files\n', sum(strcmp(processLog(:,5), 'ECG_EXPORTED')));
    fprintf('  RESP exported: %d files\n', sum(strcmp(processLog(:,6), 'RESP_EXPORTED')));
    fprintf('  No RESP channel: %d files\n', sum(strcmp(processLog(:,6), 'NO_RESP')));
end

if ~isempty(warningLog)
    warningTable = cell2table(warningLog, ...
        'VariableNames', {'Subject', 'Session', 'Warning'});
    fprintf('\n⚠ Warnings (%d):\n', size(warningLog, 1));
    disp(warningTable);
end

if ~isempty(errorLog)
    errorTable = cell2table(errorLog, ...
        'VariableNames', {'Subject', 'Session', 'Error'});
    fprintf('\n✗ Errors (%d):\n', size(errorLog, 1));
    disp(errorTable);

    % Save error log
    writetable(errorTable, fullfile(OUTPUT_PATH, 'preprocessing_errors.csv'));
end

fprintf('\nReady for ASR parameter tuning!\n');
fprintf('Next: Use ASR_parameter_tuning_QC.m to optimize ASR settings\n\n');
