% ==============================================================
% RS8_UNIFIED: MULTI-PIPELINE HEP STACKING WITH DUAL OUTPUT STRUCTURES
% ==============================================================
% This script performs FOUR HEP analysis pipelines with TWO output structures:
%
% OUTPUT 1: ft_HEPdata (PAIRED SUBJECTS ONLY - for cluster permutations)
% -----------------------------------------------------------------------
% - Only subjects with BOTH follicular AND luteal sessions
% - Perfect alignment: ft_HEPdata.*.follicular{i} ↔ ft_HEPdata.*.luteal{i}
% - Contains 4 pipelines:
%   .raw          - No adjustments
%   .IBIadjusted  - Trial-level IBI regression
%   .HRadjusted   - Subject-level HR residualization
%   .HRVadjusted  - Subject-level HRV (RMSSD) residualization
% - Compatible with depsamplesT design matrix for ft_timelockstatistics
%
% OUTPUT 2: ft_HEPdata_all (ALL SUBJECTS - for exploratory analyses)
% -------------------------------------------------------------------
% - Includes ALL subjects regardless of phase availability
% - Missing phases filled with NaN placeholders
% - Contains ONLY raw data (no residualization)
% - Useful for descriptive stats, single-phase analyses
%
% PROCESSING STRATEGY:
% --------------------
% Phase 1A: Build subject registry (two-pass scan)
% Phase 1B: Process paired subjects (strict alignment, all 4 pipelines)
% Phase 1C: Process all subjects (raw data only, with NaN placeholders)
% Phase 2:  Subject-level residualization (paired subjects only)
% Phase 3:  Save dual outputs with validation
% ==============================================================

clc; clear;

addpath /data/u_prinsen_software/eeglab2024.2.1
addpath /data/u_prinsen_software/fieldtrip-20250507
ft_defaults
eeglab nogui;

basePath = '/data/pt_02938/rsHEP_analysis/';

% --- Settings ---
icaType = 'icCFA';            % 'icCFA' or 'iclabel'
phaseTable = readtable('/data/p_02938/FEM-CYCLE/NC_group/participant_cycle_phases.xlsx');
hrvFile = fullfile(basePath, '4_results', 'RSEEG_HRV_features.xlsx');

subDirs = dir(fullfile(basePath, '0_data', 'NC*'));

fprintf('========================================\n');
fprintf('MULTI-PIPELINE HEP STACKING\n');
fprintf('========================================\n\n');

% ==============================================================
% PHASE 1A: BUILD SUBJECT REGISTRY (TWO-PASS SCAN)
% ==============================================================

fprintf('PHASE 1A: Scanning subjects and sessions...\n');
fprintf('========================================\n');

subjRegistry = struct();

for i = 1:numel(subDirs)
    subjID = subDirs(i).name;
    sessions = dir(fullfile(subDirs(i).folder, subjID, 'ses*'));

    for s = 1:numel(sessions)
        % Check if HEP file exists
        eegFile = dir(fullfile(subDirs(i).folder, subjID, ...
            sessions(s).name, ['*_HEP_' icaType '.set']));

        if isempty(eegFile)
            continue;
        end

        % Check if phase mapping exists
        matchIdx = strcmp(phaseTable.Subject, subjID) & phaseTable.Session == s;
        if ~any(matchIdx)
            continue;
        end

        phase = phaseTable.Phase{matchIdx};

        % Register this subject-phase combination
        if ~isfield(subjRegistry, subjID)
            subjRegistry.(subjID) = struct();
        end
        subjRegistry.(subjID).(phase).sessionNum = s;
        subjRegistry.(subjID).(phase).sessionDir = fullfile(subDirs(i).folder, subjID, sessions(s).name);
        subjRegistry.(subjID).(phase).eegFile = eegFile.name;
        subjRegistry.(subjID).(phase).exists = true;
    end
end

% Identify paired vs. unpaired subjects
allSubjects = fieldnames(subjRegistry);
pairedSubjects = {};
unpairedSubjects = {};

for i = 1:numel(allSubjects)
    sid = allSubjects{i};
    hasFol = isfield(subjRegistry.(sid), 'follicular');
    hasLut = isfield(subjRegistry.(sid), 'luteal');

    if hasFol && hasLut
        pairedSubjects{end+1} = sid; %#ok<SAGROW>
    else
        unpairedSubjects{end+1} = sid; %#ok<SAGROW>
    end
end

% Sort for consistent ordering
pairedSubjects = sort(pairedSubjects);
allSubjects = sort(allSubjects);

nPaired = numel(pairedSubjects);
nAll = numel(allSubjects);
nUnpaired = numel(unpairedSubjects);

fprintf('Summary:\n');
fprintf('  Total subjects scanned: %d\n', nAll);
fprintf('  Paired (FOL+LUT):       %d\n', nPaired);
fprintf('  Unpaired (FOL or LUT):  %d\n', nUnpaired);
if nUnpaired > 0
    fprintf('  Unpaired subjects: %s\n', strjoin(unpairedSubjects, ', '));
end
fprintf('\n');

% ==============================================================
% PHASE 1B: PROCESS PAIRED SUBJECTS (STRICT ALIGNMENT)
% ==============================================================

fprintf('========================================\n');
fprintf('PHASE 1B: Processing PAIRED subjects\n');
fprintf('========================================\n');
fprintf('Processing %d paired subjects (all 4 pipelines)...\n\n', nPaired);

% Pre-allocate structures for paired subjects
ft_HEPdata_raw_paired = struct();
ft_HEPdata_raw_paired.follicular = cell(nPaired, 1);
ft_HEPdata_raw_paired.luteal = cell(nPaired, 1);

ft_HEPdata_ibi_paired = struct();
ft_HEPdata_ibi_paired.follicular = cell(nPaired, 1);
ft_HEPdata_ibi_paired.luteal = cell(nPaired, 1);

% Process each paired subject in sorted order
for subjIdx = 1:nPaired
    subjID = pairedSubjects{subjIdx};

    fprintf('Subject %d/%d: %s\n', subjIdx, nPaired, subjID);

    % --------------------------------------------------------
    % Process FOLLICULAR session
    % --------------------------------------------------------
    phase = 'follicular';
    sessionNum = subjRegistry.(subjID).(phase).sessionNum;
    sessionDir = subjRegistry.(subjID).(phase).sessionDir;
    eegFileName = subjRegistry.(subjID).(phase).eegFile;

    EEG_fol = pop_loadset(eegFileName, sessionDir);

    % Pipeline 1: RAW
    [ft_avg_raw_fol] = process_raw_pipeline(EEG_fol);
    ft_HEPdata_raw_paired.follicular{subjIdx} = ft_avg_raw_fol;

    % Pipeline 2: IBI-adjusted
    [ft_avg_ibi_fol, ibi_mean_fol, nTrials_fol] = process_ibi_pipeline(EEG_fol, subjID, sessionNum);
    ft_HEPdata_ibi_paired.follicular{subjIdx} = ft_avg_ibi_fol;

    % --------------------------------------------------------
    % Process LUTEAL session
    % --------------------------------------------------------
    phase = 'luteal';
    sessionNum = subjRegistry.(subjID).(phase).sessionNum;
    sessionDir = subjRegistry.(subjID).(phase).sessionDir;
    eegFileName = subjRegistry.(subjID).(phase).eegFile;

    EEG_lut = pop_loadset(eegFileName, sessionDir);

    % Pipeline 1: RAW
    [ft_avg_raw_lut] = process_raw_pipeline(EEG_lut);
    ft_HEPdata_raw_paired.luteal{subjIdx} = ft_avg_raw_lut;

    % Pipeline 2: IBI-adjusted
    [ft_avg_ibi_lut, ibi_mean_lut, nTrials_lut] = process_ibi_pipeline(EEG_lut, subjID, sessionNum);
    ft_HEPdata_ibi_paired.luteal{subjIdx} = ft_avg_ibi_lut;

    fprintf('  FOL: %d trials, IBI=%.3f s | LUT: %d trials, IBI=%.3f s\n', ...
        nTrials_fol, ibi_mean_fol, nTrials_lut, ibi_mean_lut);
end

% Build ft_subjMeta for paired subjects
ft_subjMeta_paired = repelem(pairedSubjects(:)', 2);

fprintf('\nPHASE 1B COMPLETE: Paired subjects processed\n');
fprintf('  Paired subjects (RAW):         FOL=%d, LUT=%d\n', ...
    length(ft_HEPdata_raw_paired.follicular), length(ft_HEPdata_raw_paired.luteal));
fprintf('  Paired subjects (IBI-adjusted): FOL=%d, LUT=%d\n', ...
    length(ft_HEPdata_ibi_paired.follicular), length(ft_HEPdata_ibi_paired.luteal));
fprintf('\n');

% ==============================================================
% PHASE 1C: PROCESS ALL SUBJECTS (RAW DATA ONLY)
% ==============================================================

fprintf('========================================\n');
fprintf('PHASE 1C: Processing ALL subjects\n');
fprintf('========================================\n');
fprintf('Processing %d subjects (raw data only, with NaN placeholders)...\n\n', nAll);

% Pre-allocate structures for all subjects
ft_HEPdata_raw_all = struct();
ft_HEPdata_raw_all.follicular = cell(nAll, 1);
ft_HEPdata_raw_all.luteal = cell(nAll, 1);

% Get common dimensions from first paired subject for NaN placeholders
example_tl = ft_HEPdata_raw_paired.follicular{1};
nChan = size(example_tl.avg, 1);
commonTime = example_tl.time;
channelLabels = example_tl.label;
hasElec = isfield(example_tl, 'elec');
if hasElec
    elecStruct = example_tl.elec;
end

% Process each subject
for subjIdx = 1:nAll
    subjID = allSubjects{subjIdx};

    fprintf('Subject %d/%d: %s', subjIdx, nAll, subjID);

    % --------------------------------------------------------
    % Process FOLLICULAR session (if exists)
    % --------------------------------------------------------
    if isfield(subjRegistry.(subjID), 'follicular')
        sessionDir = subjRegistry.(subjID).follicular.sessionDir;
        eegFileName = subjRegistry.(subjID).follicular.eegFile;

        EEG = pop_loadset(eegFileName, sessionDir);
        ft_avg_raw_fol = process_raw_pipeline(EEG);
        ft_HEPdata_raw_all.follicular{subjIdx} = ft_avg_raw_fol;
        fprintf(' [FOL: OK]');
    else
        % Create NaN placeholder
        ft_HEPdata_raw_all.follicular{subjIdx} = create_nan_placeholder(nChan, commonTime, channelLabels, hasElec, elecStruct);
        fprintf(' [FOL: NaN]');
    end

    % --------------------------------------------------------
    % Process LUTEAL session (if exists)
    % --------------------------------------------------------
    if isfield(subjRegistry.(subjID), 'luteal')
        sessionDir = subjRegistry.(subjID).luteal.sessionDir;
        eegFileName = subjRegistry.(subjID).luteal.eegFile;

        EEG = pop_loadset(eegFileName, sessionDir);
        ft_avg_raw_lut = process_raw_pipeline(EEG);
        ft_HEPdata_raw_all.luteal{subjIdx} = ft_avg_raw_lut;
        fprintf(' [LUT: OK]\n');
    else
        % Create NaN placeholder
        ft_HEPdata_raw_all.luteal{subjIdx} = create_nan_placeholder(nChan, commonTime, channelLabels, hasElec, elecStruct);
        fprintf(' [LUT: NaN]\n');
    end
end

% Build ft_subjMeta_all for all subjects
ft_subjMeta_all = repelem(allSubjects(:)', 2);

fprintf('\nPHASE 1C COMPLETE: All subjects processed\n');
fprintf('  All subjects (RAW): FOL=%d, LUT=%d (includes NaN placeholders)\n', ...
    length(ft_HEPdata_raw_all.follicular), length(ft_HEPdata_raw_all.luteal));
fprintf('\n');

% ==============================================================
% PHASE 2: SUBJECT-LEVEL RESIDUALIZATION (PAIRED SUBJECTS ONLY)
% ==============================================================

fprintf('========================================\n');
fprintf('PHASE 2: Subject-Level Residualization\n');
fprintf('========================================\n');
fprintf('Operating on %d PAIRED subjects only\n\n', nPaired);

% --- Load HRV features ---
fprintf('Loading HRV features from:\n  %s\n', hrvFile);
if ~exist(hrvFile, 'file')
    error('HRV features file not found: %s', hrvFile);
end
hrvTable = readtable(hrvFile, 'Sheet', 'HRV');

% Extract subject IDs for paired subjects
subjIDs_paired = pairedSubjects(:);
nSubj = nPaired;

fprintf('  Loaded HRV features for residualization\n\n');

% --- Determine common time vector ---
tl_FOL = ft_HEPdata_raw_paired.follicular;
tl_LUT = ft_HEPdata_raw_paired.luteal;

nTimePerSubj = cellfun(@(s) numel(s.time), [tl_FOL(:); tl_LUT(:)]);
nTimeCommon  = min(nTimePerSubj);
commonTime   = tl_FOL{1}.time(1:nTimeCommon);
fprintf('Common time vector: %d samples (%.3f to %.3f s)\n', ...
    nTimeCommon, commonTime(1), commonTime(end));

% ==============================================================
% PIPELINE 3: HR-ADJUSTED
% ==============================================================
fprintf('\n--- Pipeline 3: HR Residualization ---\n');

% Build per-subject HR matrices for each phase
hrMat_fol = nan(nSubj, 1);
hrMat_lut = nan(nSubj, 1);

for s = 1:nSubj
    sid = subjIDs_paired{s};
    idx_fol = strcmp(hrvTable.subj, sid) & strcmp(hrvTable.phase, 'fol');
    idx_lut = strcmp(hrvTable.subj, sid) & strcmp(hrvTable.phase, 'lut');

    if any(idx_fol); hrMat_fol(s) = hrvTable.time_heart_rate(idx_fol); end
    if any(idx_lut); hrMat_lut(s) = hrvTable.time_heart_rate(idx_lut); end
end

fprintf('  HR availability: FOL=%d/%d, LUT=%d/%d\n', ...
    sum(~isnan(hrMat_fol)), nSubj, sum(~isnan(hrMat_lut)), nSubj);

% Residualize each phase
tl_FOL_hr = residualize_phase_data(tl_FOL, hrMat_fol, 'HR-FOL', nSubj, nChan, nTimeCommon, commonTime);
tl_LUT_hr = residualize_phase_data(tl_LUT, hrMat_lut, 'HR-LUT', nSubj, nChan, nTimeCommon, commonTime);

ft_HEPdata_hr_paired = struct();
ft_HEPdata_hr_paired.follicular = tl_FOL_hr;
ft_HEPdata_hr_paired.luteal     = tl_LUT_hr;

% ==============================================================
% PIPELINE 4: HRV-ADJUSTED (RMSSD)
% ==============================================================
fprintf('\n--- Pipeline 4: HRV (RMSSD) Residualization ---\n');

% Build per-subject RMSSD matrices for each phase
hrvMat_fol = nan(nSubj, 1);
hrvMat_lut = nan(nSubj, 1);

for s = 1:nSubj
    sid = subjIDs_paired{s};
    idx_fol = strcmp(hrvTable.subj, sid) & strcmp(hrvTable.phase, 'fol');
    idx_lut = strcmp(hrvTable.subj, sid) & strcmp(hrvTable.phase, 'lut');

    if any(idx_fol); hrvMat_fol(s) = hrvTable.time_RMSSD(idx_fol); end
    if any(idx_lut); hrvMat_lut(s) = hrvTable.time_RMSSD(idx_lut); end
end

fprintf('  RMSSD availability: FOL=%d/%d, LUT=%d/%d\n', ...
    sum(~isnan(hrvMat_fol)), nSubj, sum(~isnan(hrvMat_lut)), nSubj);

% Residualize each phase
tl_FOL_hrv = residualize_phase_data(tl_FOL, hrvMat_fol, 'RMSSD-FOL', nSubj, nChan, nTimeCommon, commonTime);
tl_LUT_hrv = residualize_phase_data(tl_LUT, hrvMat_lut, 'RMSSD-LUT', nSubj, nChan, nTimeCommon, commonTime);

ft_HEPdata_hrv_paired = struct();
ft_HEPdata_hrv_paired.follicular = tl_FOL_hrv;
ft_HEPdata_hrv_paired.luteal     = tl_LUT_hrv;

fprintf('\nPHASE 2 COMPLETE: Subject-level residualization done\n');
fprintf('\n');

% ==============================================================
% PHASE 3: SAVE DUAL OUTPUTS WITH VALIDATION
% ==============================================================

fprintf('========================================\n');
fprintf('PHASE 3: Saving Outputs\n');
fprintf('========================================\n\n');

outDir = fullfile(basePath, '04_results', '01_Phase_main', icaType);
if ~exist(outDir, 'dir'); mkdir(outDir); end

% --------------------------------------------------------
% OUTPUT 1: PAIRED SUBJECTS (all 4 pipelines)
% --------------------------------------------------------

fprintf('--- OUTPUT 1: Paired Subjects (for cluster permutations) ---\n');

% Validate pairing
nFOL_raw = length(ft_HEPdata_raw_paired.follicular);
nLUT_raw = length(ft_HEPdata_raw_paired.luteal);
nFOL_ibi = length(ft_HEPdata_ibi_paired.follicular);
nLUT_ibi = length(ft_HEPdata_ibi_paired.luteal);
nFOL_hr = length(ft_HEPdata_hr_paired.follicular);
nLUT_hr = length(ft_HEPdata_hr_paired.luteal);
nFOL_hrv = length(ft_HEPdata_hrv_paired.follicular);
nLUT_hrv = length(ft_HEPdata_hrv_paired.luteal);

if ~(nFOL_raw == nLUT_raw && nFOL_ibi == nLUT_ibi && nFOL_hr == nLUT_hr && nFOL_hrv == nLUT_hrv)
    error('Pairing validation FAILED: Follicular and luteal arrays have different lengths!');
end

fprintf('  Pairing validation: PASSED ✓\n');
fprintf('  All pipelines have %d paired subjects\n', nFOL_raw);

% Create unified structure
ft_HEPdata = struct();
ft_HEPdata.raw         = ft_HEPdata_raw_paired;
ft_HEPdata.IBIadjusted = ft_HEPdata_ibi_paired;
ft_HEPdata.HRadjusted  = ft_HEPdata_hr_paired;
ft_HEPdata.HRVadjusted = ft_HEPdata_hrv_paired;

ft_subjMeta = ft_subjMeta_paired;

outFile1 = fullfile(outDir, ['HEPs_' icaType '_all_versions.mat']);
save(outFile1, 'ft_HEPdata', 'ft_subjMeta', '-v7.3');

fprintf('\nSaved: %s\n', outFile1);
fprintf('Structure summary (PAIRED SUBJECTS):\n');
fprintf('  ft_HEPdata.raw.follicular:         %d subjects\n', length(ft_HEPdata.raw.follicular));
fprintf('  ft_HEPdata.raw.luteal:             %d subjects\n', length(ft_HEPdata.raw.luteal));
fprintf('  ft_HEPdata.IBIadjusted.follicular: %d subjects\n', length(ft_HEPdata.IBIadjusted.follicular));
fprintf('  ft_HEPdata.IBIadjusted.luteal:     %d subjects\n', length(ft_HEPdata.IBIadjusted.luteal));
fprintf('  ft_HEPdata.HRadjusted.follicular:  %d subjects\n', length(ft_HEPdata.HRadjusted.follicular));
fprintf('  ft_HEPdata.HRadjusted.luteal:      %d subjects\n', length(ft_HEPdata.HRadjusted.luteal));
fprintf('  ft_HEPdata.HRVadjusted.follicular: %d subjects\n', length(ft_HEPdata.HRVadjusted.follicular));
fprintf('  ft_HEPdata.HRVadjusted.luteal:     %d subjects\n', length(ft_HEPdata.HRVadjusted.luteal));
fprintf('  ft_subjMeta:                       %d entries (2 per subject)\n', length(ft_subjMeta));

% --------------------------------------------------------
% OUTPUT 2: ALL SUBJECTS (raw data only)
% --------------------------------------------------------

fprintf('\n--- OUTPUT 2: All Subjects (for exploratory analyses) ---\n');

ft_HEPdata_all = struct();
ft_HEPdata_all.raw = ft_HEPdata_raw_all;

ft_subjMeta_all_out = ft_subjMeta_all;

outFile2 = fullfile(outDir, ['HEPs_' icaType '_all_subjects_raw.mat']);
save(outFile2, 'ft_HEPdata_all', 'ft_subjMeta_all', '-v7.3');

fprintf('\nSaved: %s\n', outFile2);
fprintf('Structure summary (ALL SUBJECTS):\n');
fprintf('  ft_HEPdata_all.raw.follicular: %d subjects (includes NaN placeholders)\n', ...
    length(ft_HEPdata_all.raw.follicular));
fprintf('  ft_HEPdata_all.raw.luteal:     %d subjects (includes NaN placeholders)\n', ...
    length(ft_HEPdata_all.raw.luteal));
fprintf('  ft_subjMeta_all:               %d entries (2 per subject)\n', length(ft_subjMeta_all));

% --------------------------------------------------------
% Save electrode file for downstream scripts
% --------------------------------------------------------

if exist('EEG', 'var')
    exampleFile = fullfile(outDir, '64_channels_example.set');
    if ~exist(exampleFile, 'file')
        pop_saveset(EEG, '64_channels_example', outDir);
        fprintf('\n  Electrode file saved: 64_channels_example.set\n');
    end
end

fprintf('\n========================================\n');
fprintf('ALL PROCESSING COMPLETE\n');
fprintf('========================================\n');
fprintf('Summary:\n');
fprintf('  OUTPUT 1 (paired):  %d subjects × 4 pipelines → %s\n', nPaired, outFile1);
fprintf('  OUTPUT 2 (all):     %d subjects × raw only    → %s\n', nAll, outFile2);
fprintf('\nReady for downstream analyses:\n');
fprintf('  - Cluster-based permutations: Use OUTPUT 1\n');
fprintf('  - Exploratory / descriptive:  Use OUTPUT 2\n');
fprintf('========================================\n');

% ==============================================================
% NESTED FUNCTIONS
% ==============================================================

function [ft_avg_raw] = process_raw_pipeline(EEG)
    % Process raw HEP data without any adjustments
    ft_raw = eeglab2fieldtrip(EEG, 'raw');

    cfgTL = [];
    cfgTL.keeptrials = 'yes';
    ft_tl_raw = ft_timelockanalysis(cfgTL, ft_raw);

    ft_avg_raw       = [];
    ft_avg_raw.avg   = squeeze(mean(ft_tl_raw.trial, 1));
    ft_avg_raw.time  = ft_tl_raw.time;
    ft_avg_raw.label = ft_tl_raw.label;
    ft_avg_raw.dimord = 'chan_time';
    if isfield(EEG, 'chanlocs')
        ft_avg_raw.elec = ft_raw.elec;
    end
end

function [ft_avg_ibi, ibi_mean, nTrials] = process_ibi_pipeline(EEG, subjID, sessionNum)
    % Process HEP data with trial-level IBI regression

    % Extract per-epoch IBI from BrainBeats
    NN_times = EEG.brainbeats.preprocessings.NN_times;
    ibi_all  = diff(NN_times(:));

    nTrials = EEG.trials;
    if numel(ibi_all) ~= nTrials
        warning('%s ses-%d: IBI count (%d) ~= trial count (%d). Truncating to min.', ...
            subjID, sessionNum, numel(ibi_all), nTrials);
        nUse = min(numel(ibi_all), nTrials);
        ibi_all = ibi_all(1:nUse);
        EEG = pop_select(EEG, 'trial', 1:nUse);
        nTrials = nUse;
    end

    % Convert to FieldTrip
    ft_raw = eeglab2fieldtrip(EEG, 'raw');

    cfgTL = [];
    cfgTL.keeptrials = 'yes';
    ft_tl = ft_timelockanalysis(cfgTL, ft_raw);

    % Regress out IBI at trial level
    ibi_z   = zscore(ibi_all);
    confMat = [ones(nTrials, 1), ibi_z];

    cfgRC           = [];
    cfgRC.confound  = confMat;
    cfgRC.reject    = 2;
    cfgRC.normalize = 'no';
    cfgRC.output    = 'residual';
    ft_tl_clean = ft_regressconfound(cfgRC, ft_tl);

    % Average cleaned trials
    ft_avg_ibi       = [];
    ft_avg_ibi.avg   = squeeze(mean(ft_tl_clean.trial, 1));
    ft_avg_ibi.time  = ft_tl_clean.time;
    ft_avg_ibi.label = ft_tl_clean.label;
    ft_avg_ibi.dimord = 'chan_time';
    if isfield(EEG, 'chanlocs')
        ft_avg_ibi.elec = ft_raw.elec;
    end

    ibi_mean = mean(ibi_all);
end

function tl_nan = create_nan_placeholder(nChan, commonTime, channelLabels, hasElec, elecStruct)
    % Create NaN-filled placeholder timelock structure
    tl_nan = struct();
    tl_nan.avg = NaN(nChan, length(commonTime));
    tl_nan.time = commonTime;
    tl_nan.label = channelLabels;
    tl_nan.dimord = 'chan_time';
    if hasElec
        tl_nan.elec = elecStruct;
    end
end

function tl_out = residualize_phase_data(tl_in, covMat, covName, nSubj, nChan, nTimeCommon, commonTime)
    % Residualizes HEP data against a single covariate using ft_regressconfound
    % Keeps intercept to preserve signal level

    % Identify subjects with valid covariate data
    valid  = ~isnan(covMat);
    nValid = sum(valid);
    fprintf('    %s: %d/%d subjects with valid data\n', covName, nValid, nSubj);

    if nValid < 3
        warning('%s: too few valid subjects (%d). Returning NaN.', covName, nValid);
        tl_out = tl_in;
        for s = 1:nSubj
            tl_out{s}.avg  = nan(nChan, nTimeCommon);
            tl_out{s}.time = commonTime;
        end
        return;
    end

    % Z-score covariate across valid subjects
    covZ = zscore(covMat(valid));

    % Stack all subjects into 3-D array [nSubj x nChan x nTimeCommon]
    data_orig = nan(nSubj, nChan, nTimeCommon);
    for s = 1:nSubj
        data_orig(s, :, :) = tl_in{s}.avg(:, 1:nTimeCommon);
    end

    % Build a minimal FieldTrip struct for ft_regressconfound
    ga        = [];
    ga.trial  = data_orig(valid, :, :);
    ga.dimord = 'rpt_chan_time';
    ga.time   = commonTime;
    ga.label  = tl_in{1}.label;

    % Confound matrix: [intercept, covariate_z]
    confMat = [ones(nValid, 1), covZ];

    % Regress out covariate using ft_regressconfound
    cfgRC           = [];
    cfgRC.confound  = confMat;
    cfgRC.normalize = 'no';
    cfgRC.output    = 'residual';
    cfgRC.reject    = 2;             % reject covariate, keep intercept

    ga_res = ft_regressconfound(cfgRC, ga);

    % Unpack residualized data back into per-subject timelocks
    tl_out   = tl_in;
    validIdx = find(valid);
    for i = 1:nValid
        s = validIdx(i);
        tl_out{s}.avg  = squeeze(ga_res.trial(i, :, :));
        tl_out{s}.time = commonTime;
        if isfield(tl_out{s}, 'var')
            tl_out{s}.var = tl_out{s}.var(:, 1:nTimeCommon);
        end
    end

    % Set NaN for invalid subjects (missing covariates)
    for s = find(~valid)'
        tl_out{s}.avg  = nan(nChan, nTimeCommon);
        tl_out{s}.time = commonTime;
        if isfield(tl_out{s}, 'var')
            tl_out{s}.var = tl_out{s}.var(:, 1:nTimeCommon);
        end
    end

    fprintf('    %s: Residualization complete\n', covName);
end
