% =========================================================================
% BATCH AMICA PROCESSING
% =========================================================================
% Runs AMICA (Adaptive Mixture ICA) on cleaned EEG data.
%
% INPUT:  Files ending in '_clean.set'
% OUTPUT: AMICA folder in same directory as input file
%
% NOTES:
%   - ECG and vEOG channels are EXCLUDED from ICA
%   - AMICA runs on EEG channels (1-63) only
%   - This script can run locally (slow) or generate SLURM jobs (fast)
%   - For HPC: Use SLURM mode to parallelize across subjects
% =========================================================================

clear; clc;

%% =========================== PARAMETERS =================================

% --- Paths ---
DATA_PATH = '/data/pt_02938/rsHEP_analysis/0_data/';
EEGLAB_PATH = '/data/u_prinsen_software/eeglab2024.2.1';

% --- AMICA Parameters ---
NUM_MODELS = 1;             % Number of ICA models (1 for standard ICA)
MAX_ITERATIONS = 2000;      % Maximum iterations
NUM_THREADS = 8;            % Number of parallel threads (local mode)

% --- Processing Mode ---
MODE = 'generate_slurm';    % Options: 'local', 'generate_slurm'
                            % 'local' = run AMICA here (slow, one at a time)
                            % 'generate_slurm' = create SLURM job scripts

% --- SLURM Configuration (if MODE = 'generate_slurm') ---
SLURM_ACCOUNT = 'your_account';     % Your SLURM account
SLURM_PARTITION = 'medium';         % Partition name
SLURM_TIME = '04:00:00';            % Time limit per job (HH:MM:SS)
SLURM_MEM = '16G';                  % Memory per job
SLURM_CPUS = 8;                     % CPUs per job
SLURM_OUTPUT_DIR = fullfile(DATA_PATH, 'AMICA_logs');

% --- Processing Options ---
FORCE_REPROCESS = false;    % Set to true to rerun AMICA even if exists

%% =========================== INITIALIZE =================================

addpath(EEGLAB_PATH);
eeglab('nogui');

% Find all cleaned files
files = dir(fullfile(DATA_PATH, '**', '*_clean.set'));

if isempty(files)
    error('No *_clean.set files found! Run RS5_ASR_and_ICA_prep.m first.');
end

fprintf('\n========================================\n');
fprintf('AMICA BATCH PROCESSING\n');
fprintf('========================================\n');
fprintf('Mode: %s\n', MODE);
fprintf('Found %d files\n\n', length(files));

% Create log directory for SLURM if needed
if strcmp(MODE, 'generate_slurm') && ~exist(SLURM_OUTPUT_DIR, 'dir')
    mkdir(SLURM_OUTPUT_DIR);
end

%% =========================== PROCESS FILES ==============================

processLog = {};
slurm_scripts = {};

for file_i = 1:length(files)

    input_file = files(file_i).name;
    input_folder = files(file_i).folder;
    input_path = fullfile(input_folder, input_file);

    % Check if AMICA already exists
    amica_folder = fullfile(input_folder, 'AMICA');
    if exist(amica_folder, 'dir') && ~FORCE_REPROCESS
        fprintf('[%d/%d] SKIP: %s (AMICA folder exists)\n', ...
            file_i, length(files), input_file);
        continue;
    end

    fprintf('[%d/%d] Processing: %s\n', file_i, length(files), input_file);

    switch MODE
        case 'local'
            % Run AMICA locally
            try
                % Load data
                EEG = pop_loadset(input_file, input_folder);

                % Exclude ECG and vEOG channels
                ecg_idx = find(strcmpi({EEG.chanlocs.labels}, 'ECG'));
                eog_idx = find(strcmpi({EEG.chanlocs.type}, 'EOG'));
                exclude_idx = unique([ecg_idx, eog_idx]);

                if ~isempty(exclude_idx)
                    fprintf('  Excluding channels: %s\n', ...
                        strjoin({EEG.chanlocs(exclude_idx).labels}, ', '));
                    EEG_for_ica = pop_select(EEG, 'nochannel', exclude_idx);
                else
                    EEG_for_ica = EEG;
                end

                % Run AMICA
                fprintf('  Running AMICA (this may take 15-30 minutes)...\n');
                fprintf('  Output folder: %s\n', amica_folder);

                runamica15(EEG_for_ica.data, ...
                    'num_models', NUM_MODELS, ...
                    'outdir', amica_folder, ...
                    'max_iter', MAX_ITERATIONS, ...
                    'num_mix_comps', EEG_for_ica.nbchan, ...
                    'pcakeep', EEG_for_ica.nbchan, ...
                    'num_threads', NUM_THREADS);

                fprintf('  ✓ AMICA completed\n');

                processLog{end+1, 1} = input_file;
                processLog{end, 2} = 'SUCCESS';

            catch ME
                fprintf('  ✗ ERROR: %s\n', ME.message);
                processLog{end+1, 1} = input_file;
                processLog{end, 2} = ME.message;
            end

        case 'generate_slurm'
            % Generate SLURM job script
            job_name = sprintf('amica_%s', extractBefore(input_file, '_clean.set'));
            slurm_script_name = sprintf('%s.sh', job_name);
            slurm_script_path = fullfile(SLURM_OUTPUT_DIR, slurm_script_name);

            % Create SLURM script
            fid = fopen(slurm_script_path, 'w');

            fprintf(fid, '#!/bin/bash\n');
            fprintf(fid, '#SBATCH --job-name=%s\n', job_name);
            fprintf(fid, '#SBATCH --account=%s\n', SLURM_ACCOUNT);
            fprintf(fid, '#SBATCH --partition=%s\n', SLURM_PARTITION);
            fprintf(fid, '#SBATCH --time=%s\n', SLURM_TIME);
            fprintf(fid, '#SBATCH --mem=%s\n', SLURM_MEM);
            fprintf(fid, '#SBATCH --cpus-per-task=%d\n', SLURM_CPUS);
            fprintf(fid, '#SBATCH --output=%s\n', fullfile(SLURM_OUTPUT_DIR, sprintf('%s_%%j.out', job_name)));
            fprintf(fid, '#SBATCH --error=%s\n', fullfile(SLURM_OUTPUT_DIR, sprintf('%s_%%j.err', job_name)));
            fprintf(fid, '\n');
            fprintf(fid, '# Load MATLAB module (adjust for your system)\n');
            fprintf(fid, 'module load matlab/R2023b\n');
            fprintf(fid, '\n');
            fprintf(fid, '# Run MATLAB script\n');
            fprintf(fid, 'matlab -nodisplay -nosplash -r "');
            fprintf(fid, 'addpath(''%s''); ', EEGLAB_PATH);
            fprintf(fid, 'eeglab nogui; ');
            fprintf(fid, 'EEG = pop_loadset(''%s'', ''%s''); ', input_file, input_folder);
            fprintf(fid, 'ecg_idx = find(strcmpi({EEG.chanlocs.labels}, ''ECG'')); ');
            fprintf(fid, 'eog_idx = find(strcmpi({EEG.chanlocs.type}, ''EOG'')); ');
            fprintf(fid, 'exclude_idx = unique([ecg_idx, eog_idx]); ');
            fprintf(fid, 'if ~isempty(exclude_idx), EEG = pop_select(EEG, ''nochannel'', exclude_idx); end; ');
            fprintf(fid, 'runamica15(EEG.data, ''num_models'', %d, ', NUM_MODELS);
            fprintf(fid, '''outdir'', ''%s'', ', amica_folder);
            fprintf(fid, '''max_iter'', %d, ', MAX_ITERATIONS);
            fprintf(fid, '''num_mix_comps'', EEG.nbchan, ');
            fprintf(fid, '''pcakeep'', EEG.nbchan, ');
            fprintf(fid, '''num_threads'', %d); ', SLURM_CPUS);
            fprintf(fid, 'exit;"\n');

            fclose(fid);

            % Make script executable
            system(sprintf('chmod +x %s', slurm_script_path));

            fprintf('  ✓ SLURM script created: %s\n', slurm_script_name);

            slurm_scripts{end+1} = slurm_script_path;

        otherwise
            error('Invalid MODE: %s', MODE);
    end

end

%% =========================== SUMMARY ====================================

fprintf('\n========================================\n');
fprintf('PROCESSING COMPLETE\n');
fprintf('========================================\n');

switch MODE
    case 'local'
        if ~isempty(processLog)
            processTable = cell2table(processLog, ...
                'VariableNames', {'Filename', 'Status'});
            disp(processTable);
        end

        fprintf('\nNext: Run RS6_IClabel.m\n');

    case 'generate_slurm'
        fprintf('\nGenerated %d SLURM job scripts\n', length(slurm_scripts));
        fprintf('Location: %s\n\n', SLURM_OUTPUT_DIR);

        fprintf('To submit jobs:\n');
        fprintf('  cd %s\n', SLURM_OUTPUT_DIR);
        fprintf('  for script in *.sh; do sbatch "$script"; done\n\n');

        fprintf('To monitor jobs:\n');
        fprintf('  squeue -u $USER\n\n');

        % Create master submission script
        master_script = fullfile(SLURM_OUTPUT_DIR, 'submit_all_amica_jobs.sh');
        fid = fopen(master_script, 'w');
        fprintf(fid, '#!/bin/bash\n');
        fprintf(fid, '# Submit all AMICA jobs\n\n');
        for i = 1:length(slurm_scripts)
            fprintf(fid, 'sbatch %s\n', slurm_scripts{i});
        end
        fclose(fid);
        system(sprintf('chmod +x %s', master_script));

        fprintf('Or run master script:\n');
        fprintf('  bash %s\n\n', master_script);

        fprintf('After all jobs complete, run RS6_IClabel.m\n');
end

fprintf('========================================\n\n');
