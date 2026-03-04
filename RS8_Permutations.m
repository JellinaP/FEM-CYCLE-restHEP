%% Initiate
clc; clear

addpath /data/u_prinsen_software/eeglab2024.2.1
addpath /data/u_prinsen_software/fieldtrip-20250507
ft_defaults

basePath = '/data/pt_02938/rsHEP_analysis/';
cd(basePath)

subDirs = dir(fullfile(basePath, '/0_data/', 'NC*'));

%% LOAD ALL DATA AND CONVERT TO FIELDTRIP DATASET STRUCTURE
eeglab nogui
ft_HEPdata = struct();
ft_subjMeta = {};

icaType = 'iclabel'; % icCFA or iclabel
phaseTable = readtable('/data/p_02938/MINT_EEG/participant_cycle_phases.xlsx');

for i = 1:length(subDirs)
    sessions = dir(fullfile(subDirs(i).folder, subDirs(i).name, 'ses*'));
    
    if numel(sessions) == 2
        for s = 1:2
            eegFile = dir(fullfile(subDirs(i).folder, subDirs(i).name, sessions(s).name, '*_HEP_icCFA.set'));
            
            if isempty(eegFile)
                fprintf('No HEP file in %s/%s. Skipping session.\n', eegFile.name, sessions(s).name);
                continue;
            end

            EEG = pop_loadset(eegFile.name, eegFile.folder);
            subjID = subDirs(i).name;

            matchIdx = strcmp(phaseTable.Subject, subjID) & phaseTable.Session == s;
            if any(matchIdx)
                Phase = phaseTable.Phase{matchIdx};
            else
                warning('No phase mapping for %s session %d. Skipping.', subjID, sesNum);
                continue;
            end

            ft_data = eeglab2fieldtrip(EEG,'timelockanalysis','none');

            if ~isfield(ft_HEPdata,Phase)
                ft_HEPdata.(Phase) = {};
            end

            ft_HEPdata.(Phase){end+1} = ft_data;
            ft_subjMeta{end+1} = subjID;

        end
    else
        fprintf('Skipping %s: found %d session(s)\n', subDirs(i).name, numel(sessions));
    end
end

save(fullfile(basePath, '04_results', ['HEPs_' icaType '.mat']), 'ft_HEPdata', 'ft_subjMeta'); 
EEG = pop_saveset(EEG, '64_channels_example', fullfile(basePath, '04_results/01_Phase_main/icCFA'));

%% CLUSTER-BASED PERMUTATIONS IN FIELDTRIP
% Prepare neighbours
elec = ft_read_sens(fullfile(basePath, '64_channels_example.set'), 'senstype', 'eeg');

cfg_neighb          = [];
cfg_neighb.method   = 'triangulation'; % or distance
cfg_neighb.compress = 'yes';
cfg_neighb.feedback = 'no';
cfg_neighb.elec     = elec;

neighbours = ft_prepare_neighbours(cfg_neighb);

% Parameters for ft_timelockstatistics
cfg         = [];
cfg.channel = {'EEG'};
cfg.latency = [0.2 0.6];

cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.neighbours       = neighbours;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.computecritval   = 'no';
cfg.numrandomization = 4000;

% Define design matrix
nSubj = numel(ft_HEPdata.follicular);  

design = zeros(2, nSubj*2);
design(1,:) = [1:nSubj 1:nSubj]; % subject index * 2
design(2,:) = [ones(1,nSubj) ones(1,nSubj)*2]; % condition

% Assign design and variable roles
cfg.design = design;
cfg.uvar   = 1;  % unit of observation (subject)
cfg.ivar   = 2;  % independent variable (phase)

timelockFOL = ft_HEPdata.follicular;
timelockLUT = ft_HEPdata.luteal;

% Run permutations
stat = ft_timelockstatistics(cfg, timelockFOL{:}, timelockLUT{:});

% Calculate the average timecourse of the HEP per group
cfg = [];
cfg.keepindividual = 'no';
avgFOL = ft_timelockgrandaverage(cfg, timelockFOL{:});
avgLUT = ft_timelockgrandaverage(cfg, timelockLUT{:});

% Take the raw difference of the timecourses
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
diff = ft_math(cfg, avgFOL, avgLUT);

% Save
save(fullfile(basePath, '04_results', 'timelock_stats.mat'), 'avgFOL', 'avgLUT', 'diff', 'stat');